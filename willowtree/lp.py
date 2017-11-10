def lp(z, q, k, tol = 1e-9, extra_precision = False):
    '''
    Generate a time-inhomogeneous, discrete time Markov chain for the willow
    tree [1] via linear programming (LP), using the discrete density pairs
    {z(i), q(i)}, for i = 1, ..., n, output of the function 'sampling'.

    The willow tree linear programming problem is:

                                      c'x
                         subject to:
                                      A_eq * x = b_eq
                                      p(k; i,j) >= 0

    with:

     * c: the vector of coefficients of the linear objective function;
     * A_eq: a matrix of linear equality constraints;
     * b_eq: a vector of linear equality constraints;
     * p(k; i,j): transition probability at position (i,j) in the k-th
                  transition matrix;
     * x: the array solution to the problem.

    Each solution x, when reshaped, is a transition matrix P(k), for k = 1,
    ..., n-1.

    Input
    ---------------------------------------------------------------------------
    q, z: NumPy arrays, required arguments. The discrete density pairs, a
          discrete approximation of the standard normal distribution. Output
          of the function 'sampling';
    k: int, required argument. The number of time steps. k-1 is the number of
       transition matrices generated;
    tol: float, generally in scientific notation, optional argument. Set the
         precision of the solutions to the linear programming problems.
    extra_precision: bool, optional argument. If True, set the upper bound of
                     each variable p(i,j) in the LP problems to 1. Otherwise,
                     leave it to None.

    Output
    ---------------------------------------------------------------------------
    P: NumPy array. The Markov chain, whose elements are transition matrices.
       P is 2-dim if either of the following is true: k = 2, len(t_new) = 3.
       Otherwise, P is a 3-dim array with shape (k-1, len(z), len(z));
    t_new: Numpy array.

    How does the algorithm work?
    ---------------------------------------------------------------------------


    Resources
    ---------------------------------------------------------------------------
    [1] Curran, M. (2001). Willow Power: Optimizing Derivative Pricing Trees,
        ALGO Research Quarterly, Vol. 4, No. 4, p. 15, December 2001.
    '''

    # Import required libraries
    import time
    import numpy as np
    from scipy import optimize

    def objective(z, a, beta, normalize):
        F = (np.abs(a-beta*a.transpose()) ** 3).transpose()\
            .reshape(len(z)**2)
        c = F * normalize
        return c

    def beq(q, u, z, beta, Aeq):
        beq = np.array([u, beta*z, (beta**2)*r + (1-beta**2)*u,
                        q]).reshape(len(Aeq))
        return beq

    def transition_matrix(z, c, Aeq, beq, tol, extra_precision):
        if extra_precision:
            bounds = (0, 1)
        else:
            bounds = (0, None)

        options = {'maxiter': 1e4, 'tol': tol, 'disp': False}
        P = optimize.linprog(c, A_eq = Aeq, b_eq = beq,
                             bounds = bounds, method = 'simplex',
                             options = options)
        return P

    def test(n, P):
        try:
            P = P.reshape(n,n)
            return np.isclose(P.sum(axis=1), np.ones(n), 1e-6).all() == True
        except:
            return False

    def interpolate(P_min, P_max, alpha_min, alpha_max, alpha_interp):
        x1 = 1 / np.sqrt(1+alpha_min)
        x2 = 1 / np.sqrt(1+alpha_max)
        x3 = 1 / np.sqrt(1+alpha_interp)

        coeff_min = (x3-x2) / (x1-x2)
        coeff_max = (x1-x3) / (x1-x2)

        return coeff_min*P_min + coeff_max*P_max


    initial_tol = tol
    n = len(z)
    t = np.linspace(0, 1, k + 1)

    u = np.ones(len(z), dtype = np.int)
    r = z ** 2
    h = t[2:] - t[1:-1]
    alpha = h / t[1:-1]
    beta = 1 / np.sqrt(1+alpha)

    a = z[:, np.newaxis] @ np.ones(len(z))[np.newaxis]
    normalize = np.kron(q, np.ones(len(z)))

    c = np.array([objective(z, a, beta[i], normalize) \
                  for i in range(len(h))])

    Aeq = np.vstack([np.kron(np.eye(len(z)), u),
                     np.kron(np.eye(len(z)), z),
                     np.kron(np.eye(len(z)), r),
                     np.kron(q, np.eye(len(z)))])

    beq = np.array([beq(q, u, z, beta[i], Aeq) \
                    for i in range(len(h))])

    Px = np.array([np.zeros([len(z), len(z)]) \
         for i in range(len(h))])

    flag = np.zeros(len(h), dtype = np.int)

    for i in range(len(h)):
        P = transition_matrix(z, c[i], Aeq, beq[i], tol,
                              extra_precision)

        if type(P.x) != np.float:
            start = time.time()
            while (P.status != 0) | (P.fun < 0) | (P.fun > 1) \
                | ((P.x[P.x < 0]).any()) | ((P.x[P.x > 1]).any()) \
                | (test(n, P.x) != True):

                if (tol < 1e-3) & (time.time() - start < 60):
                    tol *= 10
                    P = transition_matrix(z, c[i], Aeq, beq[i], tol,
                                          extra_precision)
                else:
                    flag[i] = -1
                    break

            Px[i] = P.x.reshape(len(z), len(z))

        else:
            flag[i] = -1


        if flag[i] == -1:
            print('Warning: P[{}] wrongly specified.'.format(i))
            print('Replacing with interpolated matrix if possible.')
        else:
            print('P[{}] successfully generated.'.format(i))

        tol = initial_tol

        failure = np.nonzero(flag)[0]
        success = np.nonzero(flag + 1)[0]

    try:
        minvec = np.array([], dtype = np.int)
        maxvec = minvec

        for i in reversed(range(len(failure))):
            minvec = np.append(minvec, [max(x for x in success if x < failure[i])])
            minvec.sort()
    except ValueError:
        pass

    try:
        for i in range(len(failure)):
            maxvec = np.append(maxvec, [min(x for x in success if x > failure[i])])
    except ValueError:
        pass

    repl_min = np.full(len(flag), -1, dtype=np.int)
    repl_max = np.full(len(flag), -1, dtype=np.int)

    for i in failure:
        repl_min[i] = 1
        repl_max[i] = 1

    minvec = np.pad(minvec, ((len(failure)-len(minvec)),0),
                    mode='constant', constant_values=-1)
    repl_min[repl_min>0] = minvec

    maxvec = np.pad(maxvec, (0,len(failure) - len(maxvec)),
                    mode='constant', constant_values=-1)
    repl_max[repl_max>0] = maxvec

    succ_vec = (repl_min > -1) & (repl_min < repl_max)
    succ_vec = np.array([1 if succ_vec[i] == True else 0 for i \
                         in range(len(succ_vec))])

    try:
        threshold_low = np.argwhere(succ_vec)[0,0]
        threshold_high = np.argwhere(succ_vec)[-1,0]
    except:
        threshold_low, threshold_high = -1, -1

    failure = failure[(failure >= threshold_low) \
                    & (failure <= threshold_high)]

    minvec = repl_min * succ_vec
    maxvec = repl_max * succ_vec

    minvec = minvec[minvec > -1]
    maxvec = maxvec[maxvec > 0]

    if (flag == -1).any():
        try:
            Px[failure] = [interpolate(Px[minvec[i]], Px[maxvec[i]],
                           alpha[minvec[i]], alpha[maxvec[i]],
                           alpha[failure[i]]) for i \
                           in range(len(failure))]
        except ValueError:
            pass

        for i in failure:
            print('Interpolation of P[{}] successful.'.format(i))
        flag[failure] = 0
    else:
        pass

    success = np.nonzero(flag + 1)[0]
    Px = Px[success]

    try:
        if success[0] == 0:
            t_new = t[range(len(success)+2)]
        else:
            t_new = np.append(0,t[(t >= t[success[0]+1]) \
                                & (t <= t[success[-1]+2])])

        if t_new[1] != t[1]:
            print('Warning: t has been increased. t[1] = {:.2f}'\
                  .format(t_new[1]))
        if t_new[-1] != t[-1]:
            print('Warning: t has been shortened. T = {:.2f}'\
                  .format(t_new[-1]))
    except:
        t_new = t[:2]
        print('Warning: t has been shortened. T = {:.2f}'.format(t_new[-1]))

    return Px, t_new
