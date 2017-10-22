def maketree(n = 12, gamma = 0.1, algorithm = 'kurtosis-matching', k = 10,
             tol = 1e-12, extra_precision = False):

    import time
    import numpy as np
    from scipy import stats, optimize
    import matplotlib.pyplot as plt
    import seaborn as sns

    q, z, gamma = sampling(n, gamma)
    P, t = lp(z, q, k, tol = tol, extra_precision = extra_precision)
    graph(z, q, gamma, t, P)

    return q, z, P, t

def sampling(n, gamma, algorithm = 'kurtosis-matching'):
    '''
    Generate a sequence of discrete density pairs {z(i), q(i)}, for i = 1, ...,
    n, according to Xu, Hong, and Qin's methodology [1].

    Choose between kurtosis matching ('kurtosis-matching', 'KRT', 'krt') and
    first partial moment matching ('first-partial-moment', 'FPM', 'fpm')
    strategies. Additional algorithms will be provided in future versions.

    Which algorithm should you choose?
    ---------------------------------------------------------------------------
     - 'kurtosis-matching': best suited for deep out-of-the-money derivatives
       since it proxies for kurtosis (the parameter governing the shape of the
       tails of the density);

     - 'first-partial-moment': best suited for very near-the-money derivatives
       since it provides extra precision in computing mean and variance (while
       relaxing the condition on kurtosis).

    Input
    ---------------------------------------------------------------------------
    n: int, required argument. The number of points in the space dimension.
    gamma: float, required argument. Weighting parameter governing the shape
           of the distribution of probabilities {q(i)}.
           Admissible values are in the closed interval [0;1]. Otherwise, the
           algorithm will set gamma to the nearest integer, e.g., if the user
           inputs gamma < 0, the parameter will adjust to 0; if gamma > 1, it
           will adjust to 1.
    algorithm: str, optional argument. The algorithm used to compute the pairs.
               For kurtosis matching omit or set to 'kurtosis-matching', 'KRT',
               or 'krt'. For first partial moment set to 'first-partial-moment',
               'FPM', or 'fpm'.

    Output
    ---------------------------------------------------------------------------
    q, z: NumPy arrays, the discrete density pairs, a discrete approximation of
          the standard normal distribution.
        * z is a sequence of representative standard normal variates determined
          with stratified sampling: partition the standard normal CDF into
          subgroups (strata) and draw a single value from each of them;
        * q is the set of probabilities associated to z and governing the width
          of the strata.

    How does the algorithm work?
    ---------------------------------------------------------------------------
    Suppose you want to obtain {z(i), q(i)}, for i = 1, ..., 14, a 14-element
    partition of the standard normal CDF. This means the willow tree will have
    14 (vertical) space points.

    * You start with n = 14, and decide to choose gamma = 0. However, with such
      value of gamma there exists no solution. The algorithm will increase the
      level of gamma by 1e-9 for the first two seconds of runtime, then by 1e-6
      for the following 8 seconds (up to a total of ten seconds), and finally
      by 1e-2 until a solution is found;
    * Suppose, instead, you start with gamma = 1, and also no solution exists
      for such value of gamma. The algorithm will decrease gamma by 1e-9 for
      the first two seconds, then by 1e-6 for additional 8 seconds, and finally
      by 1e-2 until a solution is found;
    * What if, for n = 14, no solution exists for gamma in [0;1]? In this case,
      after a full gamma cycle, n will increase by 1 and gamma be reset to the
      original value input by the user. Another full gamma cycle will occur and
      if again no solution is found n is increased by 1 additional unit, and so
      on.

    Example
    ---------------------------------------------------------------------------
    q, z = sampling(n = 14, gamma = 0, algorithm = 'first-partial-moment')

    Compute the discrete density pairs with n = 14 space points, gamma = 0.
    With gamma 0 the problem is infeasible. An optimal solution is found with
    n = 14 and gamma = 0.000000032 (increased by 1e-9 32 times).

    Resources
    ---------------------------------------------------------------------------
    [1] Xu, W., Hong, Z. and Qin, C. (2013): A New Sampling Strategy Willow
    Tree Method With Application To Path-Dependent Option Pricing. Quantitative
    Finance, March 2013, 13(6): 861–872.
    '''

    # Define auxiliary functions to compute the discrete density pairs
    def aux(n, gamma):
        '''
        Generate q, the array of probabilities governing the width of the
        strata, and Z, the strata themselves.
        '''

        '''
        Generate index i, an array of linearly spaced natural numbers,
        accounting for n odd, if necessary. Compute only half, then flip and
        add the other half (plus a midpoint, if required).
        '''
        i = np.arange(1, 0.5*(n+1), dtype = int)
        i = np.hstack([i, i[::-1]] if n % 2 == 0 \
                       else [i, i[-1]+1, i[::-1]])

        '''
        Compute the array of probabilities q as a normalised, either linear (if
        gamma ≈ 0) or nonlinear transformation of i.
        '''
        q = (i-0.5)**gamma / n
        q = q / np.sum(q)

        '''
        Compute Z, the array of bounds for the standard normal variates {z(i)}
        (the endpoints of the strata). Set Z[0] and Z[-1] to, respectively,
        -Inf and Inf.
        '''
        Z = np.hstack([[-np.inf], stats.norm.ppf(np.cumsum(q[:-1])), \
                       [np.inf]])

        return q, Z

    def aux2(q, Z, algorithm):
        '''
        Compute z, the sequence of representative standard normal variates.
        '''

        '''
        Choose which algorithm to use based on user input. If the input is not
        recognised, switch to 'kurtosis-matching'.
        '''
        if algorithm in ('kurtosis-matching', 'KRT', 'krt'):
            fun = lambda z: (q @ (z**4) - 3) ** 2
        elif algorithm in ('first-partial-moment', 'FPM', 'fpm'):
            fun = lambda z: np.sum(np.abs(q @ (z[:, np.newaxis] \
                - Z[1:-1][np.newaxis, :]).clip(np.min(0)) \
                - (1/np.sqrt(2*np.pi) * np.exp(0.5 * (-Z[1:-1]**2)) \
                - Z[1:-1] * (1-stats.norm.cdf(Z[1:-1])))))
        else:
            fun = lambda z: (q @ (z**4) - 3) ** 2

        '''
        Set parameters and options for the optimisation algorithm.
         - x0: array, the vector of initial guess points for z;
         - constraints: dict, the set of constraints on q and z; (1) expected
           value equal to 0; (2) variance of z equal to 1; (3) sum of variates
           equal to zero (i.e. z is symmetric); (4) equal endpoints;
         - bounds: array of tuples, the bounds for each z, Z[low] and Z[high];
         - options: additional options, such as the solution tolerance or the
                    maximum number of iterations.
        '''
        x0 = np.full(n, 1e-6)
        constraints = ({'type': 'eq', 'fun': lambda z: q @ z},
                       {'type': 'eq', 'fun': lambda z: q @ (z**2) - 1},
                       {'type': 'eq', 'fun': lambda z: np.sum(z)},
                       {'type': 'eq', 'fun': lambda z: z[0] + z[-1]})
        bounds = np.column_stack((Z[:-1], Z[1:]))
        options = {'disp': False, 'ftol': 1e-15, 'maxiter': 1e4}

        '''
        Optimisation procedure.
        '''
        z = optimize.minimize(fun, x0, bounds = bounds, options = options,
                              tol = 1e-15, constraints = constraints)

        return z

    # Define test of consistency for the solution array z
    def test(q, z, algorithm):
        '''
        Test whether q, z are an admissible set of discrete density pairs,
        based on the chosen algorithm. For both algorithms, q, z must have
        mean = 0, and variance = 1 within a strict tolerance level (1e-15).
        For 'kurtosis-matching' only, q and z must also satisfy kurtosis = 3.
        Otherwise, the algorithm is run again, until all requirements are met.
        '''
        mean = np.isclose(q @ z, 0, 1e-15)
        variance = np.isclose(q @ z ** 2, 1, 1e-15)
        if algorithm in ('kurtosis-matching', 'KRT', 'krt'):
            kurtosis = np.isclose(q @ z ** 4, 3, 1e-15)
            return np.array([mean, variance, kurtosis]).all()
        else:
            return np.array([mean, variance]).all()

    # Time the entire algorithm
    tic = time.time()

    # Store the initial values of parameters n and gamma
    initial_n = n
    initial_gamma = gamma

    '''
    Check whether the input algorithm is valid, otherwise switch to kurtosis-
    matching.
    '''
    if algorithm not in ('kurtosis-matching', 'KRT', 'krt',
                         'first-partial-moment', 'FPM', 'fpm'):
        print("Unrecognised algorithm. Switching to 'kurtosis-matching'.")

    '''
    If gamma is within the admissible range, generate q and Z. Else, set gamma
    equal to the integer closer to the user input (either 0 if gamma < 0 or 1
    if gamma < 1), then generate q and Z.
    '''
    if 0 <= gamma <= 1:
        q, Z = aux(n, gamma)
    elif gamma < 0:
        print('Gamma out of range, must be in [0, 1]. Setting to 0.')
        gamma = 0
        q, Z = aux(n, gamma)
    else:
        print('Gamma out of range, must be in [0, 1]. Setting to 1.')
        gamma = 1
        q, Z = aux(n, gamma)

    '''
    Generate z, the array of standard normal variates, from q and Z, then test
    whether q, z are an admissible solution pair.
    z, the output of the optimisation algorithm, is a structure:
     - z.x is the solution vector (the array of standard normal variates);
     - z.fun is the value of the objective function, 'fun';
     - z.status is the output flag: 0 if successful, != 0 otherwise.
    '''
    z = aux2(q, Z, algorithm)
    test_result = test(q, z.x, algorithm)

    '''
    Set the timer and a counter for additional minimisation problems.
    - Timer: governs further variations in gamma. Plus (if initial gamma = 0)
             or minus (if initial gamma = 1) 1e-9 for the first two seconds;
             plus or minus 1e-6 for additional 8 seconds, and plus or minus
             1e-2 afterwards.
    - Counter: governs further increases in n. After a full cycle of gamma,
               either (0 -> 1) or (1 -> 0), the counter is increased by one and
               n is also raised by one. After that, the counter is reset to 0.
    '''
    start = time.time()
    counter = 0

    '''
    While the solution of the minimisation problem, z.x, is not satisfactory,
    either z.status != 0, z.fun > 0, or the constraints on mean, variance, and
    kurtosis are not strictly respected, the problem is run again.
    First, relaxing the constraints on gamma; then, increasing n by one or
    more units, if required.
    '''
    while (z.status != 0) | (z.fun > 1) | (test_result != True):

        # If the input value of gamma is smaller than one:
        if initial_gamma < 1:

            # Increase gamma while it is still smaller than 1 and counter = 0
            if (gamma < 1) & (counter <= 1):
                if time.time() < start + 2:
                    gamma += 1e-9
                elif time.time() < start + 10:
                    gamma += 1e-6
                else:
                    gamma += 1e-2
                q, Z = aux(n, gamma)
                z = aux2(q, Z, algorithm)
                test_result = test(q, z.x, algorithm)

            # Raise counter by 1 and reset gamma to 0 if gamma reaches 1
            elif (gamma == 1) & (counter <= 1):
                gamma = 0
                counter +=1
                q, Z = aux(n, gamma)
                z = aux2(q, Z, algorithm)
                test_result = test(q, z.x, algorithm)

            # Increase n by one unit after a full increasing cycle of gamma
            else:
                gamma = initial_gamma
                counter = 0
                n += 1
                q, Z = aux(n, gamma)
                z = aux2(q, Z, algorithm)
                test_result = test(q, z.x, algorithm)

        # If the input value of gamma is equal to one:
        else:

            # Decrease gamma while it is still higher than 0 and counter = 0
            if (gamma > 0) & (counter <= 1):
                if time.time() < start + 2:
                    gamma -= 1e-9
                elif time.time() < start + 10:
                    gamma -= 1e-6
                else:
                    gamma -= 1e-2
                q, Z = aux(n, gamma)
                z = aux2(q, Z, algorithm)
                test_result = test(q, z.x, algorithm)

            # Raise counter by 1 and reset gamma to 1 if gamma reaches 0
            elif (gamma == 0) & (counter <= 1):
                gamma = 1
                counter +=1
                q, Z = aux(n, gamma)
                z = aux2(q, Z, algorithm)
                test_result = test(q, z.x, algorithm)

            # Increase n by one unit after a full decreasing cycle of gamma
            else:
                gamma = initial_gamma
                counter = 0
                start = time.time()
                n += 1
                q, Z = aux(n, gamma)
                z = aux2(q, Z, algorithm)
                test_result = test(q, z.x, algorithm)

    '''
    If the values of n, gamma, or both, are different from those supplied by
    the user, print warning.
    '''
    if (n != initial_n) | (gamma != initial_gamma):
        print('Optimisation could not be terminated with supplied parameters.')
        print('Optimisation successful with n = {} and gamma = {:.9f}.' \
          .format(n, gamma))
    else:
        print('Optimisation terminated successfully with supplied parameters.')

    # Compute total running time, in seconds
    toc = time.time() - tic
    print('Total running time: {:.2f}s'.format(toc))

    return q, z.x, gamma


def lp(z, q, k, tol = 1e-12, extra_precision = False):

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

    def interpolate(P_min, P_max, alpha_min, alpha_max, alpha_interp):
        x1 = 1 / np.sqrt(1+alpha_min)
        x2 = 1 / np.sqrt(1+alpha_max)
        x3 = 1 / np.sqrt(1+alpha_interp)

        coeff_min = (x3-x2) / (x1-x2)
        coeff_max = (x1-x3) / (x1-x2)

        return coeff_min*P_min + coeff_max*P_max


    initial_tol = tol
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
                | ((P.x[P.x < 0]).any()) | ((P.x[P.x > 1]).any()):

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

    if success[0] == 0:
        t_new = t[range(len(success)+2)]
    else:
        t_new = np.append(0,t[(t >= t[success[1]]) \
                            & (t <= t[success[-1]+2])])

    if t_new[1] != t[1]:
        print('Warning: t has been increased. t[1] = {:.2f}'\
              .format(t_new[1]))
    if t_new[-1] != t[-1]:
        print('Warning: t has been shortened. T = {:.2f}'\
              .format(t_new[-1]))

    return Px, t_new

def graph(z, q, gamma, t, P):

    def aux1(a, n, t):
        res = a.transpose().reshape(n ** 2 * (len(t) - 2))
        return res

    def aux2(a, b):
        res = np.kron(a, np.ones(b))
        return res

    G = z[:, np.newaxis] @ np.sqrt(t)[np.newaxis]

    n, k = len(z), len(t)


    initial = np.array([np.zeros(n), np.full(n, t[1]), G[:,0], G[:,1], q])
    start_date = aux1(aux2(t[1:-1], (n ** 2)), n, t)
    end_date = aux1(aux2(t[2:], (n ** 2)), n, t)
    start_node = aux1(aux2(G[:,1:-1], (n, 1)), n, t)
    end_node = aux1(aux2(G[:,2:], n), n, t)
    transition = aux1(np.hstack(P), n, t)

    W = np.vstack((start_date, end_date, start_node, end_node, transition))
    W = np.hstack((initial, W)).transpose()

    steps = np.linspace(0, t[-1], 1000)
    square_root = np.sqrt(steps) * z[n - 1]

    sns.set()
    plt.figure(figsize = (10, 8))
    ax = plt.axes()
    ax.plot(t, G.transpose(), '.k')
    l1 = ax.plot(steps, square_root, 'k', steps, -square_root, 'k')
    l2 = [ax.plot(W[i,:2], W[i,2:4], 'k') for i in range(len(W)) if W[i,4] > 0]
    plt.setp((l1, l2), linewidth = 0.5)
    ax.set(title = \
    'Willow Tree, $n$ = {} space points, $k$ = {} time nodes, $\gamma$ = {}'\
        .format(n, k - 1, gamma), xlabel = 'Time', ylabel = 'Space')
    ax.invert_yaxis()
    plt.show()
