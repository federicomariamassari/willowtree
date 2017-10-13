def sampling(n, gamma, algorithm = 'kurtosis-matching'):
    
    import time
    import numpy as np
    from scipy import stats, optimize

    def aux(n, gamma):
        i = np.arange(1, 0.5*(n+1), dtype = int)
        i = np.hstack([i, i[::-1]] if n % 2 == 0 else [i, i[-1]+1, i[::-1]])

        q = (i-0.5)**gamma / n
        q = q / np.sum(q)

        Z = np.hstack([[-np.inf], stats.norm.ppf(np.cumsum(q[:-1])), [np.inf]])

        return q, Z

    def aux2(q, Z, algorithm):

        if algorithm in ('kurtosis-matching', 'KRT', 'krt'):
            fun = lambda z: (q @ (z**4) - 3) ** 2
        elif algorithm in ('first-partial-moment', 'FPM', 'fpm'):
            fun = lambda z: np.sum(np.abs(q @ (z[:, np.newaxis] \
                - Z[1:-1][np.newaxis, :]).clip(np.min(0)) \
                - (1/np.sqrt(2*np.pi) * np.exp(0.5 * (-Z[1:-1]**2)) \
                - Z[1:-1] * (1-stats.norm.cdf(Z[1:-1])))))
        else:
            fun = lambda z: (q @ (z**4) - 3) ** 2

        x0 = np.full(n, 1e-6)
        constraints = ({'type': 'eq', 'fun': lambda z: q @ z},
                       {'type': 'eq', 'fun': lambda z: q @ (z**2) - 1},
                       {'type': 'eq', 'fun': lambda z: np.sum(z)})
        bounds = np.column_stack((Z[:-1], Z[1:]))
        options = {'disp': False, 'ftol': 1e-15, 'maxiter': 1e4}

        z = optimize.minimize(fun, x0, bounds = bounds, options = options, \
                              tol = 1e-15, constraints = constraints)

        return z

    tic = time.time()

    initial_n = n
    initial_gamma = gamma
    if algorithm not in ('kurtosis-matching', 'KRT', 'krt',
                         'first-partial-moment', 'FPM', 'fpm'):
        print("Unrecognised algorithm. Switching to 'kurtosis-matching'.")

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

    z = aux2(q, Z, algorithm)

    start = time.time()
    counter = 0

    while (z.status != 0) | (z.fun > 1):
        if initial_gamma < 1:
            if (gamma < 1) & (counter <= 1):
                if time.time() < start + 2:
                    gamma += 1e-9
                elif time.time() < start + 10:
                    gamma += 1e-6
                else:
                    gamma += 1e-2
                q, Z = aux(n, gamma)
                z = aux2(q, Z, algorithm)
            elif (gamma == 1) & (counter <= 1):
                gamma = 0
                counter +=1
                q, Z = aux(n, gamma)
                z = aux2(q, Z, algorithm)
            else:
                gamma = initial_gamma
                counter = 0
                n += 1
                q, Z = aux(n, gamma)
                z = aux2(q, Z, algorithm)
        else:
            if (gamma > 0) & (counter <= 1):
                if time.time() < start + 2:
                    gamma -= 1e-9
                elif time.time() < start + 10:
                    gamma -= 1e-6
                else:
                    gamma -= 1e-2
                q, Z = aux(n, gamma)
                z = aux2(q, Z, algorithm)
            elif (gamma == 0) & (counter <= 1):
                gamma = 1
                counter +=1
                q, Z = aux(n, gamma)
                z = aux2(q, Z, algorithm)
            else:
                gamma = initial_gamma
                counter = 0
                start = time.time()
                n += 1
                q, Z = aux(n, gamma)
                z = aux2(q, Z, algorithm)

    if (n != initial_n) | (gamma != initial_gamma):
        print('Optimisation could not be terminated with supplied parameters.')
        print('Optimisation successful with n = {} and gamma = {:.9f}.' \
          .format(n, gamma))
    else:
        print('Optimisation terminated successfully with supplied parameters.')

    toc = time.time() - tic
    print('Total running time: {:.2f}s'.format(toc))

    return q, z.x
