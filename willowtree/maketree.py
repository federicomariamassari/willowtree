def maketree(n = 12, gamma = 0.1, algorithm = 'kurtosis-matching', k = 10,
             tol = 1e-12, extra_precision = False):

    import time
    import numpy as np
    from scipy import stats, optimize
    import matplotlib.pyplot as plt
    import seaborn as sns

    from willowtree import sampling, lp, graph

    q, z, gamma = sampling(n, gamma)
    P, t = lp(z, q, k, tol = tol, extra_precision = extra_precision)
    graph(z, q, gamma, t, P)

    return q, z, P, t
