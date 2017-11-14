def maketree(n=12, gamma=0.1, algorithm='kurtosis-matching', k=10,
             tol=1e-12, extra_precision=False):
    '''
    Generate and plot the willow tree in a single step.

    Input
    ---------------------------------------------------------------------------
    n: int, optional argument. The number of space nodes in the lattice;
    gamma: float, optional argument. Weighting parameter governing the shape
           of the distribution of probabilities {q(i)} in 'sampling'.
           Admissible values in [0, 1]. If input gamma < 0, the parameter will
           adjust to 0; if input gamma > 1, it will adjust to 1;
    algorithm: str, optional argument. The algorithm used to generate the
               discrete density pairs in 'sampling'.
                * For kurtosis matching, omit or set to 'kurtosis-matching',
                  'KRT', or 'krt';
                * For first partial moment set to 'first-partial-moment',
                  'FPM', or 'fpm'.
               Which algorithm should you choose? 'kurtosis-matching' is best
               suited for deep out-of-the-money derivatives, since it proxies
               for kurtosis (the parameter governing the shape of the tails of
               the density); 'first-partial-moment' is instead best suited for
               very near-the-money derivatives, since it provides extra
               precision in computing mean and variance (while relaxing the
               condition on kurtosis);
    k: int, optional argument. The number of time steps, with k-1 the number
       of transition matrices generated;
    tol: float, generally in scientific notation, optional argument. Set the
         precision of the solutions to the linear programming problems.
    extra_precision: bool, optional argument. If True, also specify the upper
                     bound for each variable p(i,j) in each of the k-1 linear
                     programming (LP) problems for the Markov chain. Otherwise
                     leave None as upper bound (extra_precision=False).

    Output
    ---------------------------------------------------------------------------
    q, z: NumPy arrays, the discrete density pairs, a discrete approximation of
          the standard normal distribution.
        * z is a sequence of representative standard normal variates determined
          with stratified sampling: partition the standard normal CDF into
          subgroups (strata) and draw a single value from each of them;
        * q is the set of probabilities associated to z and governing the width
          of the strata.
    P: NumPy array. The Markov chain, whose elements are transition matrices.
       P is 2-dim if either of the following is true: k = 2, len(t_new) = 3.
       Otherwise, P is a 3-dim array with shape (k-1, len(z), len(z));
    t: Numpy array. The vector of time nodes. Of length k+1 if the algorithm
       manages to generate the full Markov chain (including both well-defined
       and interpolated matrices); shorter, otherwise.

    Example
    ---------------------------------------------------------------------------
    q, z, P, t = maketree()
    Generate willow tree with the suggested parameters.
    '''

    # Import required libraries
    import time
    import numpy as np
    from scipy import stats, optimize
    import matplotlib.pyplot as plt
    import seaborn as sns

    # Import companion willowtree modules
    from willowtree import sampling, lp, graph

    # Generate discrete density pairs
    q, z, gamma = sampling(n, gamma)

    # Generate Markov chain
    P, t = lp(z, q, k, tol = tol, extra_precision = extra_precision)

    # Plot results
    graph(z, q, gamma, t, P)

    return q, z, P, t
