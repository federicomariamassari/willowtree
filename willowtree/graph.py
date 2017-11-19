def graph(z, q, gamma, t, P):
    '''
    Plot the willow tree.

    Input
    ---------------------------------------------------------------------------
    z, q: np.arrays, required arguments. The discrete density pairs, output of
          the function 'sampling'.
    gamma: float, required argument. Weighting parameter governing the shape of
           the distribution of probabilities {q(i)} in 'sampling'.
    t: np.array, required argument. The time array, a partition of [0, T].
    P: np.array, 2-D or 3-D, required argument. The Markov chain, output of the
       function 'lp'.
    '''

    # Import required libraries
    import numpy as np
    import matplotlib.pyplot as plt
    import seaborn as sns

    '''
    Define auxiliary functions to automate repetitive tasks like transposing,
    reshaping, and computing Kronecker product.
    '''
    def aux1(a, n, t):
        res = a.transpose().reshape(n ** 2 * (len(t) - 2))
        return res

    def aux2(a, b):
        res = np.kron(a, np.ones(b))
        return res

    # Make bidimensional grid G, the structure of the willow tree
    G = z[:, np.newaxis] @ np.sqrt(t)[np.newaxis]

    # Assign the space and time dimensions to variables n, k
    n, k = len(z), len(t)-1

    '''
    If at least one transition matrix is well-defined, also use P(k) to plot
    the willow tree.
    '''
    if k > 1:
        '''
        Define the initial ramification, from t(0) to t(1). 'initial' stacks,
        in order: the start date t(0); the end date t(1); the initial node
        z(0); the end nodes z(1, i), with 1 denoting t(1) and i the node in
        1, ..., n; the transition probabilities q(i) associated to each path
        z(0) -> z(1, i).
        '''
        initial = np.array([np.zeros(n), np.full(n, t[1]), G[:,0], G[:,1], q])

        '''
        Define additional ramifications using P(k) to plot, for each node, the
        n**2 transition probabilities.
        '''
        start_date = aux1(aux2(t[1:-1], (n ** 2)), n, t)
        end_date = aux1(aux2(t[2:], (n ** 2)), n, t)
        start_node = aux1(aux2(G[:,1:-1], (n, 1)), n, t)
        end_node = aux1(aux2(G[:,2:], n), n, t)
        transition = aux1(np.hstack(P), n, t)

        # Stack the arrays in a single structure W
        W = np.vstack((start_date, end_date, start_node, end_node, transition))
        W = np.hstack((initial, W)).transpose()

    else:
        '''
        If no well-defined transition matrix could be generated (k = 1), only
        plot the initial ramification.
        '''
        initial = np.array([np.zeros(n), np.full(n, t[1]), G[:,0], G[:,1], q])
        W = initial.transpose()

    '''
    Include a square root function to define the upper and lower limits of the
    tree.
    '''
    steps = np.linspace(0, t[-1], 1000)
    square_root = np.sqrt(steps) * z[n - 1]

    # Plot grid G together with the square root bounds
    sns.set()
    plt.figure(figsize = (10, 8))
    ax = plt.axes()
    ax.plot(t, G.transpose(), '.k')
    l1 = ax.plot(steps, square_root, 'k', steps, -square_root, 'k')
    l2 = [ax.plot(W[i,:2], W[i,2:4], 'k') for i in range(len(W)) if W[i,4] > 0]
    plt.setp((l1, l2), linewidth = 0.5)
    ax.set(title = \
    'Willow Tree, $n$ = {} space points, $k$ = {} time steps, $\gamma$ = {}'\
        .format(n, k, gamma), xlabel = 'Time', ylabel = 'Space')
    ax.invert_yaxis()
    plt.show()
