def graph(z, q, gamma, t, P):

    import numpy as np
    import matplotlib.pyplot as plt
    import seaborn as sns

    def aux1(a, n, t):
        res = a.transpose().reshape(n ** 2 * (len(t) - 2))
        return res

    def aux2(a, b):
        res = np.kron(a, np.ones(b))
        return res

    G = z[:, np.newaxis] @ np.sqrt(t)[np.newaxis]

    n, k = len(z), len(t)-1

    if k > 1:
        initial = np.array([np.zeros(n), np.full(n, t[1]), G[:,0], G[:,1], q])
        start_date = aux1(aux2(t[1:-1], (n ** 2)), n, t)
        end_date = aux1(aux2(t[2:], (n ** 2)), n, t)
        start_node = aux1(aux2(G[:,1:-1], (n, 1)), n, t)
        end_node = aux1(aux2(G[:,2:], n), n, t)
        transition = aux1(np.hstack(P), n, t)

        W = np.vstack((start_date, end_date, start_node, end_node, transition))
        W = np.hstack((initial, W)).transpose()

    else:
        initial = np.array([np.zeros(n), np.full(n, t[1]), G[:,0], G[:,1], q])
        W = initial.transpose()

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
    'Willow Tree, $n$ = {} space points, $k$ = {} time steps, $\gamma$ = {}'\
        .format(n, k, gamma), xlabel = 'Time', ylabel = 'Space')
    ax.invert_yaxis()
    plt.show()
