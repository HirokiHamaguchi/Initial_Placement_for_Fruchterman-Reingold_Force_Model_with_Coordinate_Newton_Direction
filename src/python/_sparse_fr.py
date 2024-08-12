from matplotlib.pylab import f
import networkx as nx
from networkx.utils import np_random_state
from src.python.cost import cost
from tqdm.auto import tqdm
import time


# Copied from networkx.drawing.layout.py
@np_random_state(7)
def _sparse_fruchterman_reingold(
    A,
    k=None,
    pos=None,
    fixed=None,
    iterations=50,
    threshold=1e-4,
    dim=2,
    seed=None,
    method="FR",
    verbose=True,
):
    # Position nodes in adjacency matrix A using Fruchterman-Reingold
    # Entry point for NetworkX graph is fruchterman_reingold_layout()
    # Sparse version
    import numpy as np
    import scipy as sp
    from scipy.optimize import minimize

    try:
        nnodes, _ = A.shape
    except AttributeError as err:
        msg = "fruchterman_reingold() takes an adjacency matrix as input"
        raise nx.NetworkXError(msg) from err

    if pos is None:
        # random initial positions
        pos = np.asarray(seed.rand(nnodes, dim), dtype=A.dtype)
    else:
        # make sure positions are of same type as matrix
        pos = pos.astype(A.dtype)

    # no fixed nodes
    if fixed is None:
        fixed = []

    # optimal distance between nodes
    if k is None:
        k = np.sqrt(1.0 / nnodes)

    if method == "FR":
        # make sure we have a LIst of Lists representation
        try:
            A = A.tolil()
        except AttributeError:
            A = (sp.sparse.coo_array(A)).tolil()

        # the initial "temperature" is about .1 of domain area (=1x1)
        # this is the largest step allowed in the dynamics.
        t = max(max(pos.T[0]) - min(pos.T[0]), max(pos.T[1]) - min(pos.T[1])) * 0.1
        # simple cooling scheme.
        # linearly step down by dt on each iteration so last iteration is size dt.
        dt = t / (iterations + 1)

        displacement = np.zeros((dim, nnodes))
        for iteration in range(iterations):
            displacement *= 0
            # loop over rows
            for i in range(A.shape[0]):
                if i in fixed:
                    continue
                # difference between this row's node position and all others
                delta = (pos[i] - pos).T
                # distance between points
                distance = np.sqrt((delta**2).sum(axis=0))
                # enforce minimum distance of 0.01
                distance = np.where(distance < 0.01, 0.01, distance)
                # the adjacency matrix row
                Ai = A.getrowview(i).toarray()  # TODO: revisit w/ sparse 1D container
                # displacement "force"
                displacement[:, i] += (
                    delta * (k * k / distance**2 - Ai * distance / k)
                ).sum(axis=1)
            # update positions
            length = np.sqrt((displacement**2).sum(axis=0))
            length = np.where(length < 0.01, 0.1, length)
            delta_pos = (displacement * t / length).T
            pos += delta_pos
            # cool temperature
            t -= dt
            if verbose:
                print(f"{cost(pos,A,k)=}")
            if (np.linalg.norm(delta_pos) / nnodes) < threshold:
                break
            yield pos
    elif method == "BFGS":
        # make sure we have a coo_matrix representation
        try:
            A = A.tolil()
        except AttributeError:
            A = (sp.sparse.coo_array(A)).tolil()

        k_inv = 1 / k

        def cost_fun(x):
            pos = x.reshape((nnodes, dim))
            grad = np.zeros((nnodes, dim))
            for i in range(nnodes):
                delta = pos[i] - pos
                distance = np.linalg.norm(delta, axis=1)
                distance = np.where(distance < 0.01, 0.01, distance)
                distance_inv = 1 / distance
                Ai = A.getrow(i).toarray().flatten()
                coefficient1 = Ai * distance * k_inv - (k * distance_inv) ** 2
                grad[i] = coefficient1 @ delta
            ret = cost(pos, A, k)
            print(f"{ret=}")
            return ret, grad.ravel()

        pos_hist = []
        sp.optimize.minimize(
            cost_fun,
            pos.ravel(),
            method="L-BFGS-B",
            jac=True,
            options={"maxiter": iterations, "disp": verbose},
            callback=lambda x: pos_hist.append(x.reshape((nnodes, dim))),
        )

        for pos in pos_hist:
            yield pos

    elif method == "RS":
        # make sure we have a coo_matrix representation
        try:
            A = A.tolil()
        except AttributeError:
            A = (sp.sparse.coo_array(A)).tolil()

        vertices = np.delete(np.arange(nnodes), fixed)
        k_inv = 1 / k

        for it in range(iterations):
            # t = (it / iterations) ** 0.5
            t = 1
            np.random.shuffle(vertices)
            for i in vertices:

                def cost_fun_i(x):
                    # todo iの除外
                    assert t == 1, "change np.log for t!=1"
                    delta = x - pos
                    distance = np.linalg.norm(delta, axis=1)
                    distance = np.where(distance < 0.01, 0.01, distance)
                    distance_inv = 1 / distance
                    Ai = A.getrow(i).toarray().flatten()
                    coefficient1 = Ai * distance * k_inv - (k * distance_inv) ** 2
                    coefficient2 = Ai * distance_inv * k_inv + ((1 + t) * k**2) * (
                        distance_inv ** (3 + t)
                    )
                    cost = np.sum(
                        Ai * distance**3 * (1 / 3) * k_inv - (k**2) * np.log(distance)
                    )
                    grad = coefficient1 @ delta
                    hess = np.sum(coefficient1) * np.eye(dim) + np.einsum(
                        "ij,ik->jk", coefficient2[:, np.newaxis] * delta, delta
                    )

                    return cost, grad, hess

                pos_hist = []
                sp.optimize.minimize(
                    cost_fun_i,
                    pos.ravel(),
                    method="L-BFGS-B",
                    jac=True,
                    options={"maxiter": iterations, "disp": verbose},
                    callback=lambda x: pos_hist.append(x.reshape((nnodes, dim))),
                )

                pos[i] += delta_pos
            yield pos
            t -= dt
            if verbose:
                print(f"{cost(pos,A,k)=}")
    else:
        raise ValueError("method is unknown")


if __name__ == "__main__":
    import numpy as np
    import scipy.sparse
    from src.python.spring_layout import spring_layout

    mat = np.array([[0, 1, 2], [1, 0, -1], [2, -1, 0]])
    mat = scipy.sparse.csr_matrix(mat)
    G = nx.Graph(mat)
    pos = spring_layout(G, iterations=50, seed=0, method="RS")
    print(pos)
