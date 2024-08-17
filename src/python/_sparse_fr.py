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
    elif method == "L-BFGS-B" or method == "CG":
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
                assert np.all(delta[i] == 0.0)
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
            method="L-BFGS-B" if method == "L-BFGS-B" else "CG",
            jac=True,
            options={"maxiter": iterations, "disp": verbose},
            callback=lambda x: pos_hist.append(x.reshape((nnodes, dim))),
        )

        for pos in pos_hist:
            yield pos

    elif method == "myCG":
        from scipy.optimize._optimize import _line_search_wolfe12

        # make sure we have a coo_matrix representation
        try:
            A = A.tolil()
        except AttributeError:
            A = (sp.sparse.coo_array(A)).tolil()

        k_inv = 1 / k

        def f(x: np.ndarray) -> float:
            pos = x.reshape((nnodes, dim))
            ret = cost(pos, A, k)
            print(f"{ret=}")
            return ret

        def myfprime(x: np.ndarray) -> np.ndarray:
            pos = x.reshape((nnodes, dim))
            grad = np.zeros((nnodes, dim))
            for i in range(nnodes):
                delta = pos[i] - pos
                assert np.all(delta[i] == 0.0)
                distance = np.linalg.norm(delta, axis=1)
                distance = np.where(distance < 0.01, 0.01, distance)
                distance_inv = 1 / distance
                Ai = A.getrow(i).toarray().flatten()
                coefficient1 = Ai * distance * k_inv - (k * distance_inv) ** 2
                grad[i] = coefficient1 @ delta
            return grad.ravel()

        x0 = pos.ravel()
        gtol = 1e-5
        maxiter = 200
        disp = True
        c1 = 1e-4
        c2 = 0.4

        old_fval = f(x0)
        gfk = myfprime(x0)

        k_cg = 0
        xk = x0
        old_old_fval = old_fval + np.linalg.norm(gfk) / 2

        warnflag = 0
        pk = -gfk
        gnorm = np.amax(np.abs(gfk))

        sigma_3 = 0.01

        while (gnorm > gtol) and (k_cg < maxiter):
            deltak = np.dot(gfk, gfk)

            cached_step = [None]

            def polak_ribiere_powell_step(alpha, gfkp1=None):
                xkp1 = xk + alpha * pk
                if gfkp1 is None:
                    gfkp1 = myfprime(xkp1)
                yk = gfkp1 - gfk
                beta_k = max(0, np.dot(yk, gfkp1) / deltak)
                pkp1 = -gfkp1 + beta_k * pk
                gnorm = np.amax(np.abs(gfkp1))
                return (alpha, xkp1, pkp1, gfkp1, gnorm)

            def descent_condition(alpha, xkp1, fp1, gfkp1):
                # Polak-Ribiere+ needs an explicit check of a sufficient
                # descent condition, which is not guaranteed by strong Wolfe.
                #
                # See Gilbert & Nocedal, "Global convergence properties of
                # conjugate gradient methods for optimization",
                # SIAM J. Optimization 2, 21 (1992).
                cached_step[:] = polak_ribiere_powell_step(alpha, gfkp1)
                alpha, _, pk, gfk, gnorm = cached_step

                # Accept step if it leads to convergence.
                if gnorm <= gtol:
                    return True

                # Accept step if sufficient descent condition applies.
                return np.dot(pk, gfk) <= -sigma_3 * np.dot(gfk, gfk)

            try:
                alpha_k, _, _, old_fval, old_old_fval, gfkp1 = _line_search_wolfe12(
                    f,
                    myfprime,
                    xk,
                    pk,
                    gfk,
                    old_fval,
                    old_old_fval,
                    c1=c1,
                    c2=c2,
                    amin=1e-100,
                    amax=1e100,
                    extra_condition=descent_condition,
                )
            except RuntimeError:
                # Line search failed to find a better solution.
                warnflag = 2
                break

            # Reuse already computed results if possible
            if alpha_k == cached_step[0]:
                alpha_k, xk, pk, gfk, gnorm = cached_step
            else:
                alpha_k, xk, pk, gfk, gnorm = polak_ribiere_powell_step(alpha_k, gfkp1)

            k_cg += 1

            yield xk.reshape((nnodes, dim))

        fval = old_fval
        if disp:
            if warnflag == 2:
                msg = "pr_loss"
            elif k_cg >= maxiter:
                warnflag = 1
                msg = "maxiter"
            elif np.isnan(gnorm) or np.isnan(fval) or np.isnan(xk).any():
                warnflag = 3
                msg = "nan"
            else:
                msg = "success"
            print("         WarnFlag and message: %d" % warnflag, msg)
            print("         Current function value: %f" % fval)
            print("         Iterations: %d" % k)

    elif method == "SCG":
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
                    assert t == 1, "change np.log for t!=1"
                    delta = x - pos
                    delta[i] = 0.0
                    distance = np.linalg.norm(delta, axis=1)
                    distance = np.where(distance < 0.01, 0.01, distance)
                    distance_inv = 1 / distance
                    Ai = A.getrow(i).toarray().flatten()
                    cost = np.sum(
                        Ai * distance**3 * (1 / 3) * k_inv - (k**2) * np.log(distance)
                    )
                    grad = (Ai * distance * k_inv - (k * distance_inv) ** 2) @ delta

                    return cost, grad

                pos[i] = sp.optimize.minimize(
                    cost_fun_i,
                    pos[i],
                    jac=True,
                    options={"maxiter": 3},
                ).x

            yield pos

            if verbose:
                print(f"{cost(pos,A,k)=}")
    elif method == "SCCG":
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
                    assert t == 1, "change np.log for t!=1"
                    delta = x - pos
                    delta[i] = 0.0
                    distance = np.linalg.norm(delta, axis=1)
                    distance = np.where(distance < 0.01, 0.01, distance)
                    distance_inv = 1 / distance
                    Ai = A.getrow(i).toarray().flatten()
                    cost = np.sum(
                        Ai * distance**3 * (1 / 3) * k_inv - (k**2) * np.log(distance)
                    )
                    grad = (Ai * distance * k_inv - (k * distance_inv) ** 2) @ delta

                    return cost, grad

                pos[i] = sp.optimize.minimize(
                    cost_fun_i,
                    pos[i],
                    jac=True,
                    options={"maxiter": 3},
                ).x

            yield pos

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
    pos = spring_layout(G, iterations=50, seed=0, method="BFGS")
    print(pos)
