from pathlib import Path

import networkx as nx
import numpy as np
import scipy.io
import scipy.sparse
import scipy.sparse as sp
import ssgetpy

DATA_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = DATA_DIR.parent
MATRIX_NAMES_PATH = PROJECT_ROOT / "data/_matrixNames.txt"


def load_matrix(mat):
    path = Path(mat.download(extract=True)[0]) / f"{mat.name}.mtx"
    assert path.exists(), f"Missing matrix file: {path}"
    A = scipy.io.mmread(path)
    return A if sp.issparse(A) else sp.coo_matrix(A)


def make_graph(A):
    A = A.tocsr(copy=True)
    A.setdiag(0)
    A.eliminate_zeros()
    return nx.from_scipy_sparse_array(A)


def process_matrix(mat, force=False):
    if mat.rows!=mat.cols:
        return None

    A = load_matrix(mat).tocoo() # type: ignore
    if not force and np.any(A.data<0):
        return None
    A.data = np.abs(A.data)

    if not force and not np.allclose((A-A.T).data, 0):
        return None
    A = ((A + A.T) / 2).tocoo()

    G = make_graph(A)
    if  not nx.is_connected(G):
        return None

    scipy.io.mmwrite(DATA_DIR / f"{mat.name}.mtx", A)
    return mat.name, G.number_of_nodes()

def main1():
    auto_mats = ssgetpy.search(rowbounds=(None, 1000), limit=10000)

    manual_names = [
        "can_144",
        "jagmesh1",
        "dwt_1005",
        "1138_bus",
        "bcsstk13",
        "add20",
        "dwt_2680",
        "3elt",
        "USPowerGrid",
        "bcspwr10",
    ]

    mat_names = []

    for mat in auto_mats:
        result = process_matrix(mat)
        if result is not None and mat.rows<=1000:
            mat_names.append(result)

    for name in manual_names:
        mat = ssgetpy.search(name)[0]
        result = process_matrix(mat, force=True)
        assert result is not None, f"{name} failed processing"
        # don't add these names to MATRIX_NAMES_PATH
        # since the experiments to be run are different.
        # mat_names.append(result)

    with MATRIX_NAMES_PATH.open("w") as f:
        f.write("\n".join(name for name, _ in mat_names) + "\n")

    num_nodes = dict(mat_names)

    print(f"cnt={len(mat_names)}")
    print(num_nodes)


# handmade matrices

def nx2mtx(G: nx.Graph, filename: str):
    adj = nx.adjacency_matrix(G)
    coo = scipy.sparse.coo_matrix(adj)
    scipy.io.mmwrite(DATA_DIR/filename, coo)


def main2():
    for n in range(1, 12 + 1):
        G = nx.balanced_tree(2, n)
        nx2mtx(G, f"btree{n}.mtx")
    for n in [50, 100, 200, 300, 500]:
        G = nx.circulant_graph(n, [1])
        nx2mtx(G, f"cycle{n}.mtx")
    for a, b, c in [(3, 2, 6)]:
        G = nx.Graph()
        for i in range(a):
            G.add_edge(i, (i + 1) % a)
        for i in range(a):
            for j in range(b):
                idx = a + c * (b * i + j)
                G.add_edge(i, idx)
                for k in range(c):
                    G.add_edge(idx + k, idx + (k + 1) % c)
        nx2mtx(G, f"_circleHandmade_{a}_{b}_{c}.mtx")


if __name__ == "__main__":
    main1()
    main2()
