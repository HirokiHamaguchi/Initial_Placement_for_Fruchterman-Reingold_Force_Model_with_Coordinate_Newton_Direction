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
    if mat.rows != mat.cols:
        return None

    A = load_matrix(mat).tocoo()  # type: ignore
    if not force and np.any(A.data < 0):
        return None
    A.data = np.abs(A.data)

    if not force and not np.allclose((A - A.T).data, 0):
        return None
    A = ((A + A.T) / 2).tocoo()

    G = make_graph(A)
    if not nx.is_connected(G):
        if not force:
            return None
        original_n = len(G)
        largest_cc = max(nx.connected_components(G), key=len)
        G = G.subgraph(largest_cc).copy()
        A = nx.adjacency_matrix(G)
        print(f"subgraph {original_n} -> {len(G)}")

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
        if result is not None:
            mat_names.append(result)

    for name in manual_names:
        mat = ssgetpy.search(name)[0]
        result = process_matrix(mat, force=True)
        assert result is not None, f"{name} failed processing"
        # don't add these names to MATRIX_NAMES_PATH
        # since the experiments to be run are different.

    with MATRIX_NAMES_PATH.open("w") as f:
        f.write("\n".join(name for name, _ in mat_names) + "\n")

    print(f"cnt={len(mat_names)}")


def nx2mtx(G: nx.Graph, filename: str):
    adj = nx.adjacency_matrix(G)
    coo = scipy.sparse.coo_matrix(adj)
    scipy.io.mmwrite(DATA_DIR / filename, coo)


def main2():
    # handmade matrices
    for n in range(1, 12 + 1):
        G = nx.balanced_tree(2, n)
        nx2mtx(G, f"btree{n}.mtx")
    for n in [50, 100, 200, 300, 500]:
        G = nx.circulant_graph(n, [1])
        nx2mtx(G, f"cycle{n}.mtx")


def relabel_int(G: nx.Graph) -> nx.Graph:
    return nx.convert_node_labels_to_integers(G, ordering="sorted")


def make_cylinder_grid(m: int, n: int) -> nx.Graph:
    """
    m x n cylindrical grid.
    Periodic in the second dimension, open in the first dimension.

    Nodes are initially (i, j), then relabeled to 0, ..., mn-1.
    """
    G = nx.Graph()

    for i in range(m):
        for j in range(n):
            G.add_node((i, j))

            # vertical edge
            if i + 1 < m:
                G.add_edge((i, j), (i + 1, j))

            # periodic horizontal edge
            G.add_edge((i, j), (i, (j + 1) % n))

    return relabel_int(G)


def _sierpinski_recursive(order: int):
    """
    Returns a Sierpinski gasket graph and its three outer corner nodes.
    """
    if order == 0:
        G = nx.Graph()
        G.add_edges_from([(0, 1), (1, 2), (2, 0)])
        return G, (0, 1, 2)  # top, left, right

    G0, c0 = _sierpinski_recursive(order - 1)
    n = len(G0)

    H = nx.disjoint_union_all([G0, G0.copy(), G0.copy()])

    A = (c0[0], c0[1], c0[2])
    B = tuple(x + n for x in c0)
    C = tuple(x + 2 * n for x in c0)

    # Identify adjacent corners of the three smaller triangles.
    H = nx.contracted_nodes(H, A[1], B[0], self_loops=False)
    H = nx.contracted_nodes(H, A[2], C[0], self_loops=False)
    H = nx.contracted_nodes(H, B[2], C[1], self_loops=False)

    outer_corners = (A[0], B[1], C[2])

    mapping = {old: i for i, old in enumerate(H.nodes())}
    H = nx.relabel_nodes(H, mapping, copy=True)

    return H, tuple(mapping[v] for v in outer_corners)


def make_sierpinski(order: int) -> nx.Graph:
    G, _ = _sierpinski_recursive(order)
    return relabel_int(G)


def make_grid(m: int, n: int) -> nx.Graph:
    """
    Ordinary m x n grid graph.
    Added as a simple regular-graph baseline in response to the reviewer.
    """
    G = nx.grid_2d_graph(m, n)
    return relabel_int(G)


def main3():
    # additional graphs from SuiteSparse / Florida collection
    manual_names = [
        "jagmesh8",
        "crack",
        "USpowerGrid",
        "wiki-Vote",
    ]

    for name in manual_names:
        mat = ssgetpy.search(name)[0]
        result = process_matrix(mat, force=True)
        assert result is not None, f"{name} failed processing"

    generated_graphs = {
        "cylinder_32_32.mtx": make_cylinder_grid(32, 32),
        "sierpinski_06.mtx": make_sierpinski(6),
        "grid_32_32.mtx": make_grid(32, 32),
    }

    for filename, G in generated_graphs.items():
        assert nx.is_connected(G), f"{filename} is not connected"
        nx2mtx(G, filename)


if __name__ == "__main__":
    main1()
    main2()
    main3()
