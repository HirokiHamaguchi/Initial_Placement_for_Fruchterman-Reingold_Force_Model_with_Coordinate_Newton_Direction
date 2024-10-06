import os
import networkx as nx
import scipy.io
import scipy.sparse


def nx2mtx(G: nx.Graph, filename: str):
    adj = nx.adjacency_matrix(G)
    coo = scipy.sparse.coo_matrix(adj)
    scipy.io.mmwrite(filename, coo)


def main():
    os.chdir("data")
    for n in range(1, 12 + 1):
        G = nx.balanced_tree(2, n)
        nx2mtx(G, f"balanced_tree_2_{n}.mtx")
    for n in [50, 100, 200, 500, 1000]:
        G = nx.barabasi_albert_graph(n, 2)
        nx2mtx(G, f"barabasi_albert_2_{n}.mtx")
        G = nx.circulant_graph(n, [1])
        nx2mtx(G, f"circulant_{n}.mtx")
        G = nx.circulant_graph(n, [1, 2, 3])
        nx2mtx(G, f"circulant_{n}_1_2_3.mtx")


if __name__ == "__main__":
    main()
