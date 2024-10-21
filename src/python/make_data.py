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
        nx2mtx(G, f"btree{n}.mtx")
    for n in [50, 100, 200, 300, 500]:
        G = nx.circulant_graph(n, [1])
        nx2mtx(G, f"cycle{n}.mtx")


if __name__ == "__main__":
    main()
