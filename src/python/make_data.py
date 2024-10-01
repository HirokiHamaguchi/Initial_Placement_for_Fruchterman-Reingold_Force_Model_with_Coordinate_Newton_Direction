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


if __name__ == "__main__":
    main()
