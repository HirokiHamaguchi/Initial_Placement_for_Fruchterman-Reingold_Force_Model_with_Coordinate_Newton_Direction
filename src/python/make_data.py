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
    main()
