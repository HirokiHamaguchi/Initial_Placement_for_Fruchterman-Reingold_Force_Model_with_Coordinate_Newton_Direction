import os
import networkx as nx
from src.python.spring_layout import spring_layout
from src.python.vis.visGraph import visGraph


def read_file(file_path):
    with open(file_path, "r") as file:
        n, m, _k = map(float, file.readline().strip().split())
        assert n == int(n) and m == int(m)
        n = int(n)
        m = int(m)
        row = []
        col = []
        data = []
        for _ in range(m):
            r, c, d = map(float, file.readline().strip().split())
            assert r == int(r) and c == int(c)
            r = int(r)
            c = int(c)
            row.append(r)
            col.append(c)
            data.append(d)
        positions_size = int(file.readline().strip())
        positions = []
        for _ in range(positions_size):
            position = []
            for _ in range(n):
                x, y = map(float, file.readline().strip().split())
                position.extend([x, y])
            positions.append(position)

    # print(n, m, f"{m / (n * (n - 1) / 2):.3%}")

    Gs = []
    for pos in positions:
        G = nx.Graph()
        for i in range(n):
            G.add_node(i, pos=(pos[2 * i], pos[2 * i + 1]))
        for i in range(m):
            G.add_edge(row[i], col[i], weight=data[i])
        Gs.append(G)

    return Gs


Gs = read_file(f"doc/presentation/jagmesh1/result.out")

method = "CN-L-BFGS"
os.makedirs(f"doc/presentation/jagmesh1/{method}", exist_ok=True)
for i, G in enumerate(Gs):
    pos = {i: G.nodes[i]["pos"] for i in G.nodes}
    visGraph(
        G,
        pos,
        savePath=f"doc/presentation/jagmesh1/{method}/{i}.png",
    )
