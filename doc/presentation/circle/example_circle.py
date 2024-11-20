import os
import networkx as nx
from src.python.spring_layout import spring_layout
from src.python.vis.visGraph import visGraph

G = nx.Graph()
n = 300
for i in range(n - 1):
    G.add_edge(i, i + 1)
G.add_edge(n - 1, 0)

os.chdir(__file__[: __file__.rfind("/")])

for i, pos in enumerate(spring_layout(G, seed=0, method="FR")):
    pos = {i: pos[i] for i in range(n)}
    visGraph(
        G, pos, node_size=100, title=f"iter: {i+1:02d}", savePath=f"circle-{i+1}.png"
    )
