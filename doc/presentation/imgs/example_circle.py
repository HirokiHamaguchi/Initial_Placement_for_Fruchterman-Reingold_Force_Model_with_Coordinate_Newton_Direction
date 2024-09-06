import os
import numpy as np
import networkx as nx
from src.python.spring_layout import spring_layout
from src.python.vis.visGraph import visGraph
from src.python.vis.visAnimate import visAnimate

G = nx.Graph()
n = 300
for i in range(n - 1):
    G.add_edge(i, i + 1)
G.add_edge(n - 1, 0)

os.chdir(__file__[: __file__.rfind("\\")])

# visAnimate(
#     G,
#     [(np.copy(res),) for res in spring_layout(G, seed=0, method="FR")],
#     matrixName="circle",
#     methodName="FR",
# )

os.makedirs("circle", exist_ok=True)
for i, pos in enumerate(spring_layout(G, seed=0, method="FR")):
    visGraph(G, pos, title=f"iter: {i+1:02d}", savePath=f"circle/circle-{i+1}.png")
