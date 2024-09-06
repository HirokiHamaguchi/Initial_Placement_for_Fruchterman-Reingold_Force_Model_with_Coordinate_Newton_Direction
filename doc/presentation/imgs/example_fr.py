import networkx as nx
from src.python.vis.visGraph import visGraph

G = nx.generators.random_internet_as_graph(30, seed=7)

pos = nx.spring_layout(G, seed=42, iterations=1000)

visGraph(G, pos, savePath="doc/presentation/imgs/example_fr", node_size=300)
