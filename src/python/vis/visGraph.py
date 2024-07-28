import networkx as nx
import numpy as np
import matplotlib.pyplot as plt
from typing import Union


def visGraph(
    G: Union[nx.Graph, nx.DiGraph, nx.MultiGraph, nx.MultiDiGraph],
    pos: dict,
    title: str = None,
    savePath: str = None,
) -> None:
    cmap = plt.get_cmap("jet")
    n = G.number_of_nodes()
    colorMap = np.array([cmap(i / (n - 1)) for i in range(n)])
    nx.draw(G, pos, node_size=50, node_color=colorMap)
    if title:
        plt.title(title)
    if savePath:
        plt.savefig(savePath, bbox_inches="tight")
        plt.close()
    else:
        plt.show()


if __name__ == "__main__":
    mat = np.array([[1, 2, 3], [2, 1, 4], [3, 4, 1]])
    G = nx.Graph(mat)
    pos = nx.spring_layout(G)
    visGraph(G, pos)
