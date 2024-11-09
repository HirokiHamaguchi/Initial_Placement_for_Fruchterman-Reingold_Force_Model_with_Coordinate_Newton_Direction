import networkx as nx
import numpy as np
import matplotlib.pyplot as plt
from typing import Union


def visGraph(
    G: Union[nx.Graph, nx.DiGraph, nx.MultiGraph, nx.MultiDiGraph],
    pos: dict,
    _dirs: np.ndarray = None,
    title: str = None,
    savePath: str = None,
    node_size: int = 50,
    width: float = 1,
) -> None:
    plt.figure(figsize=(8, 8))
    plt.axis("equal")
    cmap = plt.get_cmap("jet")
    n = G.number_of_nodes()
    colorMap = np.array([cmap(i / (n - 1)) for i in range(n)])
    nx.draw(G, pos, node_size=node_size, node_color=colorMap, width=width)
    if _dirs is not None and np.any(_dirs != 0.0):
        dirs = _dirs.copy()
        assert pos.shape == (n, 2)
        assert dirs.shape == (n, 2)
        for i in range(n):
            plt.arrow(
                pos[i][0],
                pos[i][1],
                dirs[i, 0],
                dirs[i, 1],
                head_width=0.05,
                head_length=0.05,
                fc="k",
                ec="k",
            )
    if title:
        plt.title(title)
    if savePath:
        plt.gca().set_axis_off()
        plt.subplots_adjust(top=1, bottom=0, right=1, left=0, hspace=0, wspace=0)
        plt.margins(0, 0)
        plt.savefig(savePath, dpi=75)
        plt.close()
    else:
        plt.show()


if __name__ == "__main__":
    mat = np.array([[1, 2, 3], [2, 1, 4], [3, 4, 1]])
    G = nx.Graph(mat)
    pos = nx.spring_layout(G)
    visGraph(G, pos)
