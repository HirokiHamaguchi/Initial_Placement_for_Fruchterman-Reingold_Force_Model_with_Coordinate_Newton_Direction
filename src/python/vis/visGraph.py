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
    posNp = np.array(list(pos.values()))
    minX = np.min(posNp[:, 0])
    maxX = np.max(posNp[:, 0])
    minY = np.min(posNp[:, 1])
    maxY = np.max(posNp[:, 1])
    centerX = (minX + maxX) / 2
    centerY = (minY + maxY) / 2
    maxDiff = max(maxX - minX, maxY - minY) / 2
    plt.gca().set_xlim(centerX - maxDiff * 1.03, centerX + maxDiff * 1.03)
    plt.gca().set_ylim(centerY - maxDiff * 1.03, centerY + maxDiff * 1.03)
    if title:
        plt.title(title)
    if savePath:
        plt.gca().set_axis_off()
        plt.subplots_adjust(top=1, bottom=0, right=1, left=0, hspace=0, wspace=0)
        plt.margins(0, 0)
        plt.savefig(savePath)
        plt.close()
    else:
        plt.show()


if __name__ == "__main__":
    mat = np.array([[1, 2, 3], [2, 1, 4], [3, 4, 1]])
    G = nx.Graph(mat)
    pos = nx.spring_layout(G)
    visGraph(G, pos)
