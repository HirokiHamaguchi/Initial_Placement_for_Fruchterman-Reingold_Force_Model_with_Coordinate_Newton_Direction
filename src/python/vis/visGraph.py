from io import BytesIO
from typing import Optional, Union

import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
from matplotlib.colors import Normalize
from PIL import Image


def visGraph(
    G: Union[nx.Graph, nx.DiGraph, nx.MultiGraph, nx.MultiDiGraph],
    pos: dict,
    _dirs: Optional[np.ndarray] = None,
    title: Optional[str] = None,
    savePath: Optional[str] = None,
    node_size: int = 50,
    width: float = 1,
    _colorMap: Optional[np.ndarray] = None,
) -> None:
    plt.figure(figsize=(8, 8))
    plt.axis("equal")
    n = G.number_of_nodes()
    nodes = list(G.nodes())
    cmap = plt.get_cmap("viridis")
    norm = Normalize(vmin=0, vmax=max(n - 1, 1))
    nodeIndices = np.arange(n)
    colorMap = cmap(norm(nodeIndices)) if _colorMap is None else _colorMap
    nx.draw(
        G,
        pos,
        nodelist=nodes,
        node_size=node_size,
        node_color=colorMap,
        width=width,
    )

    scalarMappable = plt.cm.ScalarMappable(norm=norm, cmap=cmap)
    scalarMappable.set_array([])
    posNp = np.array([pos[node] for node in nodes])
    if _dirs is not None and np.any(_dirs != 0.0):
        dirs = _dirs.copy()
        assert posNp.shape == (n, 2)
        assert dirs.shape == (n, 2)
        for i in range(n):
            plt.arrow(
                posNp[i][0],
                posNp[i][1],
                dirs[i, 0],
                dirs[i, 1],
                head_width=0.05,
                head_length=0.05,
                fc="k",
                ec="k",
            )
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

        ram = BytesIO()
        plt.savefig(
            ram, format="png", transparent=False, bbox_inches="tight", pad_inches=0
        )
        ram.seek(0)
        im = Image.open(ram).convert("RGB")
        im2 = im.quantize(
            colors=32,
            method=Image.Quantize.MEDIANCUT,
            dither=Image.Dither.NONE,
        )
        im2.save(savePath, format="PNG", optimize=True, compress_level=9)

        plt.close()
    else:
        plt.show()


if __name__ == "__main__":
    mat = np.array([[1, 2, 3], [2, 1, 4], [3, 4, 1]])
    G = nx.Graph(mat)
    pos = nx.spring_layout(G)
    visGraph(G, pos)
