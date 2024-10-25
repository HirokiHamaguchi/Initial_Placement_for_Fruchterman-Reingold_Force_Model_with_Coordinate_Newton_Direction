import math
import random

import matplotlib.pyplot as plt
import networkx as nx
import numpy as np


def main():
    G = nx.Graph()
    for i in range(10):
        if i < 4:
            G.add_node(
                i,
                x=0.2 * math.cos(i * 2 * math.pi / 4) - 0.3,
                y=0.2 * math.sin(i * 2 * math.pi / 4),
            )
        else:
            G.add_node(
                i,
                x=0.3 * math.cos((i - 1) * 2 * math.pi / 6) + 0.4,
                y=0.3 * math.sin((i - 1) * 2 * math.pi / 6),
            )
    G.add_edge(0, 4)
    for i in range(4):
        for j in range(i + 1, 4):
            G.add_edge(i, j)
    for i in range(4, 10):
        for j in range(i + 1, 10):
            G.add_edge(i, j)

    cMap1 = list(range(G.number_of_nodes()))
    random.shuffle(cMap1)
    cMap2 = list(range(G.number_of_nodes()))
    random.shuffle(cMap2)

    # nx.draw(
    #     G,
    #     with_labels=True,
    #     pos={i: (G.nodes[i]["x"], G.nodes[i]["y"]) for i in G.nodes},
    # )
    # plt.savefig("iso.png")
    # assert False

    print("\\begin{document}")
    print("\\begin{tikzpicture}")
    print(f"\\foreach ", end="")
    for cIdx in range(10):
        c = chr(ord("a") + cIdx)
        print(f"\\x{c}/\\y{c}/", end="")
    print("\\xS/\\yS/\\ord in {")
    for seed in [0, 1]:
        np.random.seed(seed)
        x = np.random.rand(10)
        y = np.random.rand(10)
        xMax, xMin = x.max(), x.min()
        yMax, yMin = y.max(), y.min()
        x = 1.2 * (x - xMin) / (xMax - xMin) - 0.5
        y = (0.6 * (y - yMin) / (yMax - yMin) - 0.3) * math.cos(math.pi / 6)
        for i in range(10):
            print(f"{x[i]:.3f}/{y[i]:.3f}/", end="")
        if seed == 0:
            print("-1.3/+0.6/0,")
        else:
            print("-1.3/-0.6/1,")

    for j in range(2):
        for i in range(10):
            x = G.nodes[i]["x"]
            y = G.nodes[i]["y"]
            print(f"{x:.3f}/{y:.3f}/", end="")
        if j == 0:
            print("+1.3/+0.6/0,")
        else:
            print("+1.3/-0.6/1")

    print("}{")

    print("\\begin{scope}[shift={{(\\xS,\\yS)}}]")

    for i in range(10):
        print(
            f"\\coordinate ({chr(ord('a')+i)}) at (\\x{chr(ord('a')+i)},\\y{chr(ord('a')+i)});"
        )

    for i, j in G.edges:
        print(f"\\draw ({chr(ord('a')+i)}) -- ({chr(ord('a')+j)});")

    print("\\ifthenelse{\\ord=0}{")
    for i in range(10):
        c = cMap1[i]
        print(
            f"\\node[draw=none,circle,inner sep=1pt,fill=c{c}] at ({chr(ord('a')+i)}) {{}};"
        )
    print("}{")
    for i in range(10):
        c = cMap2[i]
        print(
            f"\\node[draw=none,circle,inner sep=1pt,fill=c{c}] at ({chr(ord('a')+i)}) {{}};"
        )
    print("}")

    print("\\end{scope}")

    print("}")
    print("\\end{tikzpicture}")
    print("\\end{document}")


if __name__ == "__main__":
    main()
