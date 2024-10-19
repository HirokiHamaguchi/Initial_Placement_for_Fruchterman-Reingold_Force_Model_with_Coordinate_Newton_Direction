from src.python.vis.visGraph import visGraph
import networkx as nx
import numpy as np


def main():
    G = nx.Graph()
    G.add_node(0)
    G.add_node(1)
    G.add_node(2)
    G.add_node(3)
    G.add_node(4)
    G.add_node(5)
    G.add_edge(0, 1)
    G.add_edge(1, 2)
    G.add_edge(2, 3)
    G.add_edge(3, 4)
    G.add_edge(4, 5)
    G.add_edge(5, 0)

    posInit = [
        (-1.75626, 0.85037),
        (-0.519266, -0.441324),
        (1.35645, -0.660642),
        (2.39049, 0.471616),
        (0.870969, 0.764974),
        (-0.433164, 1.37621),
    ]
    posLBFGS = [
        (-1.51971, 0.594647),
        (-0.528408, -0.466651),
        (1.26239, -0.657556),
        (2.24708, 0.319419),
        (0.999081, 1.11524),
        (-0.551213, 1.45611),
    ]
    posFR = [
        (-1.39711, 0.643084),
        (-0.25981, -0.117847),
        (0.985416, -0.47547),
        (2.01626, 0.292993),
        (1.13976, 1.08073),
        (-0.832398, 1.26411),
    ]

    posInit = np.array(posInit)
    posLBFGS = np.array(posLBFGS)
    posFR = np.array(posFR)

    dirLBFGS = posInit - posLBFGS
    dirFR = posInit - posFR

    grad = np.zeros((6, 2))
    for i in range(6):
        delta = (posInit[i] - posInit).T
        distance = np.sqrt((delta**2).sum(axis=0))
        distance = np.where(distance < 0.01, 0.01, distance)
        Ai = np.zeros(6)
        Ai[(i + 1) % 6] = 1
        Ai[(i - 1) % 6] = 1
        grad[i] = (delta * (1 / distance**2 - Ai * distance)).sum(axis=1)

    # visGraph(G, posInit, _dirs=grad, savePath="testInit")
    # visGraph(G, posFR, _dirs=dirFR, savePath="testFR")
    # visGraph(G, posLBFGS, _dirs=dirLBFGS, savePath="testLBFGS")

    # compute forces

    tex = ""
    for i in range(6):
        vals = [
            posInit[i][0],
            posInit[i][1],
            posFR[i][0],
            posFR[i][1],
            posLBFGS[i][0],
            posLBFGS[i][1],
        ]
        tex += str(i) + "/"
        tex += "/".join(list(map(str, vals)))
        tex += ",\n" if i < 5 else "\n"

    print(tex)


if __name__ == "__main__":
    main()
