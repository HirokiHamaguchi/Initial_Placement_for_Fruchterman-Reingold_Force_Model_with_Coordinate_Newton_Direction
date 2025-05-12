"""
===================
Spring Layout
===================

Draw graphs using the three different spring layout algorithms.

The spring layout generally yield by the Fruchterman-Reingold force-directed algorithm.
This algorithm treats edges as springs holding nodes close while treating nodes as repelling objects, and simulates the system until it reaches an equilibrium state.
The algorithm also terminates when a maximum number of iterations is reached.

NetworkX offers three different kind of methods based on the same theory behind the spring layout algorithm:

* `nx.spring_layout(G, method="force")`

  * The default for graphs with fewer than 500 nodes in `nx.spring_layout`.
  * Direct implementation of the Fruchterman-Reingold force-directed algorithm.
  * Can handle negative edge weights.

* `nx.spring_layout(G, method="energy")`

  * The default for graphs with more than or equal to 500 nodes in `nx.spring_layout`.
  * Instead of force simulation, it solves an energy-based optimization problem, taking the absolute value of negative edge weights.
  * Generally produces higher-quality layouts compared to `force`.
  * Uses gravitational forces acting on each connected component to prevent diverging.

* `nx.nx_agraph.graphviz_layout(G, prog="sfdp")`

  * Uses `sfdp` from GraphViz to compute the layout.
  * Can also be done with `nx.nx_pydot.graphviz_layout(G, prog="sfdp")`.
  * Employs a multilevel approach to progressively apply force-directed graph drawing.
  * Requires separate installation of GraphViz, see [here](https://networkx.org/documentation/stable/reference/drawing.html#module-networkx.drawing.nx_agraph) for details.
  * Especially fast for large graphs.
"""

import time
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import networkx as nx

graphs = [
    (nx.cycle_graph(100), "cycle"),
    (nx.grid_2d_graph(25, 25), "grid_2d"),
    (nx.gnp_random_graph(100, 0.001), "gnp_random"),
]

colors = {"force": "tab:blue", "energy": "tab:orange", "sfdp": "tab:green"}

fig, axes = plt.subplots(3, 3, figsize=(9, 9))

for i, (G, name) in enumerate(graphs):
    results = []

    t0 = time.perf_counter()
    pos = nx.spring_layout(G, method="force", seed=0, iterations=100)
    dt = time.perf_counter() - t0
    results.append(("force", pos, dt))

    t0 = time.perf_counter()
    pos = nx.spring_layout(G, method="energy", seed=0, iterations=100)
    dt = time.perf_counter() - t0
    results.append(("energy", pos, dt))

    t0 = time.perf_counter()
    pos = nx.nx_agraph.graphviz_layout(G, prog="sfdp")
    dt = time.perf_counter() - t0
    results.append(("sfdp", pos, dt))

    for j, (mname, pos, dt) in enumerate(results):
        nx.draw(G, pos=pos, ax=axes[j, i], node_color=colors[mname], node_size=20)
        title = (f"{name}\n" if i == 0 else "") + f"{dt:.2f}s"
        axes[j, i].set_title(title, fontsize=20)

handles = [mpatches.Patch(color=color, label=key) for key, color in colors.items()]
fig.legend(handles=handles, loc="upper center", ncol=3, fontsize=25)

plt.tight_layout(rect=(0, 0, 1, 0.9))
plt.savefig("test.png")
