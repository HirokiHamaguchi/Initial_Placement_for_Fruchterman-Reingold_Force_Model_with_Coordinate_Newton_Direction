# Reply1

Thank you very much for your thoughtful and encouraging feedback!
I truly appreciate your time and effort in reviewing our proposal.

Regarding your questions, we address them below.
We will commit the changes to the PR once the keywords for `method` are finalized.

## Testing

We have also considered testing and agree that adding tests would be beneficial.

One possible approach is to extend the existing test at [this location](https://github.com/networkx/networkx/blob/44aab23f18f635153c63486e5ba821a8cbba835f/networkx/drawing/tests/test_layout.py#L149).

If necessary, we could add the following:

```diff
    def test_adjacency_interface_numpy(self):
        A = nx.to_numpy_array(self.Gs)
        pos = nx.drawing.layout._fruchterman_reingold(A)
        assert pos.shape == (6, 2)
        pos = nx.drawing.layout._fruchterman_reingold(A, dim=3)
        assert pos.shape == (6, 3)
        pos = nx.drawing.layout._sparse_fruchterman_reingold(A)
        assert pos.shape == (6, 2)
+       pos = nx.drawing.layout._sparse_fruchterman_reingold(A, method="L-BFGS")
+       assert pos.shape == (6, 2)

    def test_adjacency_interface_scipy(self):
        A = nx.to_scipy_sparse_array(self.Gs, dtype="d")
        pos = nx.drawing.layout._sparse_fruchterman_reingold(A)
        assert pos.shape == (6, 2)
        pos = nx.drawing.layout._sparse_spectral(A)
        assert pos.shape == (6, 2)
        pos = nx.drawing.layout._sparse_fruchterman_reingold(A, dim=3)
        assert pos.shape == (6, 3)
+       pos = nx.drawing.layout._sparse_fruchterman_reingold(A, dim=3, method="L-BFGS")
+       assert pos.shape == (6, 3)
```

Besides, although the main concern will be whether SciPy’s L-BFGS runs without errors, since SciPy already thoroughly tests this, we believe additional safeguards are unnecessary on our end.

## L-BFGS Meaning

L-BFGS stands for **Limited-memory Broyden-Fletcher-Goldfarb-Shanno**, a well-known [optimization algorithm](https://en.wikipedia.org/wiki/Limited-memory_BFGS).

We chose this name for consistency with SciPy's terminology, as it is already used in the `method` parameter of `scipy.optimize.minimize` ("L-BFGS-B", where "B" stands for Box-constrained, which is not relevant).

Reference: [SciPy Documentation](https://docs.scipy.org/doc/scipy/reference/optimize.minimize-lbfgsb.html)

## Method Naming

We considered using "force" for the existing method (FR) and "energy" for our proposal (L-BFGS) as potential alternatives.

However, these terms might be misleading in a strict technical sense. Since "L-BFGS" is already used in `scipy.optimize.minimize`, it may be somewhat familiar to users -- though admittedly, it is not the most intuitive name.

### Adding a "Notes" Section

Adding a "Notes" section sounds like a great idea! Would something like the following work?

```diff
    method : str, optional (default='auto')
        The method used to compute the layout.
        If 'FR', the force-directed Fruchterman--Reingold algorithm [1] is used.
-       If 'L-BFGS', the energy-based optimization algorithm [2] is used
-       with absolute values of edge weights and additional forces per connected component.
+       If 'L-BFGS', the energy-based optimization algorithm [2] is used.
        If 'auto', we use 'FR' when len(G) < 500 and 'L-BFGS' otherwise.

    Returns
    -------
    pos : dict
        A dictionary of positions keyed by node

    Examples
    --------
    >>> G = nx.path_graph(4)
    >>> pos = nx.spring_layout(G)
    >>> # suppress the returned dict and store on the graph directly
    >>> _ = nx.spring_layout(G, seed=123, store_pos_as="pos")
    >>> nx.get_node_attributes(G, "pos")
    {0: array([-0.61520994, -1.        ]), 1: array([-0.21840965, -0.35501755]), 2: array([0.21841264, 0.35502078]), 3: array([0.61520696, 0.99999677])}

    # The same using longer but equivalent function name
    >>> pos = nx.fruchterman_reingold_layout(G)

+    Notes
+    -----
+    method="FR" supports negative edge weights and uses only attractive and repulsive forces.
+    method="L-BFGS" uses absolute edge weights and additional forces per connected component,
+    making graph drawing a well-defined optimization problem. This enables the use of the 
+    efficient optimization algorithm: Limited-memory Broyden-Fletcher-Goldfarb-Shanno (L-BFGS).
+    You can see the differences in the [Gallery](todo???).

 References
    ----------
    .. [1] Fruchterman, Thomas MJ, and Edward M. Reingold.
           "Graph drawing by force-directed placement."
           Software: Practice and experience 21, no. 11 (1991): 1129-1164.
           http://dx.doi.org/10.1002/spe.4380211102
    .. [2] Hamaguchi, Hiroki, Naoki Marumo, and Akiko Takeda.
           "Initial Placement for Fruchterman--Reingold Force Model With Coordinate Newton Direction."
           arXiv preprint arXiv:2412.20317 (2024).
           https://arxiv.org/abs/2412.20317
    """
```

## SciPy Version

[One of the features we used](https://docs.scipy.org/doc/scipy/reference/generated/scipy.sparse.csgraph.connected_components.html) has been available since SciPy 0.11.0.

Additionally, the another feature is already used in `nx.kamada_kawai_layout`, compatibility should not be an issue.

If we understand correctly, the proposed test should confirm compatibility with the required SciPy version, reducing concerns on that front.

## Layout Gallery

The idea of contributing to the layout gallery is interesting.
But we're unsure where to start. Are you referring to adding images to [this section](https://networkx.org/documentation/stable/auto_examples/index.html)?

If you have specific ideas or requests, please let us know -- we’d be happy to help!

---

Again, thank you so much for your review and kind words!
We look forward to hearing your thoughts on these points.
