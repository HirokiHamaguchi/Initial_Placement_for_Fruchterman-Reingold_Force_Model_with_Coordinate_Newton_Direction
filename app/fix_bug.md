# Fix bug in _sparse_fruchterman_reingold

Dear NetworkX maintainers,

I am Hiroki Hamaguchi, a contributor to NetworkX and the author of [this previous pull request](https://github.com/networkx/networkx/pull/7889).

I sincerely apologize, but I have discovered a bug introduced in that PR. Specifically, when `method="auto"` is specified, the `nx.spring_layout` function does not use the new `_energy_fruchterman_reingold` function as intended.

The fix is trivial, and I would appreciate it if you could kindly review and merge the proposed correction.

Thank you very much for your time and consideration.
