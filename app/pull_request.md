# Pull Request

Hello, NetworkX developers.

I'm Hiroki Hamaguchi from the University of Tokyo, Japan. We want to contribute to NetworkX by improving the `_sparse_fruchterman_reingold` function, a.k.a. `spring_layout`.

We know that this function is one of NetworkX's core features, so we are very careful about the changes and willing to receive any feedback to refine our approach. We hope that our proposal will benefit the NetworkX community.

## Abstract

The **Fruchterman--Reingold (FR) model** is a framework that assumes both attractive and repulsive forces act between vertices. `nx.draw` uses this model and simulates this force system to perform graph drawing.
Here, `_sparse_fruchterman_reingold` is one of the key internal functions called when **the number of vertices exceeds 500**.

To generate faster and better layouts using the FR force model, we aim to **modify this function** by adopting a different approach from the traditional force-directed method. Specifically, we utilize **SciPy's L-BFGS method** to optimize the layout, similar to `nx.kamada_kawai_layout`.

We can broadly classify our proposed method as an **energy-based algorithm**, but it may be more intuitive to interpret the approach as introducing a kind of momentum into the motion simulation and thus achieving a better layout. Please note that our proposal does not introduce a new model; **the conventional interpretation of the FR model works for most cases**.

This method is based on the research by Hosobe, H. [^hosobe] and our research paper [^our].
Our paper mainly discusses the initial placement with this L-BFGS approach, but this L-BFGS's approach itself is also discussed in detail.
We can update our paper to explain our proposal more clearly if necessary.

### Key Contribution

Our proposed method can usually achieve better visualization with faster computation. The visualization result is shown in the figure.
(The blue `FR` means the current `nx.spring_layout`, and the orange `FR (L-BFGS)` means the our proposed method.)

$n \approx 500$ and 50 iterations:
![graphs_500_50](https://github.com/user-attachments/assets/17691e8f-2b11-4989-ab63-6f3aabd998b6)

$n \approx 600$ and 200 iterations:
![graphs_600_200](https://github.com/user-attachments/assets/3cc9d5ac-fa30-4ac6-92c2-5dbd3b7d9bb3)

### Other Improvements (or Side Effects of Our Proposal?)

The output of our proposed method is theoretically equivalent to the output of the existing method when the graph $G$ is connected and has positive weights. However, the output differs when the graph is not connected or has negative weights due to the current implementation's **ill-defined** extension of the FR algorithm in [^FR].

By **ill-defined**, we mean that the corresponding optimization problem can be unbounded, leading to arbitrary results. This is why our proposed method cannot perfectly reproduce the existing method.

We used the current `nx.spring_layout` in the figure below. The three graphs on the left show the results with the number of iterations set to 50, 500, and 5000, respectively. The rightmost graph contains negative edges.

![unconnected_and_negative_FR](https://github.com/user-attachments/assets/5ad0aa1f-7ed8-41b4-82c1-0ae2a0e2fa42)

As the number of iterations increases, the graphs blow away and become unreadable. The same issue arises in the rightmost graph. This is because the optimization problem behind the FR model is ill-defined in these cases, and thus we cannot reproduce the results by the optimization method.

Our proposed method introduces additional forces for each connected component and uses the absolute values of edge weights, preventing such problems. This approach often produces better results, which may not completely match the existing method. This is a side effect of our proposed method, which can be seen as either an advantage or a disadvantage.

![unconnected_and_negative_L-BFGS](https://github.com/user-attachments/assets/a70d2af7-c862-4c59-83dd-1dababcf2164)

## Algorithm

In this section, we describe our proposed algorithm.
We wrote this content in our research paper [^our], Section 8.3.

![algorithm](https://github.com/user-attachments/assets/9416405e-7ecf-4e46-b309-3e562bd29c32)

## Implementation

Please refer to the diff in the pull request for the implementation details.
Here are a few additional notes:

* We used batch processing to improve computational efficiency. The batch size is 500. As far as we tested, even if we increase the batch size more, the efficiency does not improve or sometimes decreases.
* We set `method="auto"` as the default in `nx.spring_layout` because energy-based approaches like ours are relatively rare for the FR force model. We don't want to confuse users.
* Temporary variables such as `distance2` and `Aid` improve computational speed. Although this may make the code look a bit redundant, it was a deliberate trade-off to enhance performance.
* As discussed later, defining the additional force is one of the most challenging problems, not due to mathematical difficulty but rather the philosophical challenges of keeping the implementation concise and what aesthetic principles to base the design on. After various trials, errors, and considerations, we reached the current implementation.

## Evaluation

We evaluated our proposed method. The code is available at [here](https://github.com/HirokiHamaguchi/Initial_Placement_for_Fruchterman-Reingold_Force_Model_with_Coordinate_Newton_Direction/blob/main/app/test.ipynb).

### NetworkX Graphs

We evaluated our method on various graphs generated by NetworkX with different numbers of vertices: 10, 50, 500, and 600.
We set the number of iterations to 50 or 200. (Note that `nx.kamada_kawai_layout` does not allow specifying the number of iterations.)

![graphs_10_50](https://github.com/user-attachments/assets/58b3f93d-bb41-46af-83cd-212990cd4960)
![graphs_50_50](https://github.com/user-attachments/assets/4b51958d-96cd-431d-9ec4-3d21f5b6cb76)
![graphs_500_50](https://github.com/user-attachments/assets/17691e8f-2b11-4989-ab63-6f3aabd998b6)
![graphs_600_200](https://github.com/user-attachments/assets/3cc9d5ac-fa30-4ac6-92c2-5dbd3b7d9bb3)

From these experiments, we can see that our proposed method (orange, `FR (L-BFGS)`) can, in most cases, achieve better visualization with faster computation compared to the existing method (blue, `FR`).

For small graphs with around 10 vertices, the overhead of the `L-BFGS` method can make it slightly slower, but as the number of vertices increases, this overhead becomes negligible, or our method can even outperform the existing one due to the implementation technique. Since `nx.spring_layout` will use our approach only when the number of vertices exceeds 500 by default, we believe this trade-off is acceptable.

### NetworkX Layouts with Different Methods

We also compared our method with `nx.arf_layout` and `nx.forceatlas2_layout`.
Although these methods do not directly use the FR force model, we think that `FR (L-BFGS)` provides better visualization.

$n \approx 50$ and 100 iterations:
![arf_forceatlas_50_100](https://github.com/user-attachments/assets/140c3696-9a6e-437b-a552-1f69df493a51)

### Graphs from SuiteSparse Matrix Collection

We tested our method on large-scale, real-world graphs from the SuiteSparse Matrix Collection.
The superior performance of our method is evident in these cases as well.

200 iterations:
![ssgetpy_17758_200](https://github.com/user-attachments/assets/c0f35f74-9f87-4983-9510-c8c583e936b4)

For comparision, we listed some of the optimal layouts as follows.

||jagmesh1|1138_bus|dwt_1005|dwt_2680|3elt|USPowerGrid|
|-|-|-|-|-|-|-|
| optimal layout  | ![jagmesh1](https://arxiv.org/html/2412.20317v2/extracted/6103965/individual/vis/opt_jagmesh1.png)  | ![1138_bus](https://arxiv.org/html/2412.20317v2/extracted/6103965/individual/vis/opt_1138_bus.png)  | ![dwt_1005](https://arxiv.org/html/2412.20317v2/extracted/6103965/individual/vis/opt_dwt_1005.png)  | ![dwt_2680](https://arxiv.org/html/2412.20317v2/extracted/6103965/individual/vis/opt_dwt_2680.png)  | ![3elt](https://arxiv.org/html/2412.20317v2/extracted/6103965/individual/vis/opt_3elt.png)  | ![USPowerGrid](http://yifanhu.net/GALLERY/GRAPHS/GIF_SMALL/Pajek@USpowerGrid.gif) |

Since the problem is a nonconvex optimization problem, it is easy to fall into the local optimum solutions, especially in cases such as `dwt_1005`, but depending on the seed, L-BFGS can obtain results close to those of optimal layout.

### Special Cases

We confirmed that our method works correctly in special cases.

#### Unconnected Graphs

![separated_100](https://github.com/user-attachments/assets/40bafef9-9035-4b73-9364-db4e6b619e7f)

#### Very Large Graphs

Although we do not provide specific output results here, we confirmed that our method successfully executed on a cycle graph with 100,000 vertices without causing kernel crashes.

We estimate the necessary memory size for the optimization as follows.

1. If we compute layout for a graph with 100'000 vertices for 100 iterations, it will takes about 55.5 hours with our laptop. (`(dim*|V|^2*iterations)/(computational speed per second) = (2*1e5^2*100)/1e7 = 200'000 seconds = 55.5 hours`. This method of estimation is broadly correct for smaller graphs.) Approximately this size would be the maximum size users would want to calculate.
2. We used batch processing with a batch size of 500, i,e., the size of `delta` in the `_sparse_fruchterman_reingold` is $500 \times 100000 \times 2$ when the number of vertices is 100,000 and `dim=2`.
3. `delta`'s data type is `np.float64`, which requires 8 bytes.
4. The memory size is $500 \times 100000 \times 2 \times 8 \approx 0.74$ GB.
5. In most modern laptops or runtime environments, this memory size is not a problem.

#### 3D Graphs

We confirmed that the method works correctly for graphs with `dim=3`.

![3d](https://github.com/user-attachments/assets/0834873c-b2c4-4e08-bbae-0b70338c2c1b)

#### `fixed` Parameter

Our method also works correctly when the `fixed` and `pos` parameters are specified. The following figure shows the output for a cycle graph with four vertices specified as `fixed`.

![fixed](https://github.com/user-attachments/assets/719c8a79-0fae-468d-98fb-c742a0eae385)

## Discussion

We have demonstrated that our proposed method can, in most cases, achieve better visualization with faster computation compared to the existing method.

As far as we anticipate, there are few major inconveniences that users will experience due to this change. If we must point out one, when increasing the number of nodes stepwise from 100 to 1000 and drawing a graph that includes negative edge weights, discontinuities may occur at 500 with `method="auto"`. How to balance such corner cases against improvements in most situations is a matter of values. We personally think it's not a problem, while we would like to follow the community's intentions.

In general, the most well-known algorithm of the FR force model for large-scale graphs is the multi-level approach, such as [Graphviz's sfdp](https://graphviz.org/docs/layouts/sfdp/). Instead of our approach, we could have implemented this in NetworkX.

However, a complex implementation is not necessarily the best choice for a lightweight library like NetworkX. Furthermore, it is implemented in Python, a language known for its relatively slow execution time.

We believe leveraging SciPy's powerful scientific computing capabilities is a good choice in this situation, and it aligns well with NetworkX's [Mission and Values](https://networkx.org/documentation/stable/developer/values.html#mission-and-values).

We hope our proposal will benefit the NetworkX community.

Thank you for your attention.

[^hosobe]: [Hosobe, H. (2012). Numerical optimization-based graph drawing revisited. In 2012 IEEE Pacific Visualization Symposium (pp. 81-88). IEEE.](https://ieeexplore.ieee.org/iel5/6178307/6183555/06183577.pdf)

[^our]: [Hamaguchi, H., Marumo, N., & Takeda, A. (2024). Initial Placement for Fruchterman -- Reingold Force Model with Coordinate Newton Direction. arXiv preprint arXiv:2412.20317.](https://arxiv.org/abs/2412.20317)

[^FR]: [Fruchterman, T. M., & Reingold, E. M. (1991). Graph drawing by force-directed placement. Software: Practice and experience, 21(11), 1129-1164.](https://onlinelibrary.wiley.com/doi/abs/10.1002/spe.4380211102)
