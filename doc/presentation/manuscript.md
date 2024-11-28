# manuscript

## 1

Thank you for the introduction.
Then, I'd like to start today's my presentation. The presentation title is Initial Placement for Fruchterman--Reingold Force Model with Coordinate Newton Direction.
And my supervisor is Professor Akiko Takeda.

## 2

To begin with, I'll introduce and summarize our study.

## 3

Firstly, please let me briefly introduce what the graph drawing is.
The graph is a mathematical structure representing pairwise relationships between objects, composed from vertices $V$ and edges $E$.
There are many examples of graphs, such as social network graphs, railroad graphs, computer network graphs, and so on.
In particular, graph drawing is one of the most fundamental tasks in information visualization.

You are probably most familiar with graph drawing by CONNECTED PAPERS.
You might have used this cite before.
When you specify a paper (like this purple one), it draws a graph with vertices as papers and edges as citations.
The drawn graph allows us to understand the flow and relevance of the research visually.
Thus, correct and beautiful Graph drawing is essential not only in the theoretical field but also in the practical field.

And, Force-directed graph drawing is one of the most popular methods for drawing graphs.

## 4

The force-directed graph drawing uses a force model composed of particles with forces acting between them.
Among the various force models, the Fruchterman--Reingold (FR) force model is the most prominent one, regarded as flexible, intuitive, and simple.

Here is an example by NetworkX, the popular graph library in Python.
So we import NetworkX as nx (at here).
`nx.draw`, as known as `nx.spring_layout`, is the function to draw a graph using the FR force model, and NetworkX draws the graphs by FR algorithm with 50 iterations.

Unfortunately, the FR algorithm has an evident drawback.
Let us consider a cycle graph example.
When the number of vertices is 10, it well visualizes the graph in 0.2 seconds. There is no problem.
When the number of vertices is 20, it starts to visualize the graph with a tangled structure in 0.2 seconds.
And, When the number of vertices is 500, although it takes as long as 11.5 seconds and the graph is just a simple cycle, the visualization result is too messy, too tangled, and too twisted.

This is exactly what we want to address in this study.

## 5

Now, let me explain why this problem happens.

The FR force model uses a force model with two forces: the attractive force $F^a_{i,j}$ (this one) and the repulsive force as $F_r$ (this one).

And here, $w_{i,j}$ is the weight of the edge between vertices $i$ and $j$. It has a positive value if and only if an edge between $i$ and $j$ exists.
For every pair of vertices with an edge, the attractive force pulls them together.
For every pair of vertices, whether there is an edge or not, the repulsive force pushes them apart.

When we consider only one pair of vertices, there is an equilibrium between these two forces at some distance $d$.

Then, the FR algorithm, the algorithm behind the `nx.draw`, visualizes a graph by simulating the system and seeking the equilibrium of ALL of these forces.
This is how the FR algorithm works.

## 6

Now, we can see why the problem we mentioned before can happen.
The FR algorithm has issues that make it challenging to achieve high-quality visualizations for large-scale graphs.
The most critical issue is that twist slows down the simulation process.
The twist means unnecessary folded and tangled structures in the visualized graph.

As this gif shows, the FR algorithm from a random initial placement, the most standard way to draw a graph, stagnates at the twisted structure.

Although the optimal graph structure is just a simple cycle, the twist weakens or diminishes the forces, causing stagnation.

Further, the FR algorithm suffers from high computational complexity when directly applied, $\order{\abs{V}^2}$ per iteration where $\abs{V}$ is the number of vertices.

Thus, we need to find a way to avoid the twist with a small computational complexity.

## 7

In previous studies, several ways to improve the FR algorithm have been proposed.

One of the strategies is to directly accelerate the simulation process with the forces, in other words, the optimization process of the energy function.
Recent research has accelerated the process in various ways, one of which is to utilize numerical optimization methods.

L-BFGS algorithm, a family of quasi-Newton methods, is one such approach and is reported to be effective for graph drawing.
Although this method can overcome the twist problem to some extent, it may sometimes require many iterations.
It may fail to achieve the optimal visualization within a limited number of iterations, such as 50 iterations.
Moreover, this method merely treats the problem as a general optimization problem (like this. Instead of the variable $X \in \bbR^{2\ times n}, it uses the variable $X \in \bbR^{2n}$).
Thus, directly applying the L-BFGS algorithm may not be the optimal strategy.

## 8

Providing initial placement in the pre-processing step is another strategy.
Indeed, a pre-processing step with Simulated Annealing (SA) is also known to be effective since SA can deal with twist issues and leads to a better visualization combined with the FR algorithm. It also utilizes the inherent graph structure.
However, this work is restricted to unweighted graphs, limited to a circle initial placement, inefficient due to the random swapping, the only neighborhood in the SA step, and sometimes ignoring the graph's sparsity.
It leaves significant room to improve the effectiveness and extend the applicability.

So, we aim to provide an initial placement that can avoid the twist, accelerate the subsequent optimization process, and extend the applicability to weighted graphs.

## 9

In this situation, we propose a new initial placement for the FR force model as depicted in this figure.
We provide an initial placement (this one) with fewer twists than random placement within a short time, accelerating the subsequent optimization process.

We can use both the FR algorithm (this one) and the L-BFGS algorithm (this one) to obtain the final placement.

This work extends the applicability of the initial placement idea to larger-scale, weighted, and complicated structure graphs.

This is the main contribution of our study.
From now on, we will explain the details of our proposed method.

## 10

Ok, then, let me explain the details of our proposed method.

## 11

Firstly, we formulate the problem.
Recall that the FR force model uses the attractive force $F^a_{i,j}$ and the repulsive force $F_r$.

Integrating the forces defines its scalar potential, the energy function.
The energy function between vertices $i$ and $j$ is defined as $E_{i,j}$, the sum of the attractive energy $E^a_{i,j}$ and the repulsive energy $E_r$.

In the FR algorithm, we seek the equilibrium of the energy function by simulating the system.
In contrast to this ordinary approach, we minimize the energy function to seek the equilibrium.

Note that the minimum of $f$ yields the equilibrium positions since $\nabla f(X)$ corresponds to the forces.

## 12

Next, to introduce our proposed method, we explain the simplified problem.

Even at the expense of accuracy, obtaining an approximate solution quickly is crucial for the initial placement.
So, to obtain it, we simplify the problem into a more manageable and well-behaved discrete optimization problem, which is defined as (this one).

In this problem, we assign the vertices $V$ to the fixed set of point $Q^hex$, the hexagonal lattice, and minimize this function.
(So we assign the vertices to the hexagonal lattice like this blue one, orange one, and green one. And the objective function is minimized from this one to this one. We aim to find the optimal assignment like this.)

Here, we explain why we use this simplified problem.
We start from simplifying the original problem, the problem (1).
Since graphs generally have sparsity, we separate the attractive and repulsive energy to leverage the sparsity as this one.

## 13

Then, we convert the second term into a constraint.

We can compute the objective function with $\abs{E}$ terms rather than $\abs{V}^2$ terms.
This conversion does not lose the essence of the problem too much.
Since $E^\mathrm{r}(d)=-k^2\log{d}$ is a convex function such that it decreases monotonically concerning $d$, for sufficiently large $d$, the value of $-k^2\log{d}$ does not grow excessively. For too small $d$, we can prevent the divergence of the energy function by setting $\epsilon$.

Still, the problem involves $\order{\abs{V}^2}$ constraints, which negates the advantage of computing the objective function with $\order{\abs{E}}$ complexity.
To simplify further, we incorporate the initial placement idea mentioned in the introduction.
We simplify the problem with a fixed initial placement $Q$ as this one.

We use $Q$ such that the points are separated by at least $\epsilon$.
By fixing the possible point placement in advance, we can skip the check of the $\order{\abs{V}^2}$ constraints in the problem, reducing the computational complexity to $\order{\abs{E}}$ and thus offering significant speedup.

Although this $Q$ is arbitrary, since we only consider the attractive energy $f^a(X)$, the point set $Q$ should be as dense as possible.
Thus, we use the hexagonal lattice, one of the closet packing, as $Q$.

## 14

So, this is the summary of the first half.

We want to solve this problem,
and we simplify it into this problem.

We use the hexagonal lattice as the initial placement $Q$.

Now, we will explain how to solve this problem in the second half.

## 15

To solve the problem, we use the coordinate Newton direction.
So, to begin with, please let me introduce it in this slide.

Let us consider a strictly convex function $f$.
The second order approximation of $f$ at $x_0$ is this one.
The argmin $x^*$ of this approximation must satisfy the stationary condition, the derivative of it equals zero.
Thus, we obtain this equation and this red term is called the Newton direction in general.
Here, please note that the Hessian matrix $\nabla^2 f(x_0)$ is positive definite since $f$ is strictly convex, so we can compute the inverse Hessian and this is a descent direction.
This is an important property for the Newton direction.

Although the Newton direction provides a critical step in iterative methods, it requires the computation of the inverse Hessian $\nabla^2 f(x_0)^{-1} \in \bbR^{n \times n}$, posing a high computational cost for large-scale problems.

Still, we can leverage the concept of the Newton direction in a different manner, the stochastic coordinate descent with the coordinate Newton direction.
Instead of computing the inverse Hessian $\nabla^2 f(x_0)^{-1}$ in the entire variable space $\bbR^n$, we limit the variable $x$ to its coordinate block $x_i$ with fewer dimensions, and compute this $f_i$'s Newton direction, the coordinate Newton direction.
Since the coordinate Newton direction computation is much cheaper than that of the Newton direction, we can repeat this procedure many times.
In general, this idea is known as stochastic coordinate descent or Randomized Subspace Newton (RSN) in a broader context.

In particular, this coordinate Newton direction has an apparent natural affinity to the problem.
By taking the position $x_i$ of the vertex $i$ as the coordinate block, we can compute the Newton direction.
Although directly applying this idea to the problem is challenging, we leverage this coordinate Newton direction to propose our algorithm.

## 16

Now, we will explain the proposed method.
First, by randomly taking a vertex $i$ from $V$, we compute the coordinate Newton direction for $x_i$, and its gradient and Hessian is computed as these.

Since $f^a_i$ is strictly convex, which is different form the energy function $E_{i,j}$ and $f^a$, we can utilize the coordinate Newton direction.

## 17

Then, using this coordinate Newton direction, we update the position of $i$ as $x_i^{new}$ in the hexagonal lattice $Q^hex$, and by repeating this process, we can solve the problem.
So, here, we explain how to update the position $x_i$ in an one iteration.

If our problem is just a general continuous optimization problem, we can directly apply the coordinate Newton direction.
However, $x_i^new$ may not be in the hexagonal lattice $Q^hex$. So we need to round it to the nearest point in $Q^hex$.
This figure shows this rounding operation. Select the vertex $i$ at $x_i$, and compute the coordinate Newton direction as this arrow.
Then, by rounding it to the nearest point in $Q^hex$, we determine the new position $x_i^{new}$. To keep the injectivity of the assignment $\pi$, we need to swap the vertices if necessary.

In this process, we empirically found that adding a small random noise vector to the coordinate Newton direction is effective.
This randomness is similar to the SA one, and it can help to escape from local minima and to explore the solution space more effectively.

So, this it the one step of our proposed method, and we repeat this process for some iterations.

## 18

Finally, the obtained placement (this one) could be too small or too large, since we did not case about the scale. Thus, as the final step, we rescale by $c^*$ to obtain the initial placement  (like this) and apply L-BFGS algorithm for the final placement.

So we have to the optimal scaling factor $c^*$. Recalling the original problem for the FR force model, the problem is to minimize this function $\phi(c)$. Here and here is the $c$.

Since this function is convex, we can find the optimal scaling factor $c^*$ by solving this problem as this one.

This value can be computed in the $\order{\abs{E}}$ complexity.
Thus, as far as we rescale the placement by this $c^*$, we can select any $\epsilon$ to define the hexagonal lattice $Q^hex$. So we can eliminate the parameter $\epsilon$ from the problem.

## 19

Ok, now, this slide shows the pseudo-code of the proposed method.
It is just a repeat of the previous slides.

We first randomly select vertex $i$ from $V$, then compute the coordinate Newton direction and round it to $x_i^new$ in $Q^hex$, and update the position of $i$.
We repeat this process until some iterations.
Finally, we rescale by $c^*$ to obtain our proposed initial placement.

So, this is the outline of our proposed method.

## 20

Next, we will show the experimental results.
