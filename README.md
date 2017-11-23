# Block Puzzle Solver

This program solves the 5x5x5 brick-cube-puzzle in less than one second using a sophisticated algorithmic approach.
The implementation is both highly efficient as well as easy to understand.
The presented techniques are suited for many optimization problems and also can be used for educational purposes.
Compared to traditional strategies like backtracking, this approach is more than two magnitudes faster.

While this implementation currently only solves the provided puzzle as well as the "T-Cube" (6x6x6) it can easily be extended to solve any space filling problem using blocks of arbitrary shapes.

![Wooden Block Puzzle](https://gist.githubusercontent.com/fmoessbauer/ed15ccb82cf2c4626fdbd17de7145506/raw/29f6d2a8fa59593d60f674093203e0028385af7f/cube_puzzle_small.jpg)

*Real-world wooden block puzzle. All bricks are of equal shape.*

## Compiling

This program depends on [Blaze](https://bitbucket.org/blaze-lib/blaze) and the [Boost Graph Library](http://www.boost.org/doc/libs/1_65_1/libs/graph/doc/index.html). For compiling we recommend using cmake:

```bash
export Blaze_DIR="/path/to/blaze"
export Boost_DIR="/path/to/boost"
# for 5x5x5 problem
cmake .
# for 6x6x6 problem ("T-Cube")
cmake -DTASK=PROBLEM_666 .
# compile
make
```

## Running

```bash
# solve puzzle
./block-puzzle

# convert solution to easy to understand images
python3 cube-printer.py solution.csv
```

## How it works

Algorithmically solving problems in 3D spaces is a tricky and error prone task.
This is something humans are good at, but not computers.

Hence, to reduce complexity we map the 5x5x5 cube to a bitmask consisting of 125 bits.
Each bit then denotes if the corresponding grid cell of the cube is occupied.
These bitmasks require only a small amount of memory and collision detection of two bricks is trivial using bitwise logic.
This transformes the 3D problem space into a 1D memory block, where the problem is solved iff all bits are set.

### Generate Combinations

At first all combinations of all (unique) bricks have to be calculated. In our case we have just one brick.
For easier calculations the brick is stored as a *3 x m* matrix where each block of the brick is represented as coordinate vector.
Often this format is also called point cloud.
The matrix for the brick shown above is then given as follows (column vectors can be swapped):

```
-2 -1  0  0  1
 0  0  1  0  0
 0  0  0  0  0
```

**Rotation**

To get all orientations of the brick, it has to be rotated about every axis.
This is done by multiplying the brick-matrix with [rotation matrices](https://en.wikipedia.org/wiki/Rotation_matrix).
As the brick might be symmetric with respect to an axis or a point, this rotation process can return duplicates.
These would only differ in an offset. To avoid this, the resulting bricks are shifted so that the smalles coordinate
in each axis is 0. After that duplicates are removed.

**Translation**

The next step is to shift each brick to all possible locations in the cube. This can also be calculated easily
by adding an offset matrix to each point cloud matrix. Again, to reduce the number of position-shifted bricks,
duplicates are filtered.

### Convert to Bitmasks

The coordinate format makes calculation of the rotations and transformations easy, but it is not well suited for checking if two points share the same location (grid cell). Hence, we transform the brick into a 125 bit bitmask where a bit is set if the corresponding grid cell is occupied by this block. Here, the serialization order is not important as we do rely on neighborhood relations.

### Calculate Collision Graph

Using the bitmask format, collision detection becomes trivial. To check if two bitmasks are collision free (no two true-bits collide), the following function can be used:

```cpp
bool valid_combination(bitmask a, bitmask b){
  a&=b;
  return a.none();
}
```

Now, the goal is to find 25 bitmasks out of all candidates which do not collide.
Mathematically spoken this 25 bitmasks must not collide pairwise.
As the number of combinations (25 out of ~ 960) is roughly 3e74 exhaustive search is not possible.
This also applies to the backtracking approach as shown below.

Hence, we rely on a different approach: We calculate the collision graph of all bitmasks by checking each pair of masks:
If the masks collide, an edge is added between both vertices (indices of bitmask).
This results in a graph with one vertex per mask and edges between every two vertices where the corresponding bitmasks collide.
As the collision function is symmetric, the graph is undirected.
Theoretically there would be an reflexiv edge for each vertex, but we do not store them as each vertex must occur at most once.

**Example**

As the graph is dense (it has many more edges than vertices) we store it as [adjaceny matrix](https://en.wikipedia.org/wiki/Adjacency_matrix).

```
Bitmasks                  0 1 2 3
                         ---------
0 | 0 1 0 0 1  Graph  0 | - - - - |
1 | 1 1 0 0 0  ====>  1 | x - - - |
2 | 0 0 1 0 0         2 |     - - |
3 | 0 1 1 1 0         3 | x   x - |
                         ---------
```

The collision graph (in edge-list format) of the example is then: `{(0,1),(0,3),(2,3)}`.

### Solve MVC

Using the collision graph, a valid (collision free) subset of bricks is equal to
the [minimal-vertex-cover](https://en.wikipedia.org/wiki/Vertex_cover) (MVC) of the graph.
In fact, the solution is the maximum independent set, which is the inverse problem of MVC.

![Collision Graph](https://gist.githubusercontent.com/fmoessbauer/ed15ccb82cf2c4626fdbd17de7145506/raw/baa1f075d3025c983f5e2c0ea89166693ebbfe36/sample-graph.png)

The blue colored vertices are the minimal-vertex-cover of the example graph.
Hence, the white vertices are the independent set which denotes a collision free combination.

Although MVC is an [NP-Complete](https://en.wikipedia.org/wiki/NP-completeness) problem,
there exist good heuristics which enable algorithms to solve it quite efficiently in most cases.
Here we use the [Improved-NuMVC](https://github.com/fmoessbauer/Improved-NuMVC) solver, originally developed and published by [Shaowei Cai](http://lcs.ios.ac.cn/~caisw/).

### Reconstruct Solution

The output of the MVC algorithm is the maximum independent set, given as a set of vertex indices.
This then has to be mapped to the corresponding bricks.
Depending on how the solution should be presented, it is converted to different representations:

**Point Cloud**

A 125 element point cloud where each point is colored according to the index of the brick that occupies this grid cell.
This format is useful to visualise the solution using the python script or other visualisation tools.

![Sample Solution](https://gist.githubusercontent.com/fmoessbauer/ed15ccb82cf2c4626fdbd17de7145506/raw/f262253f348004b5139e97a2d1bcd29b55bd3b07/solution_layers.png)

*Visualization of a solution, using the python script. CSV file is available in the folder `data`.*

**Layered Cube**

A cube of indices (printed in layers) where the value at each position denotes the index of the block that occupies this grid cell.
This representation can easily be read by humans.

## Other Approaches

Other approaches for solving the problem have been considered:

### Backtracking

**Idea**: Add bricks to the cube as long as possible:

Bricks which do not collide with the partial solution are added to the cube (by adding the bitmasks).
If no further non colliding brick can be added and not all bits of the cube are set, backtrack.

While this solution seems natural to solve the problem, even with highly optimized and cache friendly code only 22 to 23 bricks could be added in acceptable time (hours).
This approach already uses a good presorting of the bricks to reduce backtracking to a minimum.
While it could be further improved by using constraint propagation, we dropped this approach as the implementation of the MVC approach is much easier and the solution is found super fast.

### Machine Learning Approaches

- Specialized solutions which optimize after multiple runs might be achived through [reinforcement learning](https://arxiv.org/pdf/1704.01665.pdf).
- Genetic Algorithms / Cross Entropy Search with a density cost function or segmented solutions (solve individually as dense `1x5x5` layers).

