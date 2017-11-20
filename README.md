# Block Puzzle Solver

## Compiling

This program depends on [Blaze](https://bitbucket.org/blaze-lib/blaze) and the [Boost Graph Library](http://www.boost.org/doc/libs/1_65_1/libs/graph/doc/index.html). For compiling I recommend using cmake: `cmake .`

## Running

```bash
# solve puzzle
./block-puzzle

# convert solution to easy to understand images
python3 cube-printer.py solution.csv
```

## How it works


### Generate Permutations

- Generate all permutations of the brick
- move every permutation to all possible locations in cube
- filter duplicates

### Convert to Bitmasks

- convert bricks from coordinate-format to bitmask for efficient collision detection

### Calculate Collision Graph

- calculate collision graph of bitmasks

### Solve MVC

- solve [minimal-vertex-cover](https://en.wikipedia.org/wiki/Vertex_cover)

### Reconstruct Solution

- convert independent set to cube

## Other Approaches

- Backtracking

