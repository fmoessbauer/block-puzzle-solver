/**
 * Solve a block-puzzle
 *
 * Field: 5x5x5
 * Blocks:
 *     ____
 *    |_  _|
 *      |_|
 *
 *  Goal: place 25 blocks in this field (fill it)
 */

#include <array>
#include <bitset>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <chrono>

#include <blaze/Blaze.h>
#include <boost/graph/adjacency_matrix.hpp>

#include "Improved-NuMVC/tsewf.hpp"

//#define PROBLEM_DEMO

#ifdef PROBLEM_DEMO
constexpr int num_stones = 4;
constexpr std::array<size_t, 3> cube_ext{{3UL,2UL,2UL}};
using stone_l  = std::integral_constant<size_t, 3UL>;
#else
constexpr int num_stones = 25;
constexpr std::array<size_t, 3> cube_ext{{5UL,5UL,5UL}};
using stone_l  = std::integral_constant<size_t, 5UL>;
#endif

using mask_t   = blaze::StaticMatrix<int, 3UL, stone_l::value>;
using offset_t = blaze::StaticMatrix<int, 3UL, 1UL>;
using rot_t    = blaze::StaticMatrix<int, 3UL, 3UL>;
using bitset_t = std::bitset<std::get<0>(cube_ext) *
                             std::get<1>(cube_ext) *
                             std::get<2>(cube_ext)>;

using coord_vec_t  = std::vector<mask_t>;
using bitset_vec_t = std::vector<bitset_t>;
using graph_t  = boost::adjacency_matrix<
                  boost::undirectedS,
                  boost::property<boost::vertex_index_t, int>,
                  boost::no_property>;

const rot_t rot_x{{1,0,0}, {0,0,-1},{0,1,0}};
const rot_t rot_y{{0,0,1}, {0,1,0}, {-1,0,0}};
const rot_t rot_z{{0,-1,0},{1,0,0}, {0,0,1}};

#ifdef PROBLEM_DEMO
const mask_t stone{{0,0,1},
                   {1,0,0},
                   {0,0,0}};
#else
const mask_t stone{{ 0, 0, 1,-1,-2},
                   { 1, 0, 0, 0, 0},
                   { 0, 0, 0, 0, 0}};
#endif

inline bool valid_combination(bitset_t, bitset_t);

mask_t move(mask_t matrix, offset_t offset){
  using ones_t = blaze::StaticMatrix<int, 1UL, stone_l::value>;
  ones_t ones(1);
  return matrix + (offset * ones);
}

inline offset_t one_hot(size_t pos){
  offset_t mat(0);
  mat.at(pos, 0) = 1;
  return mat;
}

inline std::pair<int, int> minmax_dim(const mask_t & matrix, size_t dim){
  auto row = blaze::row(matrix, dim);
  auto max = blaze::max(row);
  auto min = blaze::min(row);
  return std::make_pair(min, max);
}

/**
 * Shift stone to avoid negative coordinates
 */
mask_t stone_to_block(mask_t matrix){
  for(size_t i = 0; i<matrix.rows(); ++i){
    auto min = minmax_dim(matrix, i).first;
    matrix = move(matrix, one_hot(i) * (-min));
  }
  return matrix;
}

/**
 * Generate all permuations of a block in coord format
 */
coord_vec_t generate_block_permutations(const mask_t & stone){
  std::vector<mask_t> blocks;
  blaze::IdentityMatrix<int> eye(3UL);
  rot_t rx = eye;
  for(int i=0; i<4; ++i){
    rx = rx*rot_x;
    rot_t ry = eye;
    for(int j=0; j<4; ++j){
      ry = ry*rot_y;
      rot_t rz = eye;
      for(int k=0; k<4; k+=2){
        rz = rz*rot_z;
        if(j==k){continue;} // already calculated
        auto rot_stone = rx*ry*rz*stone;
        blocks.push_back(stone_to_block(rot_stone));
      }
    }
  }
  return blocks;
}

/**
 * Convert a block in coordinate format
 * to a bitmask in cube-format
 */
bitset_t block_to_bitmask(const mask_t & block){
  const size_t numcols = block.columns(); // n

  bitset_t bitset;
  for(size_t i=0; i<numcols; ++i){
    auto col = blaze::column(block, i);
    bitset.set(col[0]
             + col[1] * std::get<0>(cube_ext)
             + col[2] * std::get<0>(cube_ext) * std::get<1>(cube_ext));
  }
  return bitset;
}

/**
 * return bitsets without duplicates
 */
bitset_vec_t filter_bitmasks(const bitset_vec_t & bitsets){
  using comp_t = std::pair<unsigned int, std::string>;
  std::vector<comp_t> bitmasks;
  bitmasks.reserve(bitsets.size());
  // pair bitmasks with index
  unsigned int i=0;
  std::transform(bitsets.begin(), bitsets.end(),
                 std::back_inserter(bitmasks),
                 [&i](const auto & m){
                  return std::make_pair(i++,m.to_string());});
  std::sort(bitmasks.begin(), bitmasks.end(),
      [&bitmasks](const auto & a, const auto & b){
        return a.second < b.second;
      });
  // bitmasks are small so move them too (however only interested in indices)
  auto endit = std::unique(bitmasks.begin(), bitmasks.end(),
               [](const auto & a, const auto & b){
                 return a.second == b.second;
               });
  // put unique bitsets into a new vector
  // nice side effect: bitmasks are sorted
  std::vector<bitset_t> result;
  result.reserve(std::distance(bitmasks.begin(), endit));
  std::for_each(bitmasks.begin(), endit, 
      [&result, &bitsets](const auto & e){
        result.push_back(bitsets[e.first]);
      });
  return result;
}

/**
 * Generate collision graph of all bitsets
 */
graph_t generate_collision_graph(const bitset_vec_t & bitsets){
  const auto num_masks = bitsets.size();
  graph_t graph(num_masks);
 
  // generate vertices
   for(int i=0; i<num_masks; ++i){
    for(int j=(i+1); j<num_masks; ++j){
      // no collision
      if(!valid_combination(bitsets[i], bitsets[j])){
        boost::add_edge(i,j, graph);
      }
    }
  }
  return graph;
}

/**
 * export boost graph to dimacs graph format
 */
void write_dimacs(const graph_t & graph, std::ostream & out){
  auto num_edges    = boost::num_edges(graph);
  auto num_vertices = boost::num_vertices(graph);
  boost::graph_traits<graph_t>::edge_iterator ei, ei_end;

  out << "p edge " << num_vertices << " " << num_edges << std::endl;
  for (tie(ei, ei_end) = boost::edges(graph); ei != ei_end; ++ei){
    auto source = boost::source(*ei, graph);
    auto target = boost::target(*ei, graph);
    auto sidx   = boost::get(boost::vertex_index, graph, source);
    auto tidx   = boost::get(boost::vertex_index, graph, target);
    out << "e " << sidx+1 << " " << tidx+1 << std::endl;
  }
}

/**
 * Check if the combination of both bitsets is collision-free.
 * That means, no two ones are added
 */
inline bool valid_combination(bitset_t a, bitset_t b){
// faster that -(a n b)
  a.flip();
  b.flip();
  a|=b;
  return a.all();
}

/**
 * Converts solution given as a path with indices shifted by one
 * to a cube of indices where each position in the cube is
 * marked with the path-index of the covering stone.
 */
std::vector<int> cube_of_indices(
    const std::vector<int> & path,
    const bitset_vec_t & all_masks)
{
  size_t size = all_masks[0].size();
  std::vector<int> cube(size);

  for(unsigned int x=0; x<path.size(); ++x){
    const auto & s = all_masks[path[x]-1];
    for(unsigned int i=0; i<size; ++i){
      if(s.test(i)){
        cube[i] = x; 
      }
    }
  }
  return cube;
}

/**
 * Pretty print cube by printing it in layers
 */
void print_cube_slices(
    const std::vector<int> & cube)
{
  const int splits_row   = std::get<0>(cube_ext);
  const int splits_layer = std::get<0>(cube_ext) * std::get<1>(cube_ext);
  const std::string delim(splits_row*3-1, '=');

  std::cout << delim << std::endl;
  for(unsigned int i=0; i<cube.size(); ++i){
    std::cout << std::setw(2)
              << cube[i] << " ";
    if(i % splits_row == (splits_row-1)){
      std::cout << std::endl;
    }
    if(i % splits_layer == (splits_layer-1)){
      std::cout << delim << std::endl;
    }
  }
}


int main(int argc, char** argv){
  // generate all permutations block
  auto blocks = std::move(generate_block_permutations(stone));
  std::cout << "Stone permutations: " << blocks.size() << std::endl;

  // convert permutations to bitmasks
  std::vector<bitset_t> bitmasks;
  std::transform(blocks.begin(), blocks.end(),
                 std::back_inserter(bitmasks), block_to_bitmask);
  // remove duplicates
  filter_bitmasks(bitmasks);

  std::cout << "Stone bitmasks (unique): " << bitmasks.size() << std::endl;

  // move every block to all possible positions in cube
  std::vector<mask_t> all_pos_blocks;
  for(const auto & block : blocks){
    std::array<int, 3> move_steps;
    for(int d=0; d<3; ++d){
      move_steps[d] = cube_ext[d] - minmax_dim(block, d).second;
    }
    for(int x=0; x<move_steps[0]; ++x){
      for(int y=0; y<move_steps[1]; ++y){
        for(int z=0; z<move_steps[2]; ++z){
          offset_t offset({{x},{y},{z}});
          auto moved_block = move(block, offset);
          all_pos_blocks.push_back(moved_block);
        }
      }
    }
  }

  // convert blocks to bitmasks
  std::vector<bitset_t> all_pos_masks;
  all_pos_masks.reserve(all_pos_blocks.size());
  std::transform(all_pos_blocks.begin(), all_pos_blocks.end(),
                 std::back_inserter(all_pos_masks), block_to_bitmask);

  
  std::cout << "Position bitmasks: " << all_pos_masks.size() << std::endl;
  auto unique_pos_masks = std::move(filter_bitmasks(all_pos_masks));
  std::cout << "Position bitmasks (unique): " << unique_pos_masks.size() << std::endl;

  // generate collision graph
  auto col_graph  = std::move(generate_collision_graph(unique_pos_masks));
  // write graph to file for external solver
  // std::ofstream os("file.mis", std::ios::out);

  // or pipe directly into solver
  std::stringstream os;
  write_dimacs(col_graph, os);
  
  // initialize solver and solve
  TSEWF solver(os, boost::num_vertices(col_graph)-num_stones, 1000);
  solver.cover_LS();
  auto solution = std::move(solver.get_independent_set());
  // print solution
  for(const int e : solution){
    std::cout << (e-1) << " ";
  }
  std::cout << std::endl;

  // check if solution is valid
  if(solver.check_solution() && solution.size() == num_stones){
    std::cout << "Found solution: " << std::endl;
    // convert solution for pretty printing
    auto solved_cube = cube_of_indices(solution, unique_pos_masks);
    print_cube_slices(solved_cube);
  } else {
    std::cout << "Found NO solution: " << std::endl;
  }

  return 0;  
}

