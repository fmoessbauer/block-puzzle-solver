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

#include "LibMVC/LibNuMVC/numvc.hpp"

// specify which problem should be solved 
//#define PROBLEM_DEMO
//#define PROBLEM_555
//#define PROBLEM_666

// export cube as point cloud
#define CSV_OUTPUT

#ifdef PROBLEM_DEMO
constexpr int num_bricks = 4;
constexpr std::array<size_t, 3> cube_ext{{3UL,2UL,2UL}};
using brick_l  = std::integral_constant<size_t, 3UL>;
#elif PROBLEM_666
constexpr int num_bricks = 54;
constexpr std::array<size_t, 3> cube_ext{{6UL,6UL,6UL}};
using brick_l  = std::integral_constant<size_t, 4UL>;
#else
constexpr int num_bricks = 25;
constexpr std::array<size_t, 3> cube_ext{{5UL,5UL,5UL}};
using brick_l  = std::integral_constant<size_t, 5UL>;
#endif

using mask_t   = blaze::StaticMatrix<int, 3UL, brick_l::value>;
using offset_t = blaze::StaticMatrix<int, 3UL, 1UL>;
using rot_t    = blaze::StaticMatrix<int, 3UL, 3UL>;
using bitset_t = std::bitset<std::get<0>(cube_ext) *
                             std::get<1>(cube_ext) *
                             std::get<2>(cube_ext)>;

using coord_vec_t  = std::vector<mask_t>;
using bitset_vec_t = std::vector<bitset_t>;
using graph_t      = boost::adjacency_matrix<
                       boost::undirectedS,
                       boost::property<boost::vertex_index_t, int>,
                       boost::no_property>;

// rotation matrices
const rot_t rot_x{{1,0,0}, {0,0,-1},{0,1,0}};
const rot_t rot_y{{0,0,1}, {0,1,0}, {-1,0,0}};
const rot_t rot_z{{0,-1,0},{1,0,0}, {0,0,1}};

#ifdef PROBLEM_DEMO
const mask_t brick{{0,0,1},
                   {1,0,0},
                   {0,0,0}};
#elif PROBLEM_666
const mask_t brick{{ 0, 0, 1,-1},
                   { 1, 0, 0, 0},
                   { 0, 0, 0, 0}};
#else
const mask_t brick{{ 0, 0, 1,-1,-2},
                   { 1, 0, 0, 0, 0},
                   { 0, 0, 0, 0, 0}};
#endif

inline bool valid_combination(bitset_t, bitset_t);

/**
 * shift brick according to offset
 */
mask_t move(mask_t matrix, offset_t offset){
  using ones_t = blaze::StaticMatrix<int, 1UL, brick_l::value>;
  ones_t ones(1);
  return matrix + (offset * ones);
}

inline offset_t one_hot(size_t pos){
  offset_t mat(0);
  mat.at(pos, 0) = 1;
  return mat;
}

/**
 * minimum maximum offset in given dimension.
 */
inline std::pair<int, int> minmax_dim(const mask_t & matrix, size_t dim){
  auto row = blaze::row(matrix, dim);
  auto max = blaze::max(row);
  auto min = blaze::min(row);
  return std::make_pair(min, max);
}

/**
 * Shift brick to avoid negative coordinates
 */
mask_t brick_to_block(mask_t matrix){
  for(size_t i = 0; i<matrix.rows(); ++i){
    auto min = minmax_dim(matrix, i).first;
    matrix = move(matrix, one_hot(i) * (-min));
  }
  return matrix;
}

/**
 * Generate all permuations of a block in coord format
 */
coord_vec_t generate_block_permutations(const mask_t & brick){
  std::vector<mask_t> blocks;
  blaze::IdentityMatrix<int> eye(3UL);
  rot_t rx = eye;
  for(int i=0; i<4; ++i){
    rx = rx*rot_x;
    rot_t ry = eye;
    for(int j=0; j<4; ++j){
      ry = ry*rot_y;
      rot_t rz = eye;
      for(int k=0; k<4; ++k){
        rz = rz*rot_z;
        auto rot_brick = rx*ry*rz*brick;
        blocks.push_back(brick_to_block(rot_brick));
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
   for(unsigned int i=0; i<num_masks; ++i){
    for(unsigned int j=(i+1); j<num_masks; ++j){
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
template<typename OutputStream>
void write_dimacs(const graph_t & graph, OutputStream & out){
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
 * export solution as point-cloud in csv format
 */
template<typename OutputStream>
void write_csv_point_cloud(
    const std::vector<int> & cube,
    OutputStream & out)
{
  const int splits_row   = std::get<0>(cube_ext);
  const int splits_col   = std::get<1>(cube_ext);
  const int splits_layer = std::get<0>(cube_ext) * std::get<1>(cube_ext);

  out << "x,y,z,v" << std::endl;
  for(unsigned int i=0; i<cube.size(); ++i){
    int x = i % splits_row;
    int y = (i/splits_row) % splits_col;
    int z = i/splits_layer;
    out << x << "," << y << "," << z << "," << cube[i] << std::endl;
  }
}

/**
 * Check if the combination of both bitsets is collision-free.
 * That means, no two ones are added
 */
inline bool valid_combination(bitset_t a, bitset_t b){
  a&=b;
  return a.none();
}

/**
 * Converts solution given as a path with indices shifted by one
 * to a cube of indices where each position in the cube is
 * marked with the path-index of the covering brick.
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
        cube[i] = x+1;
      }
    }
  }
  return cube;
}

/**
 * Pretty print cube by printing it in layers
 */
template<typename OutputStream>
void write_cube_slices(
    const std::vector<int> & cube,
    OutputStream & out)
{
  const int splits_row   = std::get<0>(cube_ext);
  const int splits_layer = std::get<0>(cube_ext) * std::get<1>(cube_ext);
  const std::string delim(splits_row*3-1, '=');

  for(unsigned int i=0; i<cube.size(); ++i){
    if(i % splits_layer == 0){
      out << delim << " Layer " << (i/splits_layer+1) << std::endl;
    }
    out << std::setw(2)
         << cube[i] << " ";
    if(i % splits_row == (splits_row-1)){
      out << std::endl;
    }
  }
  out << delim << std::endl;
}


int main(int argc, char** argv){
  // generate all permutations block
  auto blocks = std::move(generate_block_permutations(brick));
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
  std::cout << "\n=== Start NuMVC solver ===" << std::endl;
  libmvc::NuMVC solver(os,
      boost::num_vertices(col_graph)-num_bricks,
      std::chrono::seconds(100), true);
  solver.cover_LS(libmvc::NuMVC::default_stats_printer);
  auto solution = std::move(solver.get_independent_set());
  // print solution
  std::cout << "Proposed solution: ";
  for(const int e : solution){
    std::cout << (e-1) << " ";
  }
  std::cout << "\n" << std::endl;

  // check if solution is valid
  if(solver.check_solution() && solution.size() == num_bricks){
    std::cout << "Found solution: " << std::endl;
    // convert solution for pretty printing
    auto solved_cube = cube_of_indices(solution, unique_pos_masks);
    write_cube_slices(solved_cube, std::cout);
#ifdef CSV_OUTPUT
    std::fstream os("solution.csv", std::ios::out);
    write_csv_point_cloud(solved_cube, os);
    os.close();
#endif
  } else {
    std::cout << "Found NO solution: " << std::endl;
  }

  return 0;
}

