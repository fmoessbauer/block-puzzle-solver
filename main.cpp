/**
 * Solve a block-puzzle
 *
 * Field: 9x9x2
 * Blocks:
 *   __
 *  | _|
 *  |_|
 *
 *  Goal: place 6 blocks in this field (fill it)
 */

#include <vector>
#include <unordered_set>
#include <array>
#include <bitset>
#include <algorithm>
#include <iostream>
#include <chrono>
#include <future>

#include <blaze/Blaze.h>

constexpr std::array<size_t, 3> cube_ext{{3UL,2UL,2UL}};
//constexpr std::array<size_t, 3> cube_ext{{5UL,5UL,5UL}};

//using stone_l  = std::integral_constant<size_t, 3UL>;
using stone_l  = std::integral_constant<size_t, 5UL>;
using mask_t   = blaze::StaticMatrix<int, 3UL, stone_l::value>;
using offset_t = blaze::StaticMatrix<int, 3UL, 1UL>;
using rot_t    = blaze::StaticMatrix<int, 3UL, 3UL>;
using bitset_t = std::bitset<std::get<0>(cube_ext) *
                             std::get<1>(cube_ext) *
                             std::get<2>(cube_ext)>;

struct state {
  // index, stone
  using rem_t = std::pair<int, bitset_t>;

  bitset_t           cube;
  std::vector<rem_t> remaining;
  std::vector<int>   path;
  unsigned long long processed;
  unsigned long long processed_last;
  std::chrono::time_point<std::chrono::system_clock> lasthit;
};

const rot_t rot_x{{1,0,0}, {0,0,-1},{0,1,0}};
const rot_t rot_y{{0,0,1}, {0,1,0}, {-1,0,0}};
const rot_t rot_z{{0,-1,0},{1,0,0}, {0,0,1}};

const mask_t stone{{0,0,1},
                   {1,0,0},
                   {0,0,0}};
//const mask_t stone{{ 0, 0, 1,-1,-2},
//                   { 1, 0, 0, 0, 0},
//                   { 0, 0, 0, 0, 0}};

mask_t move(mask_t matrix, offset_t offset){
  using ones_t = blaze::StaticMatrix<int, 1UL, stone_l::value>;
  ones_t ones(1);
  return matrix + (offset * ones);
}

offset_t one_hot(size_t pos){
  offset_t mat(0);
  mat.at(pos, 0) = 1;
  return mat;
}

std::pair<int, int> minmax_dim(const mask_t & matrix, size_t dim){
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

std::vector<bitset_t> filter_bitmasks(const std::vector<bitset_t> & bitsets){
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

inline bool valid_combination(bitset_t a, bitset_t b){
  a.flip();
  b.flip();
  a |= b;
  return a.all();
}

// TODO
bool layer_check(const bitset_t & cube){
  static std::vector<bitset_t> layer_masks;

  // cache masks if necessary
  if(layer_masks.size() == 0){
    bitset_t mask;
    for(unsigned int l=0; l<cube_ext[2]; ++l){
      unsigned int layer_size = cube_ext[0] * cube_ext[1];
      int begin = layer_size * l;
      int end   = layer_size * (l+1);
      for(unsigned i=begin; i<end; ++i){
        mask.set(i);
      }
    }
    layer_masks.push_back(mask);
  }
//  for(unsigned int l=0; l<cube_ext[2]; ++l){
//    cube
//  }
}

bool solve_puzzle_impl(
    int depth,
    state & state,
    int begin,
    int end){
  // reduce cache misses
  decltype(state.remaining) local_remaining(state.remaining);

  for(int i=begin; i<end; ++i){
    auto & pos = local_remaining[i];
    // range is ordered, hence end at first "used" bitmask
    if(std::abs(pos.first) < 0){continue;};
    if(valid_combination(pos.second, state.cube)){
      ++state.processed;
      state.cube|=pos.second;
      // append index to path
      state.path.push_back(pos.first);
      // remove from freelist
      pos.first*=(-1);
//      std::cout << "Try valid stone:" << pos.second
//                << " Result: " << state.cube
//                << " Depth: " << depth << std::endl;

      if(depth == 0){return true;}
      if(solve_puzzle_impl(depth-1, state, 0, local_remaining.size())){return true;}

      //backtrack
      if(depth >= 22){
        auto now  = std::chrono::system_clock::now();
        auto diff = std::chrono::duration_cast<std::chrono::milliseconds>(now - state.lasthit);
        //if(state.processed >= 17189084125){
        //  exit(0);
        //}
        if(diff.count() > 1000){
          auto keys_processed   = state.processed - state.processed_last;
          double mkeys_per_sec  = static_cast<double>(keys_processed / diff.count()) / 1000.0;
          std::cout << "Backtrack at depth:" << depth
                    << ", Remaining:" << state.remaining.size()
                    << ", Processed:" << state.processed 
                    << ", MKey/s:"    << mkeys_per_sec 
                    << ", Cube:" << state.cube << std::endl;
          state.lasthit = now;
          state.processed_last = state.processed;
        }
      }
      state.cube^=pos.second;
      pos.first*=(-1);
      state.path.pop_back();
    }
  }
  return false;
}

bool solve_puzzle(int depth, state & st){
  std::vector<std::future<bool>> tasks;
  int num_tasks    = 1;
  int problem_size = st.remaining.size();
  int task_size    = problem_size / num_tasks;
  std::vector<state> states(num_tasks);

//  for(unsigned int i=0; i<num_tasks; ++i){
//    int begin = task_size*i;
//    int end   = (i==num_tasks-1) ? problem_size : task_size*(i+1);
//    states[i] = st;
//    tasks.emplace_back(
//        std::async(
//          std::launch::async,
//          solve_puzzle_impl, depth, states[i], begin, end));
//  }
//  for(auto & t : tasks){
//    t.wait();
//  }
//  for(unsigned int i=0; i<num_tasks; ++i){
//    bool result = tasks[i].get();
//    if(result){
//      st.path = states[i].path;
//      return true;
//    }
//  }
//  return false;
  return solve_puzzle_impl(depth, st, 0, problem_size);
}


/**
 * Converts solution given as a path with indices shifted by one
 * to a cube of indices where each position in the cube is
 * marked with the path-index of the covering stone.
 */
std::vector<int> cube_of_indices(
    const std::vector<int> & path,
    const std::vector<bitset_t> & all_masks)
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

template<typename It>
void shuffle_range(It a, It b){
  std::random_device rd;
  std::mt19937 g(rd());
       
  std::shuffle(a, b, g);
}

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
  std::cout << "Permutations:" << blocks.size() << std::endl;
  for(const auto & block : blocks){
    std::cout << block << std::endl;
  }

  std::vector<bitset_t> bitmasks;
  std::transform(blocks.begin(), blocks.end(),
                 std::back_inserter(bitmasks), block_to_bitmask);

  filter_bitmasks(bitmasks);

  std::cout << "Bitmasks:" << std::endl;
  for(const auto & mask : bitmasks){
    std::cout << mask << std::endl;
  }

  std::cout << "Move Stones:" << std::endl;
  std::vector<mask_t> all_pos_blocks;
  for(const auto & block : blocks){
    std::array<int, 3> move_steps;
    for(int d=0; d<3; ++d){
      move_steps[d] = cube_ext[d] - minmax_dim(block, d).second;
//      std::cout << "Move steps: " << d << "," << move_steps[d] << std::endl;
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
  std::vector<bitset_t> all_pos_masks;
  all_pos_masks.reserve(all_pos_blocks.size());
  std::transform(all_pos_blocks.begin(), all_pos_blocks.end(),
                 std::back_inserter(all_pos_masks), block_to_bitmask);


  std::cout << "Bitmasks: " << all_pos_masks.size() << std::endl;
  auto unique_pos_masks = std::move(filter_bitmasks(all_pos_masks));
  std::cout << "Unique Bitmasks: " << unique_pos_masks.size() << std::endl;

  for(unsigned int i=0; i<unique_pos_masks.size(); ++i){
    std::cout << unique_pos_masks[i] << ", " << i << std::endl;
  }

  // prepare state for solving puzzle
  state st;
  int index=0;
  // append index+1 to each mask
  std::transform(unique_pos_masks.begin(), unique_pos_masks.end(),
                 std::back_inserter(st.remaining),
                 [&index](const auto & m){
                   return std::make_pair(++index, m);
                 });

  if(solve_puzzle(3, st)){
    std::cout << "Found solution: " << std::endl;
    for(const int e : st.path){
      std::cout << (e-1) << ",";
    }
    std::cout << std::endl;
    auto solved_cube = cube_of_indices(st.path, unique_pos_masks);
    print_cube_slices(solved_cube);
  } else{
    std::cout << "Found NO solution" << std::endl;
  }

  return 0;  
}
