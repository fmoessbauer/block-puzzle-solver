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
#include <array>
#include <bitset>
#include <algorithm>
#include <iostream>

#include <blaze/Blaze.h>

constexpr std::array<size_t, 3> cube_ext{{3UL,2UL,2UL}};

using stone_l  = std::integral_constant<size_t, 3UL>;
using mask_t   = blaze::StaticMatrix<int, 3UL, stone_l::value>;
using offset_t = blaze::StaticMatrix<int, 3UL, 1UL>;
using rot_t    = blaze::StaticMatrix<int, 3UL, 3UL>;
using bitset_t = std::bitset<std::get<0>(cube_ext) *
                             std::get<1>(cube_ext) *
                             std::get<2>(cube_ext)>;

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
  const size_t numrows = block.rows(); // 3
  const size_t numcols = block.columns(); // n

  bitset_t bitset;
  for(size_t i=0; i<numcols; ++i){
    auto col = blaze::column(block, i);
    bitset.set(col[0]
             + col[1] * std::get<0>(cube_ext)
             + col[2] * std::get<1>(cube_ext) * std::get<2>(cube_ext));
  }
  return bitset;
}

int num_unique_permutations(const std::vector<mask_t> & blocks){
  std::vector<bitset_t> bitmasks;
  std::transform(blocks.begin(), blocks.end(),
                 std::back_inserter(bitmasks), block_to_bitmask);

  std::vector<unsigned long long> bits;
  std::transform(bitmasks.begin(), bitmasks.end(),
                 std::back_inserter(bits), [](const auto & m){return m.to_ullong();});

  std::sort(bits.begin(), bits.end());
  bits.erase(std::unique(bits.begin(), bits.end()), bits.end());
  return bits.size();
}

std::vector<bitset_t> filter_bitmasks(const std::vector<bitset_t> & bitsets){
  using comp_t = std::pair<unsigned int, unsigned long long>;
  std::vector<comp_t> bitmasks;
  bitmasks.reserve(bitsets.size());
  // pair bitmasks with index
  unsigned int i=0;
  std::transform(bitsets.begin(), bitsets.end(),
                 std::back_inserter(bitmasks),
                 [&i](const auto & m){
                  return std::make_pair(i++,m.to_ullong());});
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

int main(int argc, char** argv){
  std::vector<mask_t> blocks;
  auto s = stone;
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

  std::cout << "Bitmasks:" << std::endl;
  for(const auto & mask : bitmasks){
    std::cout << mask << std::endl;
  }

  std::cout << "Unique Bitmasks: " << num_unique_permutations(blocks) << std::endl;

  std::cout << "Move Stones:" << std::endl;
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
  std::vector<bitset_t> all_pos_masks;
  all_pos_masks.reserve(all_pos_blocks.size());
  std::transform(all_pos_blocks.begin(), all_pos_blocks.end(),
                 std::back_inserter(all_pos_masks), block_to_bitmask);


  for(const auto & mask : all_pos_masks){
    std::cout << mask << std::endl;
  }

  std::cout << "Bitmasks: " << all_pos_masks.size() << std::endl;
  auto unique_pos_masks = std::move(filter_bitmasks(all_pos_masks));
  std::cout << "Unique Bitmasks: " << unique_pos_masks.size() << std::endl;

  return 0;  
}
