#ifndef N5_ODD_CHARGE_H
#define N5_ODD_CHARGE_H

#include <stdint.h>

#include "n5_core.h"

#define N5_ODD_CHARGE_CELL_COUNT (N5_V2_SIZE * N5_V2_SIZE)
#define N5_ODD_CHARGE_MAX_NULLITY N5_V3_SIZE

/*
 * 一个 GL(3,2) x Stab(K) reduced translation-charge class。
 *
 * dimension：
 *   对应 scalar flow space 中 canonical subspace 的维数，范围 0..3。
 *
 * basis_coordinates：
 *   canonical GF(2) RREF basis，坐标相对于当前 K 的 scalar nullspace basis。
 *
 * subspace_orbit_size：
 *   该 scalar subspace 在顶层 Stab(K) 下的 orbit 大小。
 *
 * weight：
 *   完整 GL(3,2) x Stab(K) charge orbit 的精确 multiplicity：
 *
 *       subspace_orbit_size * spanning_triple_count(dimension)
 *
 *   其中 spanning-triple count 为 1,7,42,168。
 *
 * r_by_cell：
 *   用 subspace basis vector 0/1/2 分别映射到 V3 bit 1/2/4
 *   得到的 concrete 16-cell representative charge assignment。
 */
typedef struct
{
    uint8_t dimension;
    uint16_t basis_coordinates[3];

    unsigned int subspace_orbit_size;
    uint64_t weight;

    N5V3 r_by_cell[N5_ODD_CHARGE_CELL_COUNT];
} N5OddChargeCase;

/*
 * 一个 fixed odd-sector K 的完整 reduced charge catalog。
 *
 * cases 使用动态数组，以适应不同输入矩阵产生的不同 quotient 大小。
 *
 * scalar_nullity：
 *   12-line scalar GF(2) system 的 nullity d。
 *
 * raw_charge_count：
 *   三个 V3 bit components 独立，因此为 8^d。
 *
 * weighted_charge_count：
 *   全部 reduced case weight 的总和；正确时必须等于 raw_charge_count。
 *
 * subspace_count：
 *   GF(2)^d 中全部 dimension<=3 subspaces 的数量。
 *
 * stabilizer_size：
 *   顶层 Stab(K) 大小。
 *
 * induced_action_count：
 *   Stab(K) 在 scalar charge coordinate space 上诱导出的不同线性作用数。
 */
typedef struct
{
    uint8_t scalar_nullity;

    unsigned int stabilizer_size;
    unsigned int induced_action_count;
    unsigned int subspace_count;

    uint64_t raw_charge_count;
    uint64_t weighted_charge_count;

    unsigned int case_count;
    N5OddChargeCase* cases;
} N5OddChargeCatalog;

/*
 * 初始化一个空 odd charge catalog。
 */
void n5_odd_charge_catalog_init(
    N5OddChargeCatalog* catalog
);

/*
 * 释放 odd charge catalog 持有的动态 case 数组，并恢复为空状态。
 */
void n5_odd_charge_catalog_free(
    N5OddChargeCatalog* catalog
);

/*
 * 对任意合法 non-all-even K 构造统一 translation-charge quotient。
 *
 * 计算链：
 *
 *   12-line XOR equations
 *       -> scalar nullspace GF(2)^d
 *       -> dim<=3 canonical scalar subspaces
 *       -> quotient by Stab(K)
 *       -> restore GL(3,2) multiplicity 1/7/42/168
 *       -> concrete V3 representative r_by_cell。
 *
 * 该函数不读取历史 reduced-count table，也不读取最终 C(K)。
 *
 * 返回：
 *   1 表示成功；
 *   0 表示输入非法、内存不足或 consistency check 失败。
 */
int n5_odd_charge_reduce(
    const N5SizeMatrix* size_matrix,
    N5OddChargeCatalog* catalog
);

/*
 * 检查一个 concrete 16-cell V3 charge assignment 是否满足统一 odd/even
 * translation-charge conservation geometry。
 */
int n5_odd_charge_assignment_is_valid(
    const N5SizeMatrix* size_matrix,
    const N5V3 r_by_cell[N5_ODD_CHARGE_CELL_COUNT]
);

#endif
