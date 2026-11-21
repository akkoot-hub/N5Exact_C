#include <stddef.h>
#include <stdlib.h>
#include <string.h>

#include "n5_odd_charge.h"
#include "n5_orbit.h"

#define N5_ODD_CHARGE_LINE_COUNT (3u * N5_V2_SIZE)
#define N5_ODD_CHARGE_INVALID_INDEX 0xFFFFFFFFu

/*
 * scalar charge system 的 nullspace basis。
 * 每个 basis_masks[i] 是一个 global 16-cell bitmask。
 */
typedef struct
{
    uint8_t active_count;
    uint8_t rank;
    uint8_t nullity;

    uint8_t active_cells[N5_ODD_CHARGE_CELL_COUNT];
    uint16_t basis_masks[N5_ODD_CHARGE_CELL_COUNT];
} N5OddScalarChargeBasis;

/*
 * 一个 dimension<=3 的 canonical scalar subspace record。
 */
typedef struct
{
    uint8_t dimension;
    uint16_t rows[3];
    uint64_t key;
} N5OddChargeSubspaceRecord;

/*
 * 动态 scalar-subspace list。
 */
typedef struct
{
    N5OddChargeSubspaceRecord* data;
    unsigned int count;
    unsigned int capacity;
} N5OddChargeSubspaceList;

/*
 * 返回 12 条 top-level line 中指定 line 的四个 cell IDs。
 *
 * line 0..3  ：rows；
 * line 4..7  ：columns；
 * line 8..11 ：xor diagonals a XOR b = w。
 */
static int n5_odd_charge_line_cells(
    unsigned int line_id,
    unsigned int cells[4]
)
{
    unsigned int index;

    if (cells == NULL ||
        line_id >= N5_ODD_CHARGE_LINE_COUNT)
    {
        return 0;
    }

    if (line_id < 4)
    {
        unsigned int row;

        row = line_id;

        for (index = 0; index < 4; ++index)
        {
            cells[index] = 4u * row + index;
        }

        return 1;
    }

    if (line_id < 8)
    {
        unsigned int column;

        column = line_id - 4u;

        for (index = 0; index < 4; ++index)
        {
            cells[index] = 4u * index + column;
        }

        return 1;
    }

    {
        unsigned int w;

        w = line_id - 8u;

        for (index = 0; index < 4; ++index)
        {
            unsigned int row;
            unsigned int column;

            row = index;
            column = row ^ w;
            cells[index] = 4u * row + column;
        }
    }

    return 1;
}

/*
 * 检查一个 concrete 16-cell V3 charge assignment。
 */
int n5_odd_charge_assignment_is_valid(
    const N5SizeMatrix* size_matrix,
    const N5V3 r_by_cell[N5_ODD_CHARGE_CELL_COUNT]
)
{
    unsigned int cell;
    unsigned int line_id;

    if (size_matrix == NULL ||
        r_by_cell == NULL ||
        !n5_size_matrix_is_legal(size_matrix))
    {
        return 0;
    }

    for (cell = 0;
         cell < N5_ODD_CHARGE_CELL_COUNT;
         ++cell)
    {
        unsigned int row;
        unsigned int column;
        unsigned int k;

        if ((unsigned int)r_by_cell[cell] >= 8u)
        {
            return 0;
        }

        row = cell >> 2u;
        column = cell & 3u;
        k = size_matrix->matrix[row][column];

        if ((k == 0u || k == 8u) &&
            r_by_cell[cell] != 0)
        {
            return 0;
        }
    }

    for (line_id = 0;
         line_id < N5_ODD_CHARGE_LINE_COUNT;
         ++line_id)
    {
        unsigned int cells[4];
        unsigned int index;
        N5V3 total;

        if (!n5_odd_charge_line_cells(
                line_id,
                cells
            ))
        {
            return 0;
        }

        total = 0;

        for (index = 0; index < 4; ++index)
        {
            total ^= r_by_cell[cells[index]];
        }

        if (total != 0)
        {
            return 0;
        }
    }

    return 1;
}

/*
 * 构造任意合法 K 的 scalar GF(2) translation-charge nullspace basis。
 *
 * 与 all-even r-flow 的几何相同；唯一特殊约束仍是 k=0/8 只能取 r=0。
 */
static int n5_odd_charge_build_scalar_basis(
    const N5SizeMatrix* size_matrix,
    N5OddScalarChargeBasis* basis
)
{
    uint16_t active_mask;
    uint16_t rows[N5_ODD_CHARGE_LINE_COUNT];
    unsigned int pivot_columns[N5_ODD_CHARGE_CELL_COUNT];
    unsigned int row_count;
    unsigned int rank;
    unsigned int cell;
    unsigned int line_id;

    if (size_matrix == NULL ||
        basis == NULL ||
        !n5_size_matrix_is_legal(size_matrix))
    {
        return 0;
    }

    memset(basis, 0, sizeof(*basis));

    active_mask = 0;

    for (cell = 0;
         cell < N5_ODD_CHARGE_CELL_COUNT;
         ++cell)
    {
        unsigned int row;
        unsigned int column;
        unsigned int k;

        row = cell >> 2u;
        column = cell & 3u;
        k = size_matrix->matrix[row][column];

        if (k != 0u && k != 8u)
        {
            basis->active_cells[basis->active_count] =
                (uint8_t)cell;
            ++basis->active_count;

            active_mask |=
                (uint16_t)(1u << cell);
        }
    }

    row_count = 0;

    for (line_id = 0;
         line_id < N5_ODD_CHARGE_LINE_COUNT;
         ++line_id)
    {
        unsigned int cells[4];
        unsigned int index;
        uint16_t row_mask;

        if (!n5_odd_charge_line_cells(line_id, cells))
        {
            return 0;
        }

        row_mask = 0;

        for (index = 0; index < 4; ++index)
        {
            row_mask |=
                (uint16_t)(1u << cells[index]);
        }

        row_mask &= active_mask;

        if (row_mask != 0)
        {
            rows[row_count] = row_mask;
            ++row_count;
        }
    }

    rank = 0;

    for (cell = 0;
         cell < N5_ODD_CHARGE_CELL_COUNT;
         ++cell)
    {
        unsigned int pivot_row;
        unsigned int search_row;

        if ((active_mask &
             (uint16_t)(1u << cell)) == 0)
        {
            continue;
        }

        pivot_row = row_count;

        for (search_row = rank;
             search_row < row_count;
             ++search_row)
        {
            if ((rows[search_row] &
                 (uint16_t)(1u << cell)) != 0)
            {
                pivot_row = search_row;
                break;
            }
        }

        if (pivot_row == row_count)
        {
            continue;
        }

        {
            uint16_t temporary;

            temporary = rows[rank];
            rows[rank] = rows[pivot_row];
            rows[pivot_row] = temporary;
        }

        for (search_row = 0;
             search_row < row_count;
             ++search_row)
        {
            if (search_row != rank &&
                (rows[search_row] &
                 (uint16_t)(1u << cell)) != 0)
            {
                rows[search_row] ^= rows[rank];
            }
        }

        pivot_columns[rank] = cell;
        ++rank;

        if (rank == row_count)
        {
            break;
        }
    }

    basis->rank = (uint8_t)rank;
    basis->nullity =
        (uint8_t)(basis->active_count - rank);

    if (basis->nullity > N5_ODD_CHARGE_MAX_NULLITY)
    {
        return 0;
    }

    {
        unsigned int basis_index;

        basis_index = 0;

        for (cell = 0;
             cell < N5_ODD_CHARGE_CELL_COUNT;
             ++cell)
        {
            int is_active;
            int is_pivot;
            unsigned int pivot_index;
            uint16_t vector;

            is_active =
                (active_mask &
                 (uint16_t)(1u << cell)) != 0;

            if (!is_active)
            {
                continue;
            }

            is_pivot = 0;

            for (pivot_index = 0;
                 pivot_index < rank;
                 ++pivot_index)
            {
                if (pivot_columns[pivot_index] == cell)
                {
                    is_pivot = 1;
                    break;
                }
            }

            if (is_pivot)
            {
                continue;
            }

            vector =
                (uint16_t)(1u << cell);

            for (pivot_index = 0;
                 pivot_index < rank;
                 ++pivot_index)
            {
                if ((rows[pivot_index] &
                     (uint16_t)(1u << cell)) != 0)
                {
                    vector |=
                        (uint16_t)(
                            1u << pivot_columns[pivot_index]
                        );
                }
            }

            basis->basis_masks[basis_index] = vector;
            ++basis_index;
        }

        if (basis_index != basis->nullity)
        {
            return 0;
        }
    }

    return 1;
}

/*
 * 计算 8^d；当前 d<=8，因此安全落在 uint64_t。
 */
static uint64_t n5_odd_charge_raw_count(
    unsigned int nullity
)
{
    uint64_t result;
    unsigned int index;

    result = UINT64_C(1);

    for (index = 0; index < nullity; ++index)
    {
        result *= UINT64_C(8);
    }

    return result;
}

/*
 * dimension d 的 ordered spanning triples 数：1,7,42,168。
 */
static uint64_t n5_odd_charge_spanning_triple_count(
    unsigned int dimension
)
{
    switch (dimension)
    {
    case 0:
        return UINT64_C(1);
    case 1:
        return UINT64_C(7);
    case 2:
        return UINT64_C(42);
    case 3:
        return UINT64_C(168);
    default:
        return UINT64_C(0);
    }
}

/*
 * 将 canonical RREF basis 编码成唯一 64-bit key。
 */
static uint64_t n5_odd_charge_subspace_key(
    unsigned int dimension,
    const uint16_t rows[3]
)
{
    return
        ((uint64_t)dimension << 48) |
        ((uint64_t)rows[0]) |
        ((uint64_t)rows[1] << 16) |
        ((uint64_t)rows[2] << 32);
}

/*
 * 对最多三个 coordinate rows 做 GF(2) RREF。
 */
static unsigned int n5_odd_charge_rref_coordinates(
    const uint16_t* input_rows,
    unsigned int input_count,
    unsigned int ambient_dimension,
    uint16_t output_rows[3]
)
{
    uint16_t rows[3];
    unsigned int row_count;
    unsigned int index;
    unsigned int rank;
    unsigned int column;

    row_count = 0;

    for (index = 0;
         index < input_count && index < 3u;
         ++index)
    {
        if (input_rows[index] != 0)
        {
            rows[row_count] = input_rows[index];
            ++row_count;
        }
    }

    rank = 0;

    for (column = 0;
         column < ambient_dimension && rank < row_count;
         ++column)
    {
        unsigned int pivot_row;
        unsigned int search_row;

        pivot_row = row_count;

        for (search_row = rank;
             search_row < row_count;
             ++search_row)
        {
            if ((rows[search_row] &
                 (uint16_t)(1u << column)) != 0)
            {
                pivot_row = search_row;
                break;
            }
        }

        if (pivot_row == row_count)
        {
            continue;
        }

        {
            uint16_t temporary;

            temporary = rows[rank];
            rows[rank] = rows[pivot_row];
            rows[pivot_row] = temporary;
        }

        for (search_row = 0;
             search_row < row_count;
             ++search_row)
        {
            if (search_row != rank &&
                (rows[search_row] &
                 (uint16_t)(1u << column)) != 0)
            {
                rows[search_row] ^= rows[rank];
            }
        }

        ++rank;
    }

    for (index = 0; index < 3u; ++index)
    {
        output_rows[index] =
            index < rank ? rows[index] : 0;
    }

    return rank;
}

/*
 * qsort comparator：subspace canonical key 升序。
 */
static int n5_odd_charge_compare_subspace(
    const void* left,
    const void* right
)
{
    const N5OddChargeSubspaceRecord* a;
    const N5OddChargeSubspaceRecord* b;

    a = (const N5OddChargeSubspaceRecord*)left;
    b = (const N5OddChargeSubspaceRecord*)right;

    if (a->key < b->key)
    {
        return -1;
    }

    if (a->key > b->key)
    {
        return 1;
    }

    return 0;
}

/*
 * qsort comparator：unsigned int 升序。
 */
static int n5_odd_charge_compare_unsigned(
    const void* left,
    const void* right
)
{
    unsigned int a;
    unsigned int b;

    a = *(const unsigned int*)left;
    b = *(const unsigned int*)right;

    if (a < b)
    {
        return -1;
    }

    if (a > b)
    {
        return 1;
    }

    return 0;
}

/*
 * 向动态 subspace list 追加一个 canonical record。
 */
static int n5_odd_charge_subspace_append(
    N5OddChargeSubspaceList* list,
    unsigned int dimension,
    const uint16_t rows[3]
)
{
    N5OddChargeSubspaceRecord* enlarged;
    N5OddChargeSubspaceRecord* record;
    unsigned int new_capacity;

    if (list == NULL || rows == NULL)
    {
        return 0;
    }

    if (list->count == list->capacity)
    {
        new_capacity =
            list->capacity == 0u
                ? 128u
                : 2u * list->capacity;

        enlarged =
            (N5OddChargeSubspaceRecord*)realloc(
                list->data,
                (size_t)new_capacity *
                    sizeof(N5OddChargeSubspaceRecord)
            );

        if (enlarged == NULL)
        {
            return 0;
        }

        list->data = enlarged;
        list->capacity = new_capacity;
    }

    record = &list->data[list->count];
    record->dimension = (uint8_t)dimension;
    record->rows[0] = rows[0];
    record->rows[1] = rows[1];
    record->rows[2] = rows[2];
    record->key =
        n5_odd_charge_subspace_key(dimension, rows);

    ++list->count;

    return 1;
}

/*
 * 对固定 pivot set 枚举所有 canonical RREF subspaces。
 */
static int n5_odd_charge_append_rref_for_pivots(
    unsigned int ambient_dimension,
    unsigned int dimension,
    const unsigned int pivots[3],
    N5OddChargeSubspaceList* list
)
{
    unsigned int free_rows[24];
    unsigned int free_columns[24];
    unsigned int free_count;
    unsigned int column;
    unsigned int row;
    unsigned int bits;
    unsigned int bit_limit;

    free_count = 0;

    for (column = 0;
         column < ambient_dimension;
         ++column)
    {
        int is_pivot;

        is_pivot = 0;

        for (row = 0; row < dimension; ++row)
        {
            if (pivots[row] == column)
            {
                is_pivot = 1;
                break;
            }
        }

        if (is_pivot)
        {
            continue;
        }

        for (row = 0; row < dimension; ++row)
        {
            if (pivots[row] < column)
            {
                free_rows[free_count] = row;
                free_columns[free_count] = column;
                ++free_count;
            }
        }
    }

    if (free_count >= 31u)
    {
        return 0;
    }

    bit_limit = 1u << free_count;

    for (bits = 0; bits < bit_limit; ++bits)
    {
        uint16_t rows[3];
        unsigned int free_index;

        rows[0] = 0;
        rows[1] = 0;
        rows[2] = 0;

        for (row = 0; row < dimension; ++row)
        {
            rows[row] =
                (uint16_t)(1u << pivots[row]);
        }

        for (free_index = 0;
             free_index < free_count;
             ++free_index)
        {
            if ((bits & (1u << free_index)) != 0)
            {
                rows[free_rows[free_index]] |=
                    (uint16_t)(
                        1u << free_columns[free_index]
                    );
            }
        }

        if (!n5_odd_charge_subspace_append(
                list,
                dimension,
                rows
            ))
        {
            return 0;
        }
    }

    return 1;
}

/*
 * 枚举 GF(2)^d 中全部 dim<=3 subspaces，每个只出现一次。
 */
static int n5_odd_charge_enumerate_subspaces(
    unsigned int ambient_dimension,
    N5OddChargeSubspaceList* list
)
{
    uint16_t zero_rows[3];
    unsigned int p0;
    unsigned int p1;
    unsigned int p2;

    if (list == NULL ||
        ambient_dimension > N5_ODD_CHARGE_MAX_NULLITY)
    {
        return 0;
    }

    list->data = NULL;
    list->count = 0;
    list->capacity = 0;

    zero_rows[0] = 0;
    zero_rows[1] = 0;
    zero_rows[2] = 0;

    if (!n5_odd_charge_subspace_append(
            list,
            0,
            zero_rows
        ))
    {
        return 0;
    }

    for (p0 = 0; p0 < ambient_dimension; ++p0)
    {
        unsigned int pivots[3];

        pivots[0] = p0;
        pivots[1] = 0;
        pivots[2] = 0;

        if (!n5_odd_charge_append_rref_for_pivots(
                ambient_dimension,
                1,
                pivots,
                list
            ))
        {
            free(list->data);
            list->data = NULL;
            return 0;
        }
    }

    for (p0 = 0; p0 < ambient_dimension; ++p0)
    {
        for (p1 = p0 + 1u;
             p1 < ambient_dimension;
             ++p1)
        {
            unsigned int pivots[3];

            pivots[0] = p0;
            pivots[1] = p1;
            pivots[2] = 0;

            if (!n5_odd_charge_append_rref_for_pivots(
                    ambient_dimension,
                    2,
                    pivots,
                    list
                ))
            {
                free(list->data);
                list->data = NULL;
                return 0;
            }
        }
    }

    for (p0 = 0; p0 < ambient_dimension; ++p0)
    {
        for (p1 = p0 + 1u;
             p1 < ambient_dimension;
             ++p1)
        {
            for (p2 = p1 + 1u;
                 p2 < ambient_dimension;
                 ++p2)
            {
                unsigned int pivots[3];

                pivots[0] = p0;
                pivots[1] = p1;
                pivots[2] = p2;

                if (!n5_odd_charge_append_rref_for_pivots(
                        ambient_dimension,
                        3,
                        pivots,
                        list
                    ))
                {
                    free(list->data);
                    list->data = NULL;
                    return 0;
                }
            }
        }
    }

    qsort(
        list->data,
        list->count,
        sizeof(N5OddChargeSubspaceRecord),
        n5_odd_charge_compare_subspace
    );

    return 1;
}

/*
 * 在 key-sorted subspace list 中 binary search。
 */
static unsigned int n5_odd_charge_find_subspace(
    const N5OddChargeSubspaceList* list,
    uint64_t key
)
{
    unsigned int low;
    unsigned int high;

    low = 0;
    high = list->count;

    while (low < high)
    {
        unsigned int middle;

        middle = low + ((high - low) >> 1u);

        if (list->data[middle].key < key)
        {
            low = middle + 1u;
        }
        else
        {
            high = middle;
        }
    }

    if (low < list->count &&
        list->data[low].key == key)
    {
        return low;
    }

    return N5_ODD_CHARGE_INVALID_INDEX;
}

/*
 * 对 16-bit scalar cell mask 应用顶层 cell permutation。
 */
static uint16_t n5_odd_charge_permute_cell_mask(
    uint16_t mask,
    const N5CellPermutation* permutation
)
{
    uint16_t output;
    unsigned int old_cell;

    output = 0;

    for (old_cell = 0;
         old_cell < N5_ODD_CHARGE_CELL_COUNT;
         ++old_cell)
    {
        if ((mask &
             (uint16_t)(1u << old_cell)) != 0)
        {
            unsigned int new_cell;

            new_cell = permutation->destination[old_cell];

            output |=
                (uint16_t)(1u << new_cell);
        }
    }

    return output;
}

/*
 * 收集固定 K 的完整顶层 stabilizer。
 */
static unsigned int n5_odd_charge_collect_stabilizer(
    const N5SizeMatrix* size_matrix,
    N5CellPermutation output[N5_TOP_GROUP_SIZE]
)
{
    N5CellPermutation group[N5_TOP_GROUP_SIZE];
    unsigned int group_count;
    unsigned int count;
    unsigned int index;

    group_count = n5_generate_top_group(group);

    count = 0;

    for (index = 0; index < group_count; ++index)
    {
        N5SizeMatrix transformed;

        n5_apply_cell_permutation_to_size_matrix(
            size_matrix,
            &group[index],
            &transformed
        );

        if (n5_size_matrix_equal(
                size_matrix,
                &transformed
            ))
        {
            output[count] = group[index];
            ++count;
        }
    }

    return count;
}

/*
 * 建立 scalar nullspace coordinate <-> global cell-mask maps。
 */
static int n5_odd_charge_build_coordinate_maps(
    const N5OddScalarChargeBasis* basis,
    uint16_t coord_to_mask[256],
    int* mask_to_coord,
    unsigned int* coord_count_out
)
{
    unsigned int coord_count;
    unsigned int coord;
    unsigned int mask_value;

    if (basis == NULL ||
        mask_to_coord == NULL ||
        coord_count_out == NULL ||
        basis->nullity > N5_ODD_CHARGE_MAX_NULLITY)
    {
        return 0;
    }

    coord_count = 1u << basis->nullity;

    for (mask_value = 0;
         mask_value < (1u << 16);
         ++mask_value)
    {
        mask_to_coord[mask_value] = -1;
    }

    for (coord = 0; coord < coord_count; ++coord)
    {
        uint16_t mask;
        unsigned int basis_index;

        mask = 0;

        for (basis_index = 0;
             basis_index < basis->nullity;
             ++basis_index)
        {
            if ((coord & (1u << basis_index)) != 0)
            {
                mask ^= basis->basis_masks[basis_index];
            }
        }

        coord_to_mask[coord] = mask;

        if (mask_to_coord[mask] >= 0)
        {
            return 0;
        }

        mask_to_coord[mask] = (int)coord;
    }

    *coord_count_out = coord_count;
    return 1;
}

/*
 * 构造 Stab(K) 在 scalar coordinate space GF(2)^d 上诱导出的
 * 不同 linear actions。
 *
 * 每个 action 保存完整 lookup table：action[coord] = transformed_coord。
 */
static int n5_odd_charge_build_stabilizer_actions(
    const N5SizeMatrix* size_matrix,
    const N5OddScalarChargeBasis* basis,
    unsigned int* stabilizer_size_out,
    uint16_t** actions_out,
    unsigned int* action_count_out,
    unsigned int* coord_count_out,
    uint16_t coord_to_mask[256]
)
{
    N5CellPermutation stabilizer[N5_TOP_GROUP_SIZE];
    int* mask_to_coord;
    unsigned int coord_count;
    unsigned int stabilizer_count;
    uint16_t* actions;
    unsigned int action_count;
    unsigned int perm_index;

    if (stabilizer_size_out == NULL ||
        actions_out == NULL ||
        action_count_out == NULL ||
        coord_count_out == NULL)
    {
        return 0;
    }

    *actions_out = NULL;

    mask_to_coord =
        (int*)malloc(
            (size_t)(1u << 16) * sizeof(int)
        );

    if (mask_to_coord == NULL)
    {
        return 0;
    }

    if (!n5_odd_charge_build_coordinate_maps(
            basis,
            coord_to_mask,
            mask_to_coord,
            &coord_count
        ))
    {
        free(mask_to_coord);
        return 0;
    }

    stabilizer_count =
        n5_odd_charge_collect_stabilizer(
            size_matrix,
            stabilizer
        );

    if (stabilizer_count == 0)
    {
        free(mask_to_coord);
        return 0;
    }

    actions =
        (uint16_t*)malloc(
            (size_t)stabilizer_count *
            (size_t)coord_count *
            sizeof(uint16_t)
        );

    if (actions == NULL)
    {
        free(mask_to_coord);
        return 0;
    }

    action_count = 0;

    for (perm_index = 0;
         perm_index < stabilizer_count;
         ++perm_index)
    {
        uint16_t temporary[256];
        unsigned int coord;
        int valid;
        int duplicate;
        unsigned int existing;

        valid = 1;

        for (coord = 0; coord < coord_count; ++coord)
        {
            uint16_t transformed_mask;
            int transformed_coord;

            transformed_mask =
                n5_odd_charge_permute_cell_mask(
                    coord_to_mask[coord],
                    &stabilizer[perm_index]
                );

            transformed_coord =
                mask_to_coord[transformed_mask];

            if (transformed_coord < 0)
            {
                valid = 0;
                break;
            }

            temporary[coord] =
                (uint16_t)transformed_coord;
        }

        if (!valid)
        {
            free(actions);
            free(mask_to_coord);
            return 0;
        }

        duplicate = 0;

        for (existing = 0;
             existing < action_count;
             ++existing)
        {
            if (memcmp(
                    &actions[existing * coord_count],
                    temporary,
                    (size_t)coord_count * sizeof(uint16_t)
                ) == 0)
            {
                duplicate = 1;
                break;
            }
        }

        if (!duplicate)
        {
            memcpy(
                &actions[action_count * coord_count],
                temporary,
                (size_t)coord_count * sizeof(uint16_t)
            );

            ++action_count;
        }
    }

    free(mask_to_coord);

    *stabilizer_size_out = stabilizer_count;
    *actions_out = actions;
    *action_count_out = action_count;
    *coord_count_out = coord_count;

    return 1;
}

/*
 * 从 canonical scalar subspace basis 恢复一个 concrete V3 charge assignment。
 */
static void n5_odd_charge_representative_from_subspace(
    const N5OddChargeSubspaceRecord* subspace,
    const uint16_t coord_to_mask[256],
    N5V3 r_by_cell[N5_ODD_CHARGE_CELL_COUNT]
)
{
    unsigned int cell;
    unsigned int component;

    for (cell = 0;
         cell < N5_ODD_CHARGE_CELL_COUNT;
         ++cell)
    {
        r_by_cell[cell] = 0;
    }

    for (component = 0;
         component < subspace->dimension;
         ++component)
    {
        uint16_t cell_mask;
        N5V3 component_value;

        cell_mask = coord_to_mask[subspace->rows[component]];
        component_value = (N5V3)(1u << component);

        for (cell = 0;
             cell < N5_ODD_CHARGE_CELL_COUNT;
             ++cell)
        {
            if ((cell_mask &
                 (uint16_t)(1u << cell)) != 0)
            {
                r_by_cell[cell] ^= component_value;
            }
        }
    }
}

/*
 * 初始化一个空 dynamic catalog。
 */
void n5_odd_charge_catalog_init(
    N5OddChargeCatalog* catalog
)
{
    if (catalog != NULL)
    {
        memset(catalog, 0, sizeof(*catalog));
    }
}

/*
 * 释放 dynamic cases。
 */
void n5_odd_charge_catalog_free(
    N5OddChargeCatalog* catalog
)
{
    if (catalog == NULL)
    {
        return;
    }

    free(catalog->cases);
    memset(catalog, 0, sizeof(*catalog));
}

/*
 * 向 dynamic case array 追加一个 case。
 */
static int n5_odd_charge_append_case(
    N5OddChargeCatalog* catalog,
    unsigned int* capacity,
    const N5OddChargeCase* source
)
{
    N5OddChargeCase* enlarged;
    unsigned int new_capacity;

    if (catalog == NULL ||
        capacity == NULL ||
        source == NULL)
    {
        return 0;
    }

    if (catalog->case_count == *capacity)
    {
        new_capacity =
            *capacity == 0u
                ? 256u
                : 2u * (*capacity);

        enlarged =
            (N5OddChargeCase*)realloc(
                catalog->cases,
                (size_t)new_capacity *
                    sizeof(N5OddChargeCase)
            );

        if (enlarged == NULL)
        {
            return 0;
        }

        catalog->cases = enlarged;
        *capacity = new_capacity;
    }

    catalog->cases[catalog->case_count] = *source;
    ++catalog->case_count;

    return 1;
}

/*
 * 构造完整 odd-sector reduced translation-charge catalog。
 */
int n5_odd_charge_reduce(
    const N5SizeMatrix* size_matrix,
    N5OddChargeCatalog* catalog
)
{
    N5OddScalarChargeBasis basis;
    uint16_t* actions;
    unsigned int stabilizer_size;
    unsigned int action_count;
    unsigned int coord_count;
    uint16_t coord_to_mask[256];

    N5OddChargeSubspaceList subspaces;
    unsigned char* visited;
    unsigned int* image_indices;
    unsigned int case_capacity;
    unsigned int seed_index;

    if (catalog == NULL ||
        size_matrix == NULL ||
        !n5_size_matrix_is_legal(size_matrix) ||
        n5_size_matrix_is_all_even(size_matrix))
    {
        return 0;
    }

    n5_odd_charge_catalog_free(catalog);

    if (!n5_odd_charge_build_scalar_basis(
            size_matrix,
            &basis
        ))
    {
        return 0;
    }

    actions = NULL;

    if (!n5_odd_charge_build_stabilizer_actions(
            size_matrix,
            &basis,
            &stabilizer_size,
            &actions,
            &action_count,
            &coord_count,
            coord_to_mask
        ))
    {
        return 0;
    }

    if (!n5_odd_charge_enumerate_subspaces(
            basis.nullity,
            &subspaces
        ))
    {
        free(actions);
        return 0;
    }

    visited =
        (unsigned char*)calloc(
            subspaces.count,
            sizeof(unsigned char)
        );

    image_indices =
        (unsigned int*)malloc(
            (size_t)action_count *
                sizeof(unsigned int)
        );

    if (visited == NULL ||
        image_indices == NULL)
    {
        free(image_indices);
        free(visited);
        free(subspaces.data);
        free(actions);
        return 0;
    }

    catalog->scalar_nullity = basis.nullity;
    catalog->stabilizer_size = stabilizer_size;
    catalog->induced_action_count = action_count;
    catalog->subspace_count = subspaces.count;
    catalog->raw_charge_count =
        n5_odd_charge_raw_count(basis.nullity);
    catalog->weighted_charge_count = 0;
    catalog->case_count = 0;
    catalog->cases = NULL;
    case_capacity = 0;

    for (seed_index = 0;
         seed_index < subspaces.count;
         ++seed_index)
    {
        const N5OddChargeSubspaceRecord* seed;
        unsigned int action_index;
        unsigned int image_count;
        unsigned int unique_count;
        unsigned int image_index;
        unsigned int canonical_index;
        N5OddChargeCase charge_case;
        uint64_t spanning_count;

        if (visited[seed_index] != 0)
        {
            continue;
        }

        seed = &subspaces.data[seed_index];
        image_count = 0;

        for (action_index = 0;
             action_index < action_count;
             ++action_index)
        {
            uint16_t transformed_input[3];
            uint16_t transformed_rows[3];
            unsigned int row;
            unsigned int transformed_rank;
            uint64_t transformed_key;
            unsigned int found_index;

            transformed_input[0] = 0;
            transformed_input[1] = 0;
            transformed_input[2] = 0;

            for (row = 0;
                 row < seed->dimension;
                 ++row)
            {
                transformed_input[row] =
                    actions[
                        action_index * coord_count +
                        seed->rows[row]
                    ];
            }

            transformed_rank =
                n5_odd_charge_rref_coordinates(
                    transformed_input,
                    seed->dimension,
                    basis.nullity,
                    transformed_rows
                );

            if (transformed_rank != seed->dimension)
            {
                n5_odd_charge_catalog_free(catalog);
                free(image_indices);
                free(visited);
                free(subspaces.data);
                free(actions);
                return 0;
            }

            transformed_key =
                n5_odd_charge_subspace_key(
                    transformed_rank,
                    transformed_rows
                );

            found_index =
                n5_odd_charge_find_subspace(
                    &subspaces,
                    transformed_key
                );

            if (found_index ==
                N5_ODD_CHARGE_INVALID_INDEX)
            {
                n5_odd_charge_catalog_free(catalog);
                free(image_indices);
                free(visited);
                free(subspaces.data);
                free(actions);
                return 0;
            }

            image_indices[image_count] = found_index;
            ++image_count;
        }

        qsort(
            image_indices,
            image_count,
            sizeof(unsigned int),
            n5_odd_charge_compare_unsigned
        );

        unique_count = 0;
        canonical_index = N5_ODD_CHARGE_INVALID_INDEX;

        for (image_index = 0;
             image_index < image_count;
             ++image_index)
        {
            unsigned int current;

            current = image_indices[image_index];

            if (image_index > 0 &&
                current == image_indices[image_index - 1])
            {
                continue;
            }

            if (canonical_index ==
                N5_ODD_CHARGE_INVALID_INDEX)
            {
                canonical_index = current;
            }

            visited[current] = 1;
            ++unique_count;
        }

        if (canonical_index ==
            N5_ODD_CHARGE_INVALID_INDEX)
        {
            n5_odd_charge_catalog_free(catalog);
            free(image_indices);
            free(visited);
            free(subspaces.data);
            free(actions);
            return 0;
        }

        memset(&charge_case, 0, sizeof(charge_case));

        charge_case.dimension =
            subspaces.data[canonical_index].dimension;

        charge_case.basis_coordinates[0] =
            subspaces.data[canonical_index].rows[0];
        charge_case.basis_coordinates[1] =
            subspaces.data[canonical_index].rows[1];
        charge_case.basis_coordinates[2] =
            subspaces.data[canonical_index].rows[2];

        charge_case.subspace_orbit_size = unique_count;

        spanning_count =
            n5_odd_charge_spanning_triple_count(
                charge_case.dimension
            );

        charge_case.weight =
            (uint64_t)unique_count * spanning_count;

        n5_odd_charge_representative_from_subspace(
            &subspaces.data[canonical_index],
            coord_to_mask,
            charge_case.r_by_cell
        );

        if (!n5_odd_charge_assignment_is_valid(
                size_matrix,
                charge_case.r_by_cell
            ) ||
            !n5_odd_charge_append_case(
                catalog,
                &case_capacity,
                &charge_case
            ))
        {
            n5_odd_charge_catalog_free(catalog);
            free(image_indices);
            free(visited);
            free(subspaces.data);
            free(actions);
            return 0;
        }

        catalog->weighted_charge_count +=
            charge_case.weight;
    }

    free(image_indices);
    free(visited);
    free(subspaces.data);
    free(actions);

    return 1;
}

