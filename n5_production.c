#include <stdlib.h>
#include <string.h>

#include "n5_production.h"
#include "n5_dense.h"
#include "n5_slicing.h"
#include "n5_parallel.h"

#define N5_PRODUCTION_GL3_SIZE ((N5_V3_SIZE - 1u) * (N5_V3_SIZE - 2u) * (N5_V3_SIZE - 4u))
#define N5_PRODUCTION_MAX_CODE_SPACE (1u << (3u * N5_PRODUCTION_MAX_SLICE_VARIABLES))

/*
 * 3-bit vectors over GF(2) 的 rank。
 */
static unsigned int n5_production_rank3(
    N5V3 a,
    N5V3 b,
    N5V3 c
)
{
    unsigned int rows[3];
    unsigned int rank;
    int bit;

    rows[0] = a;
    rows[1] = b;
    rows[2] = c;
    rank = 0;

    for (bit = 2; bit >= 0; --bit)
    {
        unsigned int pivot;
        unsigned int row;

        pivot = 3;

        for (row = rank; row < 3; ++row)
        {
            if ((rows[row] & (1u << (unsigned int)bit)) != 0)
            {
                pivot = row;
                break;
            }
        }

        if (pivot == 3)
        {
            continue;
        }

        {
            unsigned int temporary;

            temporary = rows[rank];
            rows[rank] = rows[pivot];
            rows[pivot] = temporary;
        }

        for (row = 0; row < 3; ++row)
        {
            if (row != rank &&
                (rows[row] & (1u << (unsigned int)bit)) != 0)
            {
                rows[row] ^= rows[rank];
            }
        }

        ++rank;

        if (rank == 3)
        {
            break;
        }
    }

    return rank;
}

/*
 * 由 basis images (g(1),g(2),g(4)) 生成完整 linear map V3->V3。
 */
static void n5_production_fill_linear_map(
    N5V3 image1,
    N5V3 image2,
    N5V3 image4,
    uint8_t map[8]
)
{
    unsigned int x;

    for (x = 0; x < 8; ++x)
    {
        N5V3 y;

        y = 0;

        if ((x & 1u) != 0)
        {
            y ^= image1;
        }

        if ((x & 2u) != 0)
        {
            y ^= image2;
        }

        if ((x & 4u) != 0)
        {
            y ^= image4;
        }

        map[x] = y;
    }
}

/*
 * 生成 GL(3,2) 的 168 个不同 linear value maps。
 */
static unsigned int n5_production_generate_gl3(
    uint8_t maps[N5_PRODUCTION_GL3_SIZE][8]
)
{
    unsigned int a;
    unsigned int b;
    unsigned int c;
    unsigned int count;

    count = 0;

    for (a = 1; a < 8; ++a)
    {
        for (b = 1; b < 8; ++b)
        {
            for (c = 1; c < 8; ++c)
            {
                if (n5_production_rank3(
                        (N5V3)a,
                        (N5V3)b,
                        (N5V3)c
                    ) != 3)
                {
                    continue;
                }

                if (count >= N5_PRODUCTION_GL3_SIZE)
                {
                    return 0;
                }

                n5_production_fill_linear_map(
                    (N5V3)a,
                    (N5V3)b,
                    (N5V3)c,
                    maps[count]
                );

                ++count;
            }
        }
    }

    return count;
}

/*
 * 判断 sparse network 中是否实际出现一个 variable。
 */
static int n5_production_network_contains_variable(
    const N5SparseNetwork* network,
    unsigned int variable
)
{
    unsigned int factor_index;

    for (factor_index = 0;
         factor_index < network->factor_count;
         ++factor_index)
    {
        const N5SparseFactor* factor;
        unsigned int position;

        factor = &network->factors[factor_index];

        for (position = 0; position < factor->arity; ++position)
        {
            if (factor->vars[position] == variable)
            {
                return 1;
            }
        }
    }

    return 0;
}

/*
 * 判断 variable 是否已经被选为 slice variable。
 */
static int n5_production_variable_selected(
    unsigned int variable,
    const uint8_t* selected,
    unsigned int selected_count
)
{
    unsigned int index;

    for (index = 0; index < selected_count; ++index)
    {
        if (selected[index] == variable)
        {
            return 1;
        }
    }

    return 0;
}

/*
 * 收集一个 variable 的实际 concrete domain，升序输出。
 */
static int n5_production_collect_domain(
    const N5SparseNetwork* network,
    unsigned int variable,
    uint8_t values[256],
    unsigned int* value_count
)
{
    int seen[256];
    unsigned int factor_index;
    unsigned int value;
    unsigned int count;

    if (network == NULL || values == NULL || value_count == NULL)
    {
        return 0;
    }

    for (value = 0; value < 256; ++value)
    {
        seen[value] = 0;
    }

    for (factor_index = 0;
         factor_index < network->factor_count;
         ++factor_index)
    {
        const N5SparseFactor* factor;
        int position;
        unsigned int p;
        unsigned int entry_index;

        factor = &network->factors[factor_index];
        position = -1;

        for (p = 0; p < factor->arity; ++p)
        {
            if (factor->vars[p] == variable)
            {
                position = (int)p;
                break;
            }
        }

        if (position < 0)
        {
            continue;
        }

        for (entry_index = 0;
             entry_index < factor->nnz;
             ++entry_index)
        {
            value = factor->entries[entry_index]
                .assignment[position];

            seen[value] = 1;
        }
    }

    count = 0;

    for (value = 0; value < 256; ++value)
    {
        if (seen[value])
        {
            values[count] = (uint8_t)value;
            ++count;
        }
    }

    *value_count = count;

    return count > 0;
}

/*
 * 从 selected variables 的 actual domains 计算 raw Cartesian count。
 */
static int n5_production_slice_assignment_count(
    const N5SparseNetwork* network,
    const uint8_t* variables,
    unsigned int variable_count,
    uint64_t* count
)
{
    uint64_t product;
    unsigned int index;

    if (count == NULL)
    {
        return 0;
    }

    product = 1;

    for (index = 0; index < variable_count; ++index)
    {
        uint8_t values[256];
        unsigned int value_count;

        if (!n5_production_collect_domain(
                network,
                variables[index],
                values,
                &value_count
            ))
        {
            return 0;
        }

        if (product > UINT64_MAX / (uint64_t)value_count)
        {
            return 0;
        }

        product *= (uint64_t)value_count;
    }

    *count = product;

    return 1;
}

/*
 * 只从 lambda variables 48..63 做 greedy structural slicing。
 */
int n5_production_choose_lambda_slices(
    const N5SparseNetwork* network,
    uint64_t target_peak,
    unsigned int max_slices,
    N5LambdaSlicePlan* plan
)
{
    N5SparseContractionPlan direct_plan;
    uint64_t current_peak;
    unsigned int selected_count;

    if (network == NULL || plan == NULL ||
        target_peak == 0 ||
        max_slices > N5_PRODUCTION_MAX_SLICE_VARIABLES ||
        !n5_sparse_choose_plan(network, &direct_plan))
    {
        return 0;
    }

    memset(plan, 0, sizeof(*plan));

    plan->unsliced_peak = direct_plan.estimated_peak;
    plan->target_peak = target_peak;

    current_peak = direct_plan.estimated_peak;
    selected_count = 0;

    while (current_peak > target_peak &&
           selected_count < max_slices)
    {
        int found;
        unsigned int best_variable;
        uint64_t best_peak;
        unsigned int variable;

        found = 0;
        best_variable = 0;
        best_peak = current_peak;

        for (variable = 48;
             variable < 64;
             ++variable)
        {
            uint8_t candidate[N5_PRODUCTION_MAX_SLICE_VARIABLES];
            unsigned int index;
            uint64_t candidate_peak;

            if (!n5_production_network_contains_variable(
                    network,
                    variable
                ) ||
                n5_production_variable_selected(
                    variable,
                    plan->variables,
                    selected_count
                ))
            {
                continue;
            }

            for (index = 0; index < selected_count; ++index)
            {
                candidate[index] = plan->variables[index];
            }

            candidate[selected_count] = (uint8_t)variable;

            if (!n5_slicing_estimate_peak(
                    network,
                    candidate,
                    selected_count + 1u,
                    &candidate_peak
                ))
            {
                continue;
            }

            if (!found ||
                candidate_peak < best_peak ||
                (candidate_peak == best_peak &&
                 variable < best_variable))
            {
                found = 1;
                best_variable = variable;
                best_peak = candidate_peak;
            }
        }

        if (!found || best_peak >= current_peak)
        {
            break;
        }

        plan->variables[selected_count] =
            (uint8_t)best_variable;

        ++selected_count;
        current_peak = best_peak;
    }

    plan->variable_count = (uint8_t)selected_count;
    plan->estimated_peak = current_peak;
    plan->target_reached = current_peak <= target_peak;

    if (!n5_production_slice_assignment_count(
            network,
            plan->variables,
            selected_count,
            &plan->raw_assignment_count
        ))
    {
        return 0;
    }

    return 1;
}

/*
 * 3-bit-per-variable assignment code。
 */
static unsigned int n5_production_encode_assignment(
    const uint8_t* values,
    unsigned int count
)
{
    unsigned int code;
    unsigned int index;

    code = 0;

    for (index = 0; index < count; ++index)
    {
        code |=
            ((unsigned int)values[index] & 7u) << (3u * index);
    }

    return code;
}

/*
 * assignment code -> concrete V3 values。
 */
static void n5_production_decode_assignment(
    unsigned int code,
    unsigned int count,
    uint8_t* values
)
{
    unsigned int index;

    for (index = 0; index < count; ++index)
    {
        values[index] =
            (uint8_t)((code >> (3u * index)) & 7u);
    }
}

/*
 * 判断 concrete assignment 是否属于 selected lambda domains。
 */
static int n5_production_assignment_in_domains(
    const uint8_t* values,
    const uint8_t* domains,
    const unsigned int* domain_counts,
    unsigned int count
)
{
    unsigned int index;

    for (index = 0; index < count; ++index)
    {
        unsigned int d;
        int found;

        found = 0;

        for (d = 0; d < domain_counts[index]; ++d)
        {
            if (domains[index * 8u + d] == values[index])
            {
                found = 1;
                break;
            }
        }

        if (!found)
        {
            return 0;
        }
    }

    return 1;
}

/*
 * map(lambda) 后，针对 r!=0 的 quotient lambda domain
 * canonicalize 到 min(lambda, lambda XOR r)。
 */
static N5V3 n5_production_transform_lambda(
    N5V3 lambda,
    N5V3 r,
    const uint8_t map[8]
)
{
    N5V3 image;

    image = map[lambda];

    if (r != 0)
    {
        N5V3 partner;

        partner = image ^ r;

        if (partner < image)
        {
            image = partner;
        }
    }

    return image;
}

/*
 * catalog lifecycle reset。
 */
void n5_production_free_lambda_slice_catalog(
    N5LambdaSliceCatalog* catalog
)
{
    if (catalog == NULL)
    {
        return;
    }

    free(catalog->orbits);
    catalog->orbits = NULL;
    catalog->orbit_count = 0;
    catalog->symmetry_count = 0;
    catalog->weighted_assignment_count = 0;
}

/*
 * 构造 fixed-flow GL3-stabilized lambda-slice quotient catalog。
 */
int n5_production_build_lambda_slice_catalog(
    const N5SparseNetwork* network,
    const N5V3 flow[N5_CONTRACTION_CELL_COUNT],
    const N5LambdaSlicePlan* plan,
    N5LambdaSliceCatalog* catalog
)
{
    uint8_t gl3[N5_PRODUCTION_GL3_SIZE][8];
    unsigned int gl3_count;

    uint8_t stabilizer[N5_PRODUCTION_GL3_SIZE][8];
    unsigned int stabilizer_count;

    uint8_t domains[N5_PRODUCTION_MAX_SLICE_VARIABLES][8];
    unsigned int domain_counts[N5_PRODUCTION_MAX_SLICE_VARIABLES];

    uint8_t* valid_code;
    uint8_t* visited;

    unsigned int code_space;
    unsigned int code;
    unsigned int orbit_capacity;

    if (network == NULL || flow == NULL || plan == NULL || catalog == NULL ||
        plan->variable_count > N5_PRODUCTION_MAX_SLICE_VARIABLES)
    {
        return 0;
    }

    memset(catalog, 0, sizeof(*catalog));
    catalog->plan = *plan;

    gl3_count = n5_production_generate_gl3(gl3);

    if (gl3_count != N5_PRODUCTION_GL3_SIZE)
    {
        return 0;
    }

    stabilizer_count = 0;

    {
        unsigned int map_index;

        for (map_index = 0;
             map_index < gl3_count;
             ++map_index)
        {
            unsigned int cell;
            int fixes_flow;

            fixes_flow = 1;

            for (cell = 0;
                 cell < N5_CONTRACTION_CELL_COUNT;
                 ++cell)
            {
                if (gl3[map_index][flow[cell]] != flow[cell])
                {
                    fixes_flow = 0;
                    break;
                }
            }

            if (fixes_flow)
            {
                memcpy(
                    stabilizer[stabilizer_count],
                    gl3[map_index],
                    8
                );

                ++stabilizer_count;
            }
        }
    }

    if (stabilizer_count == 0)
    {
        return 0;
    }

    catalog->symmetry_count = stabilizer_count;

    /*
     * selected vars 必须全部是 lambda vars，且 domain values 必须在 V3。
     */
    {
        unsigned int index;

        for (index = 0;
             index < plan->variable_count;
             ++index)
        {
            uint8_t all_values[256];
            unsigned int all_count;
            unsigned int d;

            if (plan->variables[index] < 48 ||
                plan->variables[index] >= 64 ||
                !n5_production_collect_domain(
                    network,
                    plan->variables[index],
                    all_values,
                    &all_count
                ) ||
                all_count > 8)
            {
                return 0;
            }

            domain_counts[index] = all_count;

            for (d = 0; d < all_count; ++d)
            {
                if (all_values[d] >= 8)
                {
                    return 0;
                }

                domains[index][d] = all_values[d];
            }
        }
    }

    code_space = 1u << (3u * plan->variable_count);

    if (code_space > N5_PRODUCTION_MAX_CODE_SPACE)
    {
        return 0;
    }

    valid_code = (uint8_t*)calloc(code_space, sizeof(uint8_t));
    visited = (uint8_t*)calloc(code_space, sizeof(uint8_t));

    if (valid_code == NULL || visited == NULL)
    {
        free(valid_code);
        free(visited);
        return 0;
    }

    /*
     * Mark actual Cartesian-domain assignments。
     */
    for (code = 0; code < code_space; ++code)
    {
        uint8_t values[N5_PRODUCTION_MAX_SLICE_VARIABLES];

        n5_production_decode_assignment(
            code,
            plan->variable_count,
            values
        );

        if (n5_production_assignment_in_domains(
                values,
                &domains[0][0],
                domain_counts,
                plan->variable_count
            ))
        {
            valid_code[code] = 1;
        }
    }

    orbit_capacity =
        plan->raw_assignment_count > 0
            ? (unsigned int)plan->raw_assignment_count
            : 1u;

    catalog->orbits = (N5LambdaSliceOrbit*)calloc(
        orbit_capacity,
        sizeof(N5LambdaSliceOrbit)
    );

    if (catalog->orbits == NULL)
    {
        free(valid_code);
        free(visited);
        return 0;
    }

    for (code = 0; code < code_space; ++code)
    {
        uint8_t values[N5_PRODUCTION_MAX_SLICE_VARIABLES];
        unsigned int orbit_code_count;
        unsigned int canonical_code;
        unsigned int map_index;

        if (!valid_code[code] || visited[code])
        {
            continue;
        }

        n5_production_decode_assignment(
            code,
            plan->variable_count,
            values
        );

        orbit_code_count = 0;
        canonical_code = code;

        for (map_index = 0;
             map_index < stabilizer_count;
             ++map_index)
        {
            uint8_t transformed[N5_PRODUCTION_MAX_SLICE_VARIABLES];
            unsigned int index;
            unsigned int transformed_code;

            for (index = 0;
                 index < plan->variable_count;
                 ++index)
            {
                unsigned int cell;
                N5V3 r;

                cell =
                    (unsigned int)plan->variables[index] - 48u;

                r = flow[cell];

                transformed[index] =
                    n5_production_transform_lambda(
                        values[index],
                        r,
                        stabilizer[map_index]
                    );
            }

            if (!n5_production_assignment_in_domains(
                    transformed,
                    &domains[0][0],
                    domain_counts,
                    plan->variable_count
                ))
            {
                n5_production_free_lambda_slice_catalog(catalog);
                free(valid_code);
                free(visited);
                return 0;
            }

            transformed_code =
                n5_production_encode_assignment(
                    transformed,
                    plan->variable_count
                );
            /* Exact group action: the image orbit is disjoint from previously
               visited orbits. Mark during enumeration to deduplicate in O(1)
               rather than an O(|G|^2) linear search per representative. */
            if (!visited[transformed_code])
            {
                visited[transformed_code] = 1u;
                ++orbit_code_count;

                if (transformed_code < canonical_code)
                {
                    canonical_code = transformed_code;
                }
            }
        }

        if (orbit_code_count == 0)
        {
            n5_production_free_lambda_slice_catalog(catalog);
            free(valid_code);
            free(visited);
            return 0;
        }

        if (catalog->orbit_count >= orbit_capacity)
        {
            n5_production_free_lambda_slice_catalog(catalog);
            free(valid_code);
            free(visited);
            return 0;
        }

        n5_production_decode_assignment(
            canonical_code,
            plan->variable_count,
            catalog->orbits[catalog->orbit_count].values
        );

        catalog->orbits[catalog->orbit_count].orbit_size =
            orbit_code_count;

        catalog->weighted_assignment_count +=
            orbit_code_count;

        ++catalog->orbit_count;
    }

    free(valid_code);
    free(visited);

    if (catalog->weighted_assignment_count !=
        plan->raw_assignment_count)
    {
        n5_production_free_lambda_slice_catalog(catalog);
        return 0;
    }

    return 1;
}

/*
 * conditioned network 是否显式为零。
 */
static int n5_production_network_is_zero(
    const N5SparseNetwork* network
)
{
    unsigned int factor_index;

    for (factor_index = 0;
         factor_index < network->factor_count;
         ++factor_index)
    {
        if (network->factors[factor_index].nnz == 0)
        {
            return 1;
        }
    }

    return 0;
}

/* ===================================================================== */
/* Exact affine lambda symmetry + shared six-prime contraction. */
/* ===================================================================== */

#define N5_PRODUCTION_AFFINE_MAX_SYMMETRIES (N5_V3_SIZE * N5_PRODUCTION_GL3_SIZE)

/*
 * 一个已经通过 exact network-invariance test 的 affine lambda action。
 */
typedef struct
{
    uint8_t linear[8];
    uint8_t translation;
} N5ProductionAffineLambdaSymmetry;

/*
 * sparse factor assignment lexicographic compare。
 *
 * SparseFactor 在 compress 后按完整 canonical assignment buffer 排序，
 * unused bytes 均为 0。
 */
static int n5_production_assignment_compare_full(
    const uint8_t* a,
    const uint8_t* b
)
{
    return memcmp(
        a,
        b,
        N5_SPARSE_MAX_ARITY
    );
}

/*
 * 在 canonical sparse factor 中 binary-search 一个 assignment。
 */
static const N5SparseEntry*
n5_production_factor_find_assignment(
    const N5SparseFactor* factor,
    const uint8_t* assignment
)
{
    unsigned int low;
    unsigned int high;

    low = 0;
    high = factor->nnz;

    while (low < high)
    {
        unsigned int middle;
        int comparison;

        middle = low + ((high - low) >> 1u);

        comparison =
            n5_production_assignment_compare_full(
                factor->entries[middle].assignment,
                assignment
            );

        if (comparison < 0)
        {
            low = middle + 1u;
        }
        else
        {
            high = middle;
        }
    }

    if (low < factor->nnz &&
        n5_production_assignment_compare_full(
            factor->entries[low].assignment,
            assignment
        ) == 0)
    {
        return &factor->entries[low];
    }

    return NULL;
}

/*
 * affine quotient action：
 *
 *   lambda -> L(lambda) XOR t
 *
 * 若 r != 0，则 lambda 是 V3/<r> 的 canonical representative，
 * 所以最后 canonicalize 到：
 *
 *   min(y, y XOR r)。
 *
 * L 必须固定 r，才能在 quotient 上定义良好。
 */
static N5V3 n5_production_transform_lambda_affine(
    N5V3 lambda,
    N5V3 r,
    const uint8_t linear[8],
    N5V3 translation
)
{
    N5V3 image;

    image =
        (N5V3)(linear[lambda] ^ translation);

    if (r != 0)
    {
        N5V3 partner;

        partner = (N5V3)(image ^ r);

        if (partner < image)
        {
            image = partner;
        }
    }

    return image;
}

/*
 * 判断一个 sparse network 是否为 pure-lambda network：
 *
 * 除 scalar factor 外，所有 variable 必须属于 48..63。
 *
 * uniform k=2 second-Walsh network 正是这种结构。
 */
static int n5_production_network_is_pure_lambda(
    const N5SparseNetwork* network
)
{
    unsigned int factor_index;

    if (network == NULL)
    {
        return 0;
    }

    for (factor_index = 0;
         factor_index < network->factor_count;
         ++factor_index)
    {
        const N5SparseFactor* factor;
        unsigned int position;

        factor = &network->factors[factor_index];

        for (position = 0;
             position < factor->arity;
             ++position)
        {
            unsigned int variable;

            variable = factor->vars[position];

            if (variable < 48u || variable >= 64u)
            {
                return 0;
            }
        }
    }

    return 1;
}

/*
 * 检查 linear map 是否逐 cell 固定当前 flow charge：
 *
 *   L(r_c)=r_c。
 */
static int n5_production_linear_fixes_flow(
    const uint8_t linear[8],
    const N5V3 flow[N5_CONTRACTION_CELL_COUNT]
)
{
    unsigned int cell;

    for (cell = 0;
         cell < N5_CONTRACTION_CELL_COUNT;
         ++cell)
    {
        if (linear[flow[cell]] != flow[cell])
        {
            return 0;
        }
    }

    return 1;
}

/*
 * 对一个 pure-lambda factor 做 exact affine invariance test。
 *
 * 对每个 nonzero entry x：
 *
 *   F(x) == F(gx)
 *
 * 必须 exact rational 完全相等。
 *
 * sparse 中不存在的 entry 按 0 解释。
 */
static int n5_production_factor_invariant_affine_lambda(
    const N5SparseFactor* factor,
    const N5V3 flow[N5_CONTRACTION_CELL_COUNT],
    const uint8_t linear[8],
    N5V3 translation
)
{
    unsigned int entry_index;

    if (factor->arity == 0)
    {
        return 1;
    }

    for (entry_index = 0;
         entry_index < factor->nnz;
         ++entry_index)
    {
        uint8_t transformed[N5_SPARSE_MAX_ARITY];
        unsigned int position;
        const N5SparseEntry* matched;

        memset(
            transformed,
            0,
            sizeof(transformed)
        );

        for (position = 0;
             position < factor->arity;
             ++position)
        {
            unsigned int variable;
            unsigned int cell;
            N5V3 value;

            variable = factor->vars[position];

            if (variable < 48u || variable >= 64u)
            {
                return 0;
            }

            cell = variable - 48u;
            value =
                factor->entries[entry_index]
                    .assignment[position];

            if (value >= 8)
            {
                return 0;
            }

            transformed[position] =
                n5_production_transform_lambda_affine(
                    value,
                    flow[cell],
                    linear,
                    translation
                );
        }

        matched =
            n5_production_factor_find_assignment(
                factor,
                transformed
            );

        if (matched == NULL ||
            !n5_small_rational_equal(
                matched->value,
                factor->entries[entry_index].value
            ))
        {
            return 0;
        }
    }

    return 1;
}

/*
 * exact whole-network affine invariance test。
 */
static int n5_production_network_invariant_affine_lambda(
    const N5SparseNetwork* network,
    const N5V3 flow[N5_CONTRACTION_CELL_COUNT],
    const uint8_t linear[8],
    N5V3 translation
)
{
    unsigned int factor_index;

    for (factor_index = 0;
         factor_index < network->factor_count;
         ++factor_index)
    {
        if (!n5_production_factor_invariant_affine_lambda(
                &network->factors[factor_index],
                flow,
                linear,
                translation
            ))
        {
            return 0;
        }
    }

    return 1;
}

/*
 * 收集 pure-lambda network 的全部 exact affine value symmetries。
 *
 * 候选数量最多：
 *
 *   |GL(3,2)| times the translation group size.
 *
 * 先用 flow stabilizer 条件过滤 linear part，
 * 再对完整 factor network 做 exact invariance test。
 */
static int n5_production_collect_affine_lambda_symmetries(
    const N5SparseNetwork* network,
    const N5V3 flow[N5_CONTRACTION_CELL_COUNT],
    N5ProductionAffineLambdaSymmetry*
        symmetries,
    unsigned int* symmetry_count
)
{
    uint8_t gl3[N5_PRODUCTION_GL3_SIZE][8];
    unsigned int gl3_count;
    unsigned int count;
    unsigned int map_index;

    if (network == NULL || flow == NULL ||
        symmetries == NULL || symmetry_count == NULL ||
        !n5_production_network_is_pure_lambda(network))
    {
        return 0;
    }

    gl3_count = n5_production_generate_gl3(gl3);

    if (gl3_count != N5_PRODUCTION_GL3_SIZE)
    {
        return 0;
    }

    count = 0;

    for (map_index = 0;
         map_index < gl3_count;
         ++map_index)
    {
        unsigned int translation;

        if (!n5_production_linear_fixes_flow(
                gl3[map_index],
                flow
            ))
        {
            continue;
        }

        for (translation = 0;
             translation < 8;
             ++translation)
        {
            if (!n5_production_network_invariant_affine_lambda(
                    network,
                    flow,
                    gl3[map_index],
                    (N5V3)translation
                ))
            {
                continue;
            }

            if (count >=
                N5_PRODUCTION_AFFINE_MAX_SYMMETRIES)
            {
                return 0;
            }

            memcpy(
                symmetries[count].linear,
                gl3[map_index],
                8
            );

            symmetries[count].translation =
                (uint8_t)translation;

            ++count;
        }
    }

    if (count == 0)
    {
        return 0;
    }

    *symmetry_count = count;
    return 1;
}

/*
 * 对 pure-lambda network 构造 exact affine slice quotient。
 */
int n5_production_build_affine_lambda_slice_catalog(
    const N5SparseNetwork* network,
    const N5V3 flow[N5_CONTRACTION_CELL_COUNT],
    const N5LambdaSlicePlan* plan,
    N5LambdaSliceCatalog* catalog
)
{
    N5ProductionAffineLambdaSymmetry*
        symmetries;
    unsigned int symmetry_count;

    uint8_t domains[N5_PRODUCTION_MAX_SLICE_VARIABLES][8];
    unsigned int domain_counts[N5_PRODUCTION_MAX_SLICE_VARIABLES];

    uint8_t* valid_code;
    uint8_t* visited;

    unsigned int code_space;
    unsigned int code;
    unsigned int orbit_capacity;
    unsigned int index;

    if (network == NULL || flow == NULL ||
        plan == NULL || catalog == NULL ||
        plan->variable_count == 0 ||
        plan->variable_count >
            N5_PRODUCTION_MAX_SLICE_VARIABLES ||
        !n5_production_network_is_pure_lambda(network))
    {
        return 0;
    }

    memset(catalog, 0, sizeof(*catalog));
    catalog->plan = *plan;

    symmetries =
        (N5ProductionAffineLambdaSymmetry*)calloc(
            N5_PRODUCTION_AFFINE_MAX_SYMMETRIES,
            sizeof(N5ProductionAffineLambdaSymmetry)
        );

    if (symmetries == NULL)
    {
        return 0;
    }

    if (!n5_production_collect_affine_lambda_symmetries(
            network,
            flow,
            symmetries,
            &symmetry_count
        ))
    {
        free(symmetries);
        return 0;
    }

    catalog->symmetry_count = symmetry_count;

    for (index = 0;
         index < plan->variable_count;
         ++index)
    {
        uint8_t values[256];
        unsigned int value_count;
        unsigned int d;

        if (!n5_production_collect_domain(
                network,
                plan->variables[index],
                values,
                &value_count
            ) ||
            value_count == 0 ||
            value_count > 8)
        {
            free(symmetries);
            return 0;
        }

        domain_counts[index] = value_count;

        for (d = 0; d < value_count; ++d)
        {
            domains[index][d] = values[d];
        }
    }

    code_space =
        1u << (3u * (unsigned int)plan->variable_count);

    valid_code = (uint8_t*)calloc(
        code_space,
        sizeof(uint8_t)
    );

    visited = (uint8_t*)calloc(
        code_space,
        sizeof(uint8_t)
    );

    if (valid_code == NULL || visited == NULL)
    {
        free(symmetries);
        free(valid_code);
        free(visited);
        return 0;
    }

    for (code = 0; code < code_space; ++code)
    {
        uint8_t values[N5_PRODUCTION_MAX_SLICE_VARIABLES];

        n5_production_decode_assignment(
            code,
            plan->variable_count,
            values
        );

        if (n5_production_assignment_in_domains(
                values,
                &domains[0][0],
                domain_counts,
                plan->variable_count
            ))
        {
            valid_code[code] = 1;
        }
    }

    orbit_capacity =
        plan->raw_assignment_count > 0
            ? (unsigned int)plan->raw_assignment_count
            : 1u;

    catalog->orbits =
        (N5LambdaSliceOrbit*)calloc(
            orbit_capacity,
            sizeof(N5LambdaSliceOrbit)
        );

    if (catalog->orbits == NULL)
    {
        free(symmetries);
        free(valid_code);
        free(visited);
        return 0;
    }

    for (code = 0; code < code_space; ++code)
    {
        uint8_t values[N5_PRODUCTION_MAX_SLICE_VARIABLES];
        unsigned int orbit_code_count;
        unsigned int canonical_code;
        unsigned int symmetry_index;

        if (!valid_code[code] || visited[code])
        {
            continue;
        }

        n5_production_decode_assignment(
            code,
            plan->variable_count,
            values
        );

        orbit_code_count = 0;
        canonical_code = code;

        for (symmetry_index = 0;
             symmetry_index < symmetry_count;
             ++symmetry_index)
        {
            uint8_t transformed[
                N5_PRODUCTION_MAX_SLICE_VARIABLES
            ];
            unsigned int transformed_code;
            unsigned int position;

            for (position = 0;
                 position < plan->variable_count;
                 ++position)
            {
                unsigned int cell;
                N5V3 r;

                cell =
                    (unsigned int)plan->variables[position] -
                    48u;

                r = flow[cell];

                transformed[position] =
                    n5_production_transform_lambda_affine(
                        values[position],
                        r,
                        symmetries[symmetry_index].linear,
                        symmetries[symmetry_index].translation
                    );
            }

            if (!n5_production_assignment_in_domains(
                    transformed,
                    &domains[0][0],
                    domain_counts,
                    plan->variable_count
                ))
            {
                n5_production_free_lambda_slice_catalog(catalog);
                free(symmetries);
                free(valid_code);
                free(visited);
                return 0;
            }

            transformed_code =
                n5_production_encode_assignment(
                    transformed,
                    plan->variable_count
                );
            /* Exact group action: the image orbit is disjoint from previously
               visited orbits. Mark during enumeration to deduplicate in O(1)
               rather than an O(|G|^2) linear search per representative. */
            if (!visited[transformed_code])
            {
                visited[transformed_code] = 1u;

                ++orbit_code_count;

                if (transformed_code < canonical_code)
                {
                    canonical_code = transformed_code;
                }
            }
        }

        if (orbit_code_count == 0)
        {
            n5_production_free_lambda_slice_catalog(catalog);
            free(symmetries);
            free(valid_code);
            free(visited);
            return 0;
        }

        if (catalog->orbit_count >= orbit_capacity)
        {
            n5_production_free_lambda_slice_catalog(catalog);
            free(symmetries);
            free(valid_code);
            free(visited);
            return 0;
        }

        n5_production_decode_assignment(
            canonical_code,
            plan->variable_count,
            catalog->orbits[catalog->orbit_count].values
        );

        catalog->orbits[catalog->orbit_count].orbit_size =
            orbit_code_count;

        catalog->weighted_assignment_count +=
            orbit_code_count;

        ++catalog->orbit_count;
    }

    free(symmetries);
    free(valid_code);
    free(visited);

    if (catalog->weighted_assignment_count !=
        plan->raw_assignment_count)
    {
        n5_production_free_lambda_slice_catalog(catalog);
        return 0;
    }

    return 1;
}

/*
 * symmetry-sliced dense multi-modulus contraction。
 */
int n5_production_contract_lambda_slice_catalog_dense_multi(
    const N5SparseNetwork* network,
    const N5LambdaSliceCatalog* catalog,
    const uint32_t* moduli,
    unsigned int modulus_count,
    N5SymmetrySlicedDenseMultiResult* result
)
{
    unsigned int orbit_index;
    unsigned int lane;
    double phase_start;
    N5PreparedCondition prepared;

    if (network == NULL || catalog == NULL ||
        moduli == NULL || result == NULL ||
        modulus_count == 0 ||
        modulus_count > N5_DENSE_MAX_MODULI ||
        catalog->weighted_assignment_count !=
            catalog->plan.raw_assignment_count)
    {
        return 0;
    }

    memset(result, 0, sizeof(*result));
    result->modulus_count = (uint8_t)modulus_count;

    for (lane = 0; lane < modulus_count; ++lane)
    {
        if (moduli[lane] <= 1)
        {
            return 0;
        }

        result->values[lane] = 0;
    }

    if (!n5_slicing_prepare_condition(network, catalog->plan.variables,
             catalog->plan.variable_count, &prepared)) return 0;

    for (orbit_index = 0;
         orbit_index < catalog->orbit_count;
         ++orbit_index)
    {
        const N5LambdaSliceOrbit* orbit;
        N5SparseNetwork conditioned;
        N5SparseContractionPlan contraction_plan;
        N5DenseMultiModularContractionResult contraction;

        orbit = &catalog->orbits[orbit_index];

        phase_start = n5_parallel_wall_time_seconds();
        if (!n5_slicing_condition_prepared(
                &prepared,
                orbit->values,
                &conditioned
            ))
        {
            return 0;
        }

        result->condition_seconds += n5_parallel_wall_time_seconds() - phase_start;

        if (n5_production_network_is_zero(&conditioned))
        {
            ++result->zero_representatives;
            n5_sparse_network_free(&conditioned);
            continue;
        }

        phase_start = n5_parallel_wall_time_seconds();
        if (!n5_sparse_choose_plan(
                &conditioned,
                &contraction_plan
            ))
        {
            n5_sparse_network_free(&conditioned);
            return 0;
        }
        result->plan_seconds += n5_parallel_wall_time_seconds() - phase_start;

        phase_start = n5_parallel_wall_time_seconds();
        if (!n5_dense_modular_contract_sparse_network_multi(
                &conditioned,
                &contraction_plan,
                moduli,
                modulus_count,
                &contraction
            ))
        {
            n5_sparse_network_free(&conditioned);
            return 0;
        }
        result->dense_seconds += n5_parallel_wall_time_seconds() - phase_start;

        for (lane = 0; lane < modulus_count; ++lane)
        {
            uint32_t weighted;

            weighted = (uint32_t)(
                ((uint64_t)contraction.values[lane] *
                 (uint64_t)(
                    orbit->orbit_size % moduli[lane]
                 )) %
                moduli[lane]
            );

            result->values[lane] = (uint32_t)(
                ((uint64_t)result->values[lane] +
                 (uint64_t)weighted) %
                moduli[lane]
            );
        }

        if (contraction.max_table_entries >
            result->max_table_entries)
        {
            result->max_table_entries =
                contraction.max_table_entries;
        }

        ++result->representative_contractions;
        n5_sparse_network_free(&conditioned);
    }

    result->represented_raw_assignments =
        catalog->weighted_assignment_count;

    return 1;
}

