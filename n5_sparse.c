#include <stdlib.h>
#include <string.h>

#include "n5_sparse.h"

/*
 * 计算 |x|，仅用于本模块的小整数 rational。
 */
static int64_t n5_sparse_abs_i64(int64_t x)
{
    return x < 0 ? -x : x;
}

/*
 * Euclid gcd。
 */
static int64_t n5_sparse_gcd_i64(int64_t a, int64_t b)
{
    a = n5_sparse_abs_i64(a);
    b = n5_sparse_abs_i64(b);

    while (b != 0)
    {
        int64_t t;

        t = a % b;
        a = b;
        b = t;
    }

    return a == 0 ? 1 : a;
}

/*
 * canonicalize rational。
 */
N5SmallRational n5_small_rational_make(
    int64_t numerator,
    int64_t denominator
)
{
    N5SmallRational result;
    int64_t gcd;

    if (denominator == 0)
    {
        result.numerator = 0;
        result.denominator = 0;
        return result;
    }

    if (numerator == 0)
    {
        result.numerator = 0;
        result.denominator = 1;
        return result;
    }

    if (denominator < 0)
    {
        numerator = -numerator;
        denominator = -denominator;
    }

    gcd = n5_sparse_gcd_i64(
        numerator,
        denominator
    );

    result.numerator = numerator / gcd;
    result.denominator = denominator / gcd;

    return result;
}

/*
 * exact rational addition。
 *
 * 本阶段 local factors 数值很小，int64_t 足够。
 */
N5SmallRational n5_small_rational_add(
    N5SmallRational a,
    N5SmallRational b
)
{
    return n5_small_rational_make(
        a.numerator * b.denominator +
        b.numerator * a.denominator,
        a.denominator * b.denominator
    );
}

/*
 * exact rational multiplication。
 */
N5SmallRational n5_small_rational_multiply(
    N5SmallRational a,
    N5SmallRational b
)
{
    int64_t g1;
    int64_t g2;

    int64_t a_num;
    int64_t a_den;
    int64_t b_num;
    int64_t b_den;

    if (a.numerator == 0 ||
        b.numerator == 0)
    {
        return n5_small_rational_make(0, 1);
    }

    a_num = a.numerator;
    a_den = a.denominator;
    b_num = b.numerator;
    b_den = b.denominator;

    /* 先交叉约分，减少中间乘积。 */
    g1 = n5_sparse_gcd_i64(a_num, b_den);
    g2 = n5_sparse_gcd_i64(b_num, a_den);

    a_num /= g1;
    b_den /= g1;

    b_num /= g2;
    a_den /= g2;

    return n5_small_rational_make(
        a_num * b_num,
        a_den * b_den
    );
}

int n5_small_rational_is_zero(
    N5SmallRational value
)
{
    return value.numerator == 0;
}

int n5_small_rational_equal(
    N5SmallRational a,
    N5SmallRational b
)
{
    return
        a.numerator == b.numerator &&
        a.denominator == b.denominator;
}

/*
 * qsort comparator。
 *
 * Stage 18 的 prototype 使用 64-bit packed assignment key，
 * 因而实际上把 factor arity 隐式限制在 8。
 *
 * Stage 19 改为直接比较完整 assignment byte buffer。
 *
 * n5_sparse_factor_append() 会把 unused bytes 统一清零，
 * 所以对完整 N5_SPARSE_MAX_ARITY bytes 做 lexicographic
 * comparison 可以同时满足：
 *
 *   1. deterministic ordering；
 *   2. duplicate assignment detection；
 *   3. arity > 8。
 */
static int n5_sparse_entry_compare(
    const void* left,
    const void* right
)
{
    const N5SparseEntry* a;
    const N5SparseEntry* b;

    a = (const N5SparseEntry*)left;
    b = (const N5SparseEntry*)right;

    return memcmp(
        a->assignment,
        b->assignment,
        N5_SPARSE_MAX_ARITY
    );
}

/*
 * 判断两个 canonical assignment buffers 是否完全相同。
 *
 * unused bytes 已经在 append 时清零，
 * 因此无需额外传入 factor arity。
 *
 * 返回：
 *   1 表示完全相同；
 *   0 表示不同。
 */
static int n5_sparse_assignment_equal(
    const uint8_t* a,
    const uint8_t* b
)
{
    return memcmp(
        a,
        b,
        N5_SPARSE_MAX_ARITY
    ) == 0;
}

/*
 * 初始化 factor。
 */
int n5_sparse_factor_init(
    N5SparseFactor* factor,
    const uint8_t* vars,
    unsigned int arity,
    unsigned int initial_capacity
)
{
    unsigned int index;

    if (factor == NULL ||
        arity > N5_SPARSE_MAX_ARITY ||
        (arity > 0 && vars == NULL))
    {
        return 0;
    }

    factor->arity = (uint8_t)arity;
    factor->nnz = 0;
    factor->capacity = 0;
    factor->entries = NULL;

    for (index = 0; index < N5_SPARSE_MAX_ARITY; ++index)
    {
        factor->vars[index] = 0;
    }

    for (index = 0; index < arity; ++index)
    {
        unsigned int j;

        if (vars[index] >= N5_SPARSE_MAX_VARIABLES)
        {
            return 0;
        }

        for (j = 0; j < index; ++j)
        {
            if (vars[j] == vars[index])
            {
                return 0;
            }
        }

        factor->vars[index] = vars[index];
    }

    if (initial_capacity == 0)
    {
        return 1;
    }

    factor->entries =
        (N5SparseEntry*)malloc(
            initial_capacity * sizeof(N5SparseEntry)
        );

    if (factor->entries == NULL)
    {
        return 0;
    }

    factor->capacity = initial_capacity;

    return 1;
}

/*
 * 释放 factor。
 */
void n5_sparse_factor_free(
    N5SparseFactor* factor
)
{
    if (factor == NULL)
    {
        return;
    }

    free(factor->entries);

    factor->entries = NULL;
    factor->nnz = 0;
    factor->capacity = 0;
    factor->arity = 0;
}

/*
 * 扩容。
 */
static int n5_sparse_factor_reserve(
    N5SparseFactor* factor,
    unsigned int required
)
{
    unsigned int new_capacity;
    N5SparseEntry* new_entries;

    if (required <= factor->capacity)
    {
        return 1;
    }

    new_capacity =
        factor->capacity == 0
            ? 16u
            : factor->capacity;

    while (new_capacity < required)
    {
        if (new_capacity > 0x7FFFFFFFu / 2u)
        {
            return 0;
        }

        new_capacity <<= 1u;
    }

    new_entries =
        (N5SparseEntry*)realloc(
            factor->entries,
            new_capacity * sizeof(N5SparseEntry)
        );

    if (new_entries == NULL)
    {
        return 0;
    }

    factor->entries = new_entries;
    factor->capacity = new_capacity;

    return 1;
}

/*
 * append，暂不合并重复 key。
 */
int n5_sparse_factor_append(
    N5SparseFactor* factor,
    const uint8_t* assignment,
    N5SmallRational value
)
{
    N5SparseEntry* entry;
    unsigned int index;

    if (factor == NULL ||
        (factor->arity > 0 && assignment == NULL))
    {
        return 0;
    }

    if (n5_small_rational_is_zero(value))
    {
        return 1;
    }

    if (!n5_sparse_factor_reserve(
            factor,
            factor->nnz + 1u
        ))
    {
        return 0;
    }

    entry =
        &factor->entries[factor->nnz];

    for (index = 0; index < N5_SPARSE_MAX_ARITY; ++index)
    {
        entry->assignment[index] = 0;
    }

    for (index = 0; index < factor->arity; ++index)
    {
        entry->assignment[index] = assignment[index];
    }

    entry->value = value;

    ++factor->nnz;

    return 1;
}

/*
 * 排序、合并重复 assignment、删除 exact zero。
 */
int n5_sparse_factor_compress(
    N5SparseFactor* factor
)
{
    unsigned int read_index;
    unsigned int write_index;

    if (factor == NULL)
    {
        return 0;
    }

    if (factor->nnz <= 1)
    {
        return 1;
    }

    /* Most cell/line factors are emitted already in strict lexicographic
       order, with unique keys. Verify this in linear time before qsort. */
    {
        int strictly_sorted = 1;
        for (read_index = 1u; read_index < factor->nnz; ++read_index)
        {
            if (n5_sparse_entry_compare(&factor->entries[read_index-1u],
                                        &factor->entries[read_index]) >= 0)
            {
                strictly_sorted = 0;
                break;
            }
        }
        if (strictly_sorted)
        {
            return 1;
        }
    }

    qsort(
        factor->entries,
        factor->nnz,
        sizeof(N5SparseEntry),
        n5_sparse_entry_compare
    );

    write_index = 0;
    read_index = 0;

    while (read_index < factor->nnz)
    {
        N5SmallRational total;

        unsigned int next;

        total =
            factor->entries[read_index].value;

        next = read_index + 1u;

        while (next < factor->nnz &&
               n5_sparse_assignment_equal(
                   factor->entries[next].assignment,
                   factor->entries[read_index].assignment
               ))
        {
            total =
                n5_small_rational_add(
                    total,
                    factor->entries[next].value
                );

            ++next;
        }

        if (!n5_small_rational_is_zero(total))
        {
            factor->entries[write_index] =
                factor->entries[read_index];

            factor->entries[write_index].value =
                total;

            ++write_index;
        }

        read_index = next;
    }

    factor->nnz = write_index;

    return 1;
}

/*
 * 找变量在 factor 中的位置。
 */
static int n5_sparse_find_var_position(
    const N5SparseFactor* factor,
    unsigned int variable
)
{
    unsigned int index;

    for (index = 0; index < factor->arity; ++index)
    {
        if (factor->vars[index] == variable)
        {
            return (int)index;
        }
    }

    return -1;
}

/*
 * exact sparse factor multiplication。
 *
 * 当前 local rewrite factor 很小，因此使用直接 sparse join。
 */
int n5_sparse_factor_multiply(
    const N5SparseFactor* a,
    const N5SparseFactor* b,
    N5SparseFactor* output
)
{
    uint8_t union_vars[N5_SPARSE_MAX_ARITY];
    unsigned int union_arity;

    int pos_a[N5_SPARSE_MAX_ARITY];
    int pos_b[N5_SPARSE_MAX_ARITY];

    unsigned int index;
    unsigned int ia;
    unsigned int ib;

    if (a == NULL ||
        b == NULL ||
        output == NULL)
    {
        return 0;
    }

    union_arity = 0;

    for (index = 0; index < a->arity; ++index)
    {
        union_vars[union_arity] =
            a->vars[index];

        ++union_arity;
    }

    for (index = 0; index < b->arity; ++index)
    {
        if (n5_sparse_find_var_position(
                a,
                b->vars[index]
            ) < 0)
        {
            if (union_arity >=
                N5_SPARSE_MAX_ARITY)
            {
                return 0;
            }

            union_vars[union_arity] =
                b->vars[index];

            ++union_arity;
        }
    }

    if (!n5_sparse_factor_init(
            output,
            union_vars,
            union_arity,
            16
        ))
    {
        return 0;
    }

    for (index = 0; index < union_arity; ++index)
    {
        pos_a[index] =
            n5_sparse_find_var_position(
                a,
                union_vars[index]
            );

        pos_b[index] =
            n5_sparse_find_var_position(
                b,
                union_vars[index]
            );
    }

    for (ia = 0; ia < a->nnz; ++ia)
    {
        for (ib = 0; ib < b->nnz; ++ib)
        {
            int compatible;
            uint8_t assignment[N5_SPARSE_MAX_ARITY];

            compatible = 1;

            for (index = 0; index < union_arity; ++index)
            {
                if (pos_a[index] >= 0 &&
                    pos_b[index] >= 0 &&
                    a->entries[ia].assignment[pos_a[index]] !=
                    b->entries[ib].assignment[pos_b[index]])
                {
                    compatible = 0;
                    break;
                }
            }

            if (!compatible)
            {
                continue;
            }

            for (index = 0; index < union_arity; ++index)
            {
                if (pos_a[index] >= 0)
                {
                    assignment[index] =
                        a->entries[ia]
                            .assignment[pos_a[index]];
                }
                else
                {
                    assignment[index] =
                        b->entries[ib]
                            .assignment[pos_b[index]];
                }
            }

            if (!n5_sparse_factor_append(
                    output,
                    assignment,
                    n5_small_rational_multiply(
                        a->entries[ia].value,
                        b->entries[ib].value
                    )
                ))
            {
                n5_sparse_factor_free(output);
                return 0;
            }
        }
    }

    if (!n5_sparse_factor_compress(output))
    {
        n5_sparse_factor_free(output);
        return 0;
    }

    return 1;
}

/*
 * exact sum-out。
 */
int n5_sparse_factor_sum_out(
    const N5SparseFactor* input,
    unsigned int variable,
    N5SparseFactor* output
)
{
    int eliminated_position;
    uint8_t vars[N5_SPARSE_MAX_ARITY];
    unsigned int output_arity;

    unsigned int input_index;
    unsigned int position;

    if (input == NULL ||
        output == NULL)
    {
        return 0;
    }

    eliminated_position =
        n5_sparse_find_var_position(
            input,
            variable
        );

    if (eliminated_position < 0)
    {
        return 0;
    }

    output_arity = 0;

    for (position = 0;
         position < input->arity;
         ++position)
    {
        if ((int)position ==
            eliminated_position)
        {
            continue;
        }

        vars[output_arity] =
            input->vars[position];

        ++output_arity;
    }

    if (!n5_sparse_factor_init(
            output,
            vars,
            output_arity,
            input->nnz
        ))
    {
        return 0;
    }

    for (input_index = 0;
         input_index < input->nnz;
         ++input_index)
    {
        uint8_t assignment[N5_SPARSE_MAX_ARITY];
        unsigned int out_position;

        out_position = 0;

        for (position = 0;
             position < input->arity;
             ++position)
        {
            if ((int)position ==
                eliminated_position)
            {
                continue;
            }

            assignment[out_position] =
                input->entries[input_index]
                    .assignment[position];

            ++out_position;
        }

        if (!n5_sparse_factor_append(
                output,
                assignment,
                input->entries[input_index].value
            ))
        {
            n5_sparse_factor_free(output);
            return 0;
        }
    }

    if (!n5_sparse_factor_compress(output))
    {
        n5_sparse_factor_free(output);
        return 0;
    }

    return 1;
}

/*
 * scalar factor。
 */
int n5_sparse_factor_make_scalar(
    N5SmallRational value,
    N5SparseFactor* factor
)
{
    uint8_t dummy_assignment[1];

    if (factor == NULL)
    {
        return 0;
    }

    if (!n5_sparse_factor_init(
            factor,
            NULL,
            0,
            1
        ))
    {
        return 0;
    }

    dummy_assignment[0] = 0;

    if (!n5_small_rational_is_zero(value))
    {
        if (!n5_sparse_factor_append(
                factor,
                dummy_assignment,
                value
            ))
        {
            n5_sparse_factor_free(factor);
            return 0;
        }
    }

    return 1;
}

/*
 * clone data，但替换 factor variable IDs。
 * assignment local values 保持不变。
 */
int n5_sparse_factor_clone_with_vars(
    const N5SparseFactor* source,
    const uint8_t* vars,
    N5SparseFactor* output
)
{
    unsigned int index;

    if (source == NULL ||
        output == NULL)
    {
        return 0;
    }

    if (!n5_sparse_factor_init(
            output,
            vars,
            source->arity,
            source->nnz
        ))
    {
        return 0;
    }

    for (index = 0; index < source->nnz; ++index)
    {
        if (!n5_sparse_factor_append(
                output,
                source->entries[index].assignment,
                source->entries[index].value
            ))
        {
            n5_sparse_factor_free(output);
            return 0;
        }
    }

    return 1;
}

/*
 * network init/free。
 */
void n5_sparse_network_init(
    N5SparseNetwork* network
)
{
    unsigned int index;

    if (network == NULL)
    {
        return;
    }

    network->factor_count = 0;

    for (index = 0;
         index < N5_SPARSE_NETWORK_MAX_FACTORS;
         ++index)
    {
        network->factors[index].arity = 0;
        network->factors[index].nnz = 0;
        network->factors[index].capacity = 0;
        network->factors[index].entries = NULL;
    }
}

void n5_sparse_network_free(
    N5SparseNetwork* network
)
{
    unsigned int index;

    if (network == NULL)
    {
        return;
    }

    for (index = 0;
         index < network->factor_count;
         ++index)
    {
        n5_sparse_factor_free(
            &network->factors[index]
        );
    }

    network->factor_count = 0;
}

/*
 * move factor ownership 到 network。
 */
int n5_sparse_network_add_move(
    N5SparseNetwork* network,
    N5SparseFactor* factor
)
{
    if (network == NULL ||
        factor == NULL ||
        network->factor_count >=
            N5_SPARSE_NETWORK_MAX_FACTORS)
    {
        return 0;
    }

    network->factors[network->factor_count] =
        *factor;

    ++network->factor_count;

    factor->arity = 0;
    factor->nnz = 0;
    factor->capacity = 0;
    factor->entries = NULL;

    return 1;
}

unsigned int n5_sparse_network_primitive_nnz(
    const N5SparseNetwork* network
)
{
    unsigned int total;
    unsigned int index;

    if (network == NULL)
    {
        return 0;
    }

    total = 0;

    for (index = 0;
         index < network->factor_count;
         ++index)
    {
        total +=
            network->factors[index].nnz;
    }

    return total;
}

/*
 * scope arity。
 */
static unsigned int n5_sparse_scope_arity(
    uint64_t scope
)
{
    unsigned int count;

    count = 0;

    while (scope != 0)
    {
        count +=
            (unsigned int)(scope & UINT64_C(1));

        scope >>= 1;
    }

    return count;
}

/*
 * dense upper bound，overflow 时饱和。
 */
static uint64_t n5_sparse_scope_size(
    uint64_t scope,
    const unsigned int domain_size[N5_SPARSE_MAX_VARIABLES]
)
{
    uint64_t size;
    unsigned int var;

    size = 1;

    for (var = 0;
         var < N5_SPARSE_MAX_VARIABLES;
         ++var)
    {
        if ((scope &
             (UINT64_C(1) << var)) == 0)
        {
            continue;
        }

        if (domain_size[var] == 0)
        {
            return 0;
        }

        if (size >
            UINT64_MAX / domain_size[var])
        {
            return UINT64_MAX;
        }

        size *= domain_size[var];
    }

    return size;
}

/*
 * 从 actual sparse entries 推导 variable domains，
 * 然后执行与 reference planner 相同的 deterministic heuristic。
 */
int n5_sparse_choose_plan(
    const N5SparseNetwork* network,
    N5SparseContractionPlan* plan
)
{
    unsigned char seen
        [N5_SPARSE_MAX_VARIABLES]
        [N5_SPARSE_MAX_VALUE + 1];

    unsigned int domain_size[N5_SPARSE_MAX_VARIABLES];

    uint64_t scopes[192];
    unsigned int scope_count;

    unsigned int active_variable_count;
    unsigned int factor_index;
    unsigned int var;

    uint64_t peak;
    unsigned int step;

    if (network == NULL ||
        plan == NULL)
    {
        return 0;
    }

    memset(seen, 0, sizeof(seen));
    memset(domain_size, 0, sizeof(domain_size));

    scope_count = 0;

    for (factor_index = 0;
         factor_index < network->factor_count;
         ++factor_index)
    {
        const N5SparseFactor* factor;
        uint64_t scope;
        unsigned int position;
        unsigned int entry_index;

        factor =
            &network->factors[factor_index];

        scope = 0;

        for (position = 0;
             position < factor->arity;
             ++position)
        {
            unsigned int variable;

            variable =
                factor->vars[position];

            if (variable >=
                N5_SPARSE_MAX_VARIABLES)
            {
                return 0;
            }

            scope |=
                UINT64_C(1) << variable;
        }

        if (scope != 0)
        {
            if (scope_count >= 192)
            {
                return 0;
            }

            scopes[scope_count] = scope;
            ++scope_count;
        }

        for (entry_index = 0;
             entry_index < factor->nnz;
             ++entry_index)
        {
            for (position = 0;
                 position < factor->arity;
                 ++position)
            {
                unsigned int variable;
                unsigned int value;

                variable =
                    factor->vars[position];

                value =
                    factor->entries[entry_index]
                        .assignment[position];

                if (!seen[variable][value])
                {
                    seen[variable][value] = 1;
                    ++domain_size[variable];
                }
            }
        }
    }

    active_variable_count = 0;

    for (var = 0;
         var < N5_SPARSE_MAX_VARIABLES;
         ++var)
    {
        if (domain_size[var] != 0)
        {
            ++active_variable_count;
        }
    }

    if (active_variable_count >
        N5_SPARSE_MAX_VARIABLES)
    {
        return 0;
    }

    peak = 1;

    for (step = 0;
         step < active_variable_count;
         ++step)
    {
        uint64_t variable_union;
        unsigned int scope_index;

        int found;
        uint64_t best_cost;
        unsigned int best_arity;
        unsigned int best_var;
        uint64_t best_new_scope;

        variable_union = 0;

        for (scope_index = 0;
             scope_index < scope_count;
             ++scope_index)
        {
            variable_union |= scopes[scope_index];
        }

        if (variable_union == 0)
        {
            return 0;
        }

        found = 0;
        best_cost = UINT64_MAX;
        best_arity = 0;
        best_var = 0;
        best_new_scope = 0;

        for (var = 0;
             var < N5_SPARSE_MAX_VARIABLES;
             ++var)
        {
            uint64_t bit;
            uint64_t merged;
            uint64_t new_scope;

            uint64_t cost;
            unsigned int arity;

            bit = UINT64_C(1) << var;

            if ((variable_union & bit) == 0)
            {
                continue;
            }

            merged = 0;

            for (scope_index = 0;
                 scope_index < scope_count;
                 ++scope_index)
            {
                if ((scopes[scope_index] & bit) != 0)
                {
                    merged |= scopes[scope_index];
                }
            }

            new_scope = merged & ~bit;

            cost =
                n5_sparse_scope_size(
                    new_scope,
                    domain_size
                );

            arity =
                n5_sparse_scope_arity(
                    new_scope
                );

            if (!found ||
                cost < best_cost ||
                (cost == best_cost &&
                 arity < best_arity) ||
                (cost == best_cost &&
                 arity == best_arity &&
                 var < best_var))
            {
                found = 1;
                best_cost = cost;
                best_arity = arity;
                best_var = var;
                best_new_scope = new_scope;
            }
        }

        if (!found)
        {
            return 0;
        }

        plan->elimination_order[step] =
            (uint8_t)best_var;

        if (best_cost > peak)
        {
            peak = best_cost;
        }

        {
            uint64_t bit;
            unsigned int write_index;

            bit = UINT64_C(1) << best_var;
            write_index = 0;

            for (scope_index = 0;
                 scope_index < scope_count;
                 ++scope_index)
            {
                if ((scopes[scope_index] & bit) == 0)
                {
                    scopes[write_index] =
                        scopes[scope_index];

                    ++write_index;
                }
            }

            scope_count = write_index;
        }

        if (best_new_scope != 0)
        {
            if (scope_count >= 192)
            {
                return 0;
            }

            scopes[scope_count] =
                best_new_scope;

            ++scope_count;
        }
    }

    plan->variable_count =
        (uint8_t)active_variable_count;

    plan->estimated_peak = peak;

    return 1;
}

int n5_sparse_plan_is_valid(
    const N5SparseContractionPlan* plan
)
{
    unsigned int seen[N5_SPARSE_MAX_VARIABLES];
    unsigned int index;

    if (plan == NULL ||
        plan->variable_count >
            N5_SPARSE_MAX_VARIABLES)
    {
        return 0;
    }

    for (index = 0;
         index < N5_SPARSE_MAX_VARIABLES;
         ++index)
    {
        seen[index] = 0;
    }

    for (index = 0;
         index < plan->variable_count;
         ++index)
    {
        unsigned int variable;

        variable =
            plan->elimination_order[index];

        if (variable >=
            N5_SPARSE_MAX_VARIABLES ||
            seen[variable] != 0)
        {
            return 0;
        }

        seen[variable] = 1;
    }

    return 1;
}

/*
 * brute-force evaluator 使用的递归上下文。
 */
typedef struct
{
    const N5SparseNetwork* network;

    unsigned int variable_count;
    uint8_t variables[N5_SPARSE_MAX_VARIABLES];

    unsigned int domain_count[N5_SPARSE_MAX_VARIABLES];
    uint8_t domain_values
        [N5_SPARSE_MAX_VARIABLES]
        [N5_SPARSE_MAX_VALUE + 1];

    uint8_t global_assignment[N5_SPARSE_MAX_VARIABLES];

    N5SmallRational total;
} N5SparseBruteForceContext;

