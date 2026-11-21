#include <string.h>

#include "n5_slicing.h"

/*
 * 从 sparse entries 恢复每个 variable 的完整 value domain。
 *
 * seen[v][x] = 1 表示 value x 在 variable v 的至少一个
 * factor entry 中实际出现。
 */
static int n5_slicing_collect_domains(
    const N5SparseNetwork* network,
    unsigned int domain_size[N5_SPARSE_MAX_VARIABLES],
    uint8_t domain_values
        [N5_SPARSE_MAX_VARIABLES]
        [N5_SPARSE_MAX_VALUE + 1]
)
{
    unsigned char seen
        [N5_SPARSE_MAX_VARIABLES]
        [N5_SPARSE_MAX_VALUE + 1];

    unsigned int factor_index;
    unsigned int variable;

    if (network == NULL ||
        domain_size == NULL ||
        domain_values == NULL)
    {
        return 0;
    }

    memset(seen, 0, sizeof(seen));
    memset(
        domain_size,
        0,
        sizeof(unsigned int) * N5_SPARSE_MAX_VARIABLES
    );

    for (factor_index = 0;
         factor_index < network->factor_count;
         ++factor_index)
    {
        const N5SparseFactor* factor;
        unsigned int entry_index;
        unsigned int position;

        factor =
            &network->factors[factor_index];

        for (entry_index = 0;
             entry_index < factor->nnz;
             ++entry_index)
        {
            for (position = 0;
                 position < factor->arity;
                 ++position)
            {
                unsigned int value;

                variable =
                    factor->vars[position];

                value =
                    factor->entries[entry_index]
                        .assignment[position];

                if (variable >= N5_SPARSE_MAX_VARIABLES)
                {
                    return 0;
                }

                if (!seen[variable][value])
                {
                    unsigned int slot;

                    slot =
                        domain_size[variable];

                    if (slot > N5_SPARSE_MAX_VALUE)
                    {
                        return 0;
                    }

                    seen[variable][value] = 1;

                    domain_values[variable][slot] =
                        (uint8_t)value;

                    ++domain_size[variable];
                }
            }
        }
    }

    /*
     * sparse entries 通常已经 lexicographic 排序，
     * 但这里显式把每个 value domain 做升序排序，
     * 保证 slice enumeration deterministic。
     */
    for (variable = 0;
         variable < N5_SPARSE_MAX_VARIABLES;
         ++variable)
    {
        unsigned int i;

        for (i = 1;
             i < domain_size[variable];
             ++i)
        {
            uint8_t current;
            unsigned int j;

            current =
                domain_values[variable][i];

            j = i;

            while (j > 0 &&
                   current <
                       domain_values[variable][j - 1])
            {
                domain_values[variable][j] =
                    domain_values[variable][j - 1];

                --j;
            }

            domain_values[variable][j] =
                current;
        }
    }

    return 1;
}

/*
 * 统计 64-bit scope mask 的 arity。
 */
static unsigned int n5_slicing_scope_arity(
    uint64_t scope
)
{
    unsigned int result;

    result = 0;

    while (scope != 0)
    {
        result +=
            (unsigned int)(scope & UINT64_C(1));

        scope >>= 1;
    }

    return result;
}

/*
 * 使用 variable domains 计算一个 scope 的 dense table size。
 *
 * overflow 时饱和到 UINT64_MAX。
 */
static uint64_t n5_slicing_scope_size(
    uint64_t scope,
    const unsigned int domain_size[N5_SPARSE_MAX_VARIABLES]
)
{
    uint64_t result;
    unsigned int variable;

    result = 1;

    for (variable = 0;
         variable < N5_SPARSE_MAX_VARIABLES;
         ++variable)
    {
        if ((scope &
             (UINT64_C(1) << variable)) == 0)
        {
            continue;
        }

        if (domain_size[variable] == 0)
        {
            return 0;
        }

        if (result >
            UINT64_MAX / domain_size[variable])
        {
            return UINT64_MAX;
        }

        result *=
            domain_size[variable];
    }

    return result;
}

/*
 * 判断 variable 是否包含在 selected slice list 中。
 */
static int n5_slicing_variable_selected(
    unsigned int variable,
    const uint8_t* selected,
    unsigned int selected_count
)
{
    unsigned int index;

    for (index = 0;
         index < selected_count;
         ++index)
    {
        if (selected[index] == variable)
        {
            return 1;
        }
    }

    return 0;
}

/*
 * 对已经删除 selected variables 的 scopes，
 * 使用 reference-compatible min-cost heuristic 生成 plan 和 peak。
 */
static int n5_slicing_build_remaining_plan(
    const N5SparseNetwork* network,
    const uint8_t* selected,
    unsigned int selected_count,
    N5SparseContractionPlan* plan,
    uint64_t* estimated_peak
)
{
    unsigned int domain_size[N5_SPARSE_MAX_VARIABLES];
    uint8_t domain_values
        [N5_SPARSE_MAX_VARIABLES]
        [N5_SPARSE_MAX_VALUE + 1];

    uint64_t scopes[192];
    unsigned int scope_count;

    unsigned int active_variable_count;
    unsigned int factor_index;
    unsigned int variable;
    unsigned int step;

    uint64_t peak;

    if (network == NULL ||
        plan == NULL ||
        estimated_peak == NULL ||
        selected_count > N5_SLICING_MAX_VARIABLES ||
        (selected_count > 0 && selected == NULL))
    {
        return 0;
    }

    if (!n5_slicing_collect_domains(
            network,
            domain_size,
            domain_values
        ))
    {
        return 0;
    }

    (void)domain_values;

    for (variable = 0;
         variable < N5_SPARSE_MAX_VARIABLES;
         ++variable)
    {
        if (n5_slicing_variable_selected(
                variable,
                selected,
                selected_count
            ))
        {
            if (domain_size[variable] == 0)
            {
                return 0;
            }

            domain_size[variable] = 0;
        }
    }

    scope_count = 0;

    for (factor_index = 0;
         factor_index < network->factor_count;
         ++factor_index)
    {
        const N5SparseFactor* factor;
        uint64_t scope;
        unsigned int position;

        factor =
            &network->factors[factor_index];

        scope = 0;

        for (position = 0;
             position < factor->arity;
             ++position)
        {
            variable =
                factor->vars[position];

            if (n5_slicing_variable_selected(
                    variable,
                    selected,
                    selected_count
                ))
            {
                continue;
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
    }

    active_variable_count = 0;

    for (variable = 0;
         variable < N5_SPARSE_MAX_VARIABLES;
         ++variable)
    {
        if (domain_size[variable] != 0)
        {
            ++active_variable_count;
        }
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
        unsigned int best_variable;
        uint64_t best_new_scope;

        variable_union = 0;

        for (scope_index = 0;
             scope_index < scope_count;
             ++scope_index)
        {
            variable_union |=
                scopes[scope_index];
        }

        if (variable_union == 0)
        {
            return 0;
        }

        found = 0;
        best_cost = UINT64_MAX;
        best_arity = 0;
        best_variable = 0;
        best_new_scope = 0;

        for (variable = 0;
             variable < N5_SPARSE_MAX_VARIABLES;
             ++variable)
        {
            uint64_t bit;
            uint64_t merged;
            uint64_t new_scope;
            uint64_t cost;
            unsigned int arity;

            bit =
                UINT64_C(1) << variable;

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
                    merged |=
                        scopes[scope_index];
                }
            }

            new_scope =
                merged & ~bit;

            cost =
                n5_slicing_scope_size(
                    new_scope,
                    domain_size
                );

            arity =
                n5_slicing_scope_arity(
                    new_scope
                );

            if (!found ||
                cost < best_cost ||
                (cost == best_cost &&
                 arity < best_arity) ||
                (cost == best_cost &&
                 arity == best_arity &&
                 variable < best_variable))
            {
                found = 1;
                best_cost = cost;
                best_arity = arity;
                best_variable = variable;
                best_new_scope = new_scope;
            }
        }

        if (!found)
        {
            return 0;
        }

        plan->elimination_order[step] =
            (uint8_t)best_variable;

        if (best_cost > peak)
        {
            peak = best_cost;
        }

        {
            uint64_t bit;
            unsigned int write_index;

            bit =
                UINT64_C(1) << best_variable;

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

            scope_count =
                write_index;
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

    plan->estimated_peak =
        peak;

    *estimated_peak =
        peak;

    return 1;
}

/*
 * 计算 selected slice variables 的 raw Cartesian assignment count。
 */
static int n5_slicing_assignment_count(
    const N5SparseNetwork* network,
    const uint8_t* variables,
    unsigned int variable_count,
    uint64_t* count
)
{
    unsigned int domain_size[N5_SPARSE_MAX_VARIABLES];
    uint8_t domain_values
        [N5_SPARSE_MAX_VARIABLES]
        [N5_SPARSE_MAX_VALUE + 1];

    uint64_t result;
    unsigned int index;

    if (count == NULL ||
        variable_count > N5_SLICING_MAX_VARIABLES ||
        (variable_count > 0 && variables == NULL))
    {
        return 0;
    }

    if (!n5_slicing_collect_domains(
            network,
            domain_size,
            domain_values
        ))
    {
        return 0;
    }

    result = 1;

    for (index = 0;
         index < variable_count;
         ++index)
    {
        unsigned int variable;

        variable =
            variables[index];

        if (variable >= N5_SPARSE_MAX_VARIABLES ||
            domain_size[variable] == 0)
        {
            return 0;
        }

        if (result >
            UINT64_MAX / domain_size[variable])
        {
            return 0;
        }

        result *=
            domain_size[variable];
    }

    *count = result;

    return 1;
}

/*
 * 公开 structural peak estimator。
 */
int n5_slicing_estimate_peak(
    const N5SparseNetwork* network,
    const uint8_t* variables,
    unsigned int variable_count,
    uint64_t* estimated_peak
)
{
    N5SparseContractionPlan plan;

    return
        n5_slicing_build_remaining_plan(
            network,
            variables,
            variable_count,
            &plan,
            estimated_peak
        );
}

/*
 * generic greedy slicer。
 */
int n5_slicing_choose_greedy(
    const N5SparseNetwork* network,
    uint64_t target_peak,
    unsigned int max_slices,
    N5SlicePlan* plan
)
{
    uint64_t current_peak;
    unsigned int selected_count;

    if (network == NULL ||
        plan == NULL ||
        target_peak == 0)
    {
        return 0;
    }

    if (max_slices > N5_SLICING_MAX_VARIABLES)
    {
        max_slices =
            N5_SLICING_MAX_VARIABLES;
    }

    selected_count = 0;

    if (!n5_slicing_estimate_peak(
            network,
            NULL,
            0,
            &current_peak
        ))
    {
        return 0;
    }

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

        for (variable = 0;
             variable < N5_SPARSE_MAX_VARIABLES;
             ++variable)
        {
            uint8_t candidate[N5_SLICING_MAX_VARIABLES];
            uint64_t candidate_peak;
            unsigned int index;

            if (n5_slicing_variable_selected(
                    variable,
                    plan->variables,
                    selected_count
                ))
            {
                continue;
            }

            for (index = 0;
                 index < selected_count;
                 ++index)
            {
                candidate[index] =
                    plan->variables[index];
            }

            candidate[selected_count] =
                (uint8_t)variable;

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

        if (!found ||
            best_peak >= current_peak)
        {
            break;
        }

        plan->variables[selected_count] =
            (uint8_t)best_variable;

        ++selected_count;

        /*
         * 保持 variable list 升序，
         * 使 enumeration convention deterministic。
         */
        {
            unsigned int i;

            for (i = 1;
                 i < selected_count;
                 ++i)
            {
                uint8_t current;
                unsigned int j;

                current =
                    plan->variables[i];

                j = i;

                while (j > 0 &&
                       current <
                           plan->variables[j - 1])
                {
                    plan->variables[j] =
                        plan->variables[j - 1];

                    --j;
                }

                plan->variables[j] =
                    current;
            }
        }

        current_peak =
            best_peak;
    }

    plan->variable_count =
        (uint8_t)selected_count;

    plan->estimated_peak =
        current_peak;

    plan->target_peak =
        target_peak;

    plan->target_reached =
        current_peak <= target_peak;

    if (!n5_slicing_assignment_count(
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
 * 对一个 factor 应用 fixed assignment 并删除 fixed vars。
 */
static int n5_slicing_condition_factor(
    const N5SparseFactor* factor,
    const uint8_t* fixed_variables,
    const uint8_t* fixed_values,
    unsigned int fixed_count,
    N5SparseFactor* output
)
{
    uint8_t output_vars[N5_SPARSE_MAX_ARITY];
    int fixed_value_by_position[N5_SPARSE_MAX_ARITY];

    unsigned int output_arity;
    unsigned int position;
    unsigned int entry_index;

    if (factor == NULL ||
        output == NULL ||
        fixed_count > N5_SLICING_MAX_VARIABLES ||
        (fixed_count > 0 &&
         (fixed_variables == NULL ||
          fixed_values == NULL)))
    {
        return 0;
    }

    output_arity = 0;

    for (position = 0;
         position < factor->arity;
         ++position)
    {
        unsigned int fixed_index;
        int found;

        fixed_value_by_position[position] = -1;
        found = 0;

        for (fixed_index = 0;
             fixed_index < fixed_count;
             ++fixed_index)
        {
            if (factor->vars[position] ==
                fixed_variables[fixed_index])
            {
                fixed_value_by_position[position] =
                    fixed_values[fixed_index];

                found = 1;
                break;
            }
        }

        if (!found)
        {
            output_vars[output_arity] =
                factor->vars[position];

            ++output_arity;
        }
    }

    if (!n5_sparse_factor_init(
            output,
            output_vars,
            output_arity,
            factor->nnz
        ))
    {
        return 0;
    }

    for (entry_index = 0;
         entry_index < factor->nnz;
         ++entry_index)
    {
        uint8_t assignment[N5_SPARSE_MAX_ARITY];
        unsigned int out_position;
        int compatible;

        compatible = 1;
        out_position = 0;

        for (position = 0;
             position < factor->arity;
             ++position)
        {
            if (fixed_value_by_position[position] >= 0)
            {
                if (factor->entries[entry_index]
                        .assignment[position] !=
                    (uint8_t)fixed_value_by_position[position])
                {
                    compatible = 0;
                    break;
                }
            }
            else
            {
                assignment[out_position] =
                    factor->entries[entry_index]
                        .assignment[position];

                ++out_position;
            }
        }

        if (!compatible)
        {
            continue;
        }

        if (!n5_sparse_factor_append(
                output,
                assignment,
                factor->entries[entry_index].value
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

/* Step 3: compile positions once per flow-class network and slice scope.
 * For canonical source factors the restriction is an order-preserving,
 * injective projection of the matching entries; thus compression is redundant.
 * For a noncanonical factor retain the full compress fallback.
 */
int n5_slicing_prepare_condition(
    const N5SparseNetwork* network,
    const uint8_t* fixed_variables,
    unsigned int fixed_count,
    N5PreparedCondition* compiled
)
{
    unsigned int i, f, j;
    if (!network || !compiled ||
        fixed_count > N5_SLICING_MAX_VARIABLES ||
        (fixed_count && !fixed_variables) ||
        network->factor_count > N5_SPARSE_NETWORK_MAX_FACTORS)
        return 0;
    memset(compiled, 0, sizeof(*compiled));
    compiled->source = network;
    compiled->fixed_count = fixed_count;
    for (i = 0; i < fixed_count; ++i) {
        if (fixed_variables[i] >= N5_SPARSE_MAX_VARIABLES) return 0;
        for (j = 0; j < i; ++j)
            if (fixed_variables[i] == fixed_variables[j]) return 0;
        compiled->fixed_variables[i] = fixed_variables[i];
    }
    for (f = 0; f < network->factor_count; ++f) {
        const N5SparseFactor* factor = &network->factors[f];
        if (factor->arity > N5_SPARSE_MAX_ARITY) return 0;
        for (i = 0; i < factor->arity; ++i) {
            unsigned int found = fixed_count;
            for (j = 0; j < fixed_count; ++j)
                if (factor->vars[i] == fixed_variables[j]) { found = j; break; }
            if (found != fixed_count) {
                unsigned int k = compiled->factor[f].fixed_count++;
                compiled->factor[f].fixed_positions[k] = (uint8_t)i;
                compiled->factor[f].fixed_value_indices[k] = (uint8_t)found;
            } else {
                unsigned int k = compiled->factor[f].keep_count++;
                compiled->factor[f].keep_positions[k] = (uint8_t)i;
                compiled->factor[f].output_vars[k] = factor->vars[i];
            }
        }
        compiled->factor[f].canonical_input = 1;
        for (i = 1; i < factor->nnz; ++i) {
            if (memcmp(factor->entries[i-1].assignment,
                       factor->entries[i].assignment,
                       N5_SPARSE_MAX_ARITY) >= 0) {
                compiled->factor[f].canonical_input = 0;
                break;
            }
        }
    }
    return 1;
}

/* Step11: a single compiled factor restriction.  This is also used by the
 * per-class shared factor cache; all numerical entries retain exact rational
 * values, and noncanonical source factors still use the reference compressor. */
int n5_slicing_condition_prepared_factor(
    const N5PreparedCondition* compiled,
    unsigned int factor_index,
    const uint8_t* fixed_values,
    N5SparseFactor* conditioned
)
{
    const N5SparseFactor* source;
    unsigned int kept, nfixed, i, j;
    if (!compiled || !compiled->source || !conditioned ||
        factor_index >= compiled->source->factor_count ||
        (compiled->fixed_count && !fixed_values)) return 0;
    source = &compiled->source->factors[factor_index];
    kept = compiled->factor[factor_index].keep_count;
    nfixed = compiled->factor[factor_index].fixed_count;
    if (!n5_sparse_factor_init(conditioned,
            compiled->factor[factor_index].output_vars, kept, source->nnz))
        return 0;
    for (i = 0; i < source->nnz; ++i) {
        const N5SparseEntry* old = &source->entries[i];
        N5SparseEntry* dest;
        for (j = 0; j < nfixed; ++j) {
            if (old->assignment[compiled->factor[factor_index].fixed_positions[j]] !=
                fixed_values[compiled->factor[factor_index].fixed_value_indices[j]])
                break;
        }
        if (j != nfixed) continue;
        dest = &conditioned->entries[conditioned->nnz++];
        memset(dest->assignment, 0, N5_SPARSE_MAX_ARITY);
        for (j = 0; j < kept; ++j)
            dest->assignment[j] =
                old->assignment[compiled->factor[factor_index].keep_positions[j]];
        dest->value = old->value;
    }
    if (!compiled->factor[factor_index].canonical_input &&
        !n5_sparse_factor_compress(conditioned)) {
        n5_sparse_factor_free(conditioned);
        return 0;
    }
    return 1;
}

int n5_slicing_condition_prepared(
    const N5PreparedCondition* compiled,
    const uint8_t* fixed_values,
    N5SparseNetwork* output
)
{
    unsigned int f;
    if (!compiled || !compiled->source || !output ||
        (compiled->fixed_count && !fixed_values)) return 0;
    n5_sparse_network_init(output);
    for (f = 0; f < compiled->source->factor_count; ++f) {
        N5SparseFactor conditioned;
        if (!n5_slicing_condition_prepared_factor(
                compiled, f, fixed_values, &conditioned)) {
            n5_sparse_network_free(output);
            return 0;
        }
        if (!n5_sparse_network_add_move(output, &conditioned)) {
            n5_sparse_factor_free(&conditioned);
            n5_sparse_network_free(output);
            return 0;
        }
        if (output->factors[output->factor_count - 1u].nnz == 0) return 1;
    }
    return 1;
}

/*
 * 公开 network conditioning API。
 */
int n5_slicing_condition_sparse_network(
    const N5SparseNetwork* network,
    const uint8_t* fixed_variables,
    const uint8_t* fixed_values,
    unsigned int fixed_count,
    N5SparseNetwork* output
)
{
    unsigned int factor_index;

    if (network == NULL ||
        output == NULL ||
        fixed_count > N5_SLICING_MAX_VARIABLES ||
        (fixed_count > 0 &&
         (fixed_variables == NULL ||
          fixed_values == NULL)))
    {
        return 0;
    }

    n5_sparse_network_init(output);

    for (factor_index = 0;
         factor_index < network->factor_count;
         ++factor_index)
    {
        N5SparseFactor conditioned;

        if (!n5_slicing_condition_factor(
                &network->factors[factor_index],
                fixed_variables,
                fixed_values,
                fixed_count,
                &conditioned
            ))
        {
            n5_sparse_network_free(output);
            return 0;
        }

        if (!n5_sparse_network_add_move(
                output,
                &conditioned
            ))
        {
            n5_sparse_factor_free(&conditioned);
            n5_sparse_network_free(output);
            return 0;
        }

        /*
         * Exact absorbing-zero shortcut:
         * a single conditioned factor with nnz == 0 makes the
         * entire product network identically zero.  Remaining
         * factors cannot change that result, so do not condition them.
         */
        if (output->factors[output->factor_count - 1u].nnz == 0)
        {
            return 1;
        }
    }

    return 1;
}

