#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <inttypes.h>

#include "n5_all_even.h"
#include "n5_orbit.h"
#include "n5_sparse.h"
#include "n5_contraction.h"
#include "n5_slicing.h"
#include "n5_dense.h"
#include "n5_production.h"
#include "n5_parallel.h"

/*
 * 一个 reduced flow class 的 contraction mode。
 */
typedef enum
{
    N5_ALL_EVEN_MODE_UNKNOWN = 0,
    N5_ALL_EVEN_MODE_DIRECT,
    N5_ALL_EVEN_MODE_LAMBDA_SYMMETRY_SLICED,
    N5_ALL_EVEN_MODE_GENERIC_SLICED
} N5AllEvenClassMode;

/*
 * 一个 reduced flow class 的独立 worker result。
 *
 * 每个 task 只写自己下标对应的 slot，
 * 因此 worker 间不需要 arithmetic lock。
 */
typedef struct
{
    int success;
    int rewrite_selected;
    int affine_slice_quotient;
    unsigned int slice_symmetry_count;
    N5AllEvenClassMode mode;

    uint32_t residues[N5_CRT_PRIME_COUNT];

    uint64_t raw_lambda_slices;
    uint64_t slice_representatives;
    uint64_t dense_contractions;
    uint64_t shared_dense_eliminations;
    uint64_t zero_representatives;
    uint64_t max_dense_table_entries;

    double wall_seconds;
    double network_build_seconds;
    double base_plan_seconds;
    double slice_catalog_seconds;
    double slice_plan_seconds;
    double contraction_seconds;
    double condition_seconds;
    double conditioned_plan_seconds;
    double dense_seconds;
} N5AllEvenClassResult;

/*
 * parallel flow-class worker 的只读共享输入 + 独立 result slots。
 */

typedef struct
{
    const N5SizeMatrix* size_matrix;
    const N5ReducedFlowCatalog* flow_catalog;
    N5AllEvenClassResult* class_results;
    int report_progress;
} N5AllEvenWorkerContext;

/*
 * 判断 K 是否为 pure uniform-k=2 network。
 *
 * 这是纯结构判定，不使用 orbit ID。
 *
 * 对这种 K，16 个 cell 全部可被 exact second-Walsh rewrite，
 * rewritten network 中只剩 lambda variables 48..63。
 */
static int n5_all_even_matrix_is_pure_k2(
    const N5SizeMatrix* size_matrix
)
{
    unsigned int row;
    unsigned int column;

    for (row = 0; row < 4; ++row)
    {
        for (column = 0; column < 4; ++column)
        {
            if (size_matrix->matrix[row][column] != 2u)
            {
                return 0;
            }
        }
    }

    return 1;
}


/*
 * 判断 network 是否只含 second-Walsh lambda variables。
 */
static int n5_all_even_network_is_pure_lambda(
    const N5SparseNetwork* network
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
            if (factor->vars[position] < 48u ||
                factor->vars[position] >= 64u)
            {
                return 0;
            }
        }
    }

    return 1;
}

/*
 * 对一个 network 直接执行六模 dense contraction。
 */
static int n5_all_even_contract_direct_six_prime(
    const N5SparseNetwork* network,
    const N5SparseContractionPlan* plan,
    N5AllEvenClassResult* output
)
{
    N5DenseMultiModularContractionResult dense_result;
    unsigned int prime_index;

    if (!n5_dense_modular_contract_sparse_network_multi(
            network,
            plan,
            n5_crt_prime_basis(),
            N5_CRT_PRIME_COUNT,
            &dense_result
        ))
    {
        return 0;
    }

    for (prime_index = 0;
         prime_index < N5_CRT_PRIME_COUNT;
         ++prime_index)
    {
        output->residues[prime_index] =
            dense_result.values[prime_index];
    }

    output->dense_contractions =
        N5_CRT_PRIME_COUNT;
    output->shared_dense_eliminations = 1;
    output->max_dense_table_entries =
        dense_result.max_table_entries;

    output->raw_lambda_slices = 1;
    output->slice_representatives = 1;
    output->mode = N5_ALL_EVEN_MODE_DIRECT;

    return 1;
}

/*
 * 对 second-Walsh network 使用 lambda-only symmetry slicing。
 *
 * 只有当 lambda slicer 真正达到 target peak 时才调用。
 */
static int n5_all_even_contract_lambda_symmetry_six_prime(
    const N5SparseNetwork* network,
    const N5V3 flow[N5_CONTRACTION_CELL_COUNT],
    const N5LambdaSlicePlan* slice_plan,
    N5AllEvenClassResult* output
)
{
    N5LambdaSliceCatalog slice_catalog;
    N5SymmetrySlicedDenseMultiResult dense_result;
    unsigned int prime_index;
    int use_affine;

    double phase_start;

    phase_start = n5_parallel_wall_time_seconds();
    use_affine =
        n5_all_even_network_is_pure_lambda(network);

    if (use_affine)
    {
        if (!n5_production_build_affine_lambda_slice_catalog(
                network,
                flow,
                slice_plan,
                &slice_catalog
            ))
        {
            return 0;
        }

        output->affine_slice_quotient = 1;
    }
    else
    {
        if (!n5_production_build_lambda_slice_catalog(
                network,
                flow,
                slice_plan,
                &slice_catalog
            ))
        {
            return 0;
        }
    }

    output->slice_catalog_seconds +=
        n5_parallel_wall_time_seconds() - phase_start;

    output->raw_lambda_slices =
        slice_plan->raw_assignment_count;
    output->slice_representatives =
        slice_catalog.orbit_count;
    output->slice_symmetry_count =
        slice_catalog.symmetry_count;

    phase_start = n5_parallel_wall_time_seconds();
    if (!n5_production_contract_lambda_slice_catalog_dense_multi(
            network,
            &slice_catalog,
            n5_crt_prime_basis(),
            N5_CRT_PRIME_COUNT,
            &dense_result
        ) ||
        dense_result.represented_raw_assignments !=
            slice_plan->raw_assignment_count)
    {
        n5_production_free_lambda_slice_catalog(
            &slice_catalog
        );

        return 0;
    }

    output->contraction_seconds +=
        n5_parallel_wall_time_seconds() - phase_start;
    output->condition_seconds += dense_result.condition_seconds;
    output->conditioned_plan_seconds += dense_result.plan_seconds;
    output->dense_seconds += dense_result.dense_seconds;

    for (prime_index = 0;
         prime_index < N5_CRT_PRIME_COUNT;
         ++prime_index)
    {
        output->residues[prime_index] =
            dense_result.values[prime_index];
    }

    /*
     * 一次 multi-modulus contraction 就同时得到六个 residues。
     * 因此这里统计的是实际 conditioned dense eliminations，
     * 不再把同一个 representative 机械乘以 6。
     */
    output->dense_contractions =
        (uint64_t)dense_result.representative_contractions *
        N5_CRT_PRIME_COUNT;

    output->shared_dense_eliminations =
        dense_result.representative_contractions;

    output->zero_representatives =
        dense_result.zero_representatives;

    output->max_dense_table_entries =
        dense_result.max_table_entries;

    n5_production_free_lambda_slice_catalog(
        &slice_catalog
    );

    output->mode =
        N5_ALL_EVEN_MODE_LAMBDA_SYMMETRY_SLICED;

    return 1;
}

/*
 * generic raw-assignment slicing fallback。
 *
 * 用于：
 *   - 没有 k=2 lambda variables 的 hard network；
 *   - lambda-only slicing 无法达到 target peak 的 network。
 *
 * 该路径不做 assignment symmetry quotient，
 * 但 arithmetic 仍是 exact modular dense contraction。
 */
static int n5_all_even_contract_generic_sliced_six_prime(
    const N5SparseNetwork* network,
    N5AllEvenClassResult* output
)
{
    N5SlicePlan slice_plan;
    unsigned int prime_index;
    N5DenseSlicedMultiResult dense_result;

    double phase_start;

    phase_start = n5_parallel_wall_time_seconds();
    if (!n5_slicing_choose_greedy(
            network,
            UINT64_C(32768),
            N5_SLICING_MAX_VARIABLES,
            &slice_plan
        ) ||
        !slice_plan.target_reached ||
        slice_plan.estimated_peak > UINT64_C(32768))
    {
        return 0;
    }

    output->slice_plan_seconds +=
        n5_parallel_wall_time_seconds() - phase_start;
    output->raw_lambda_slices = slice_plan.raw_assignment_count;
    output->slice_representatives = slice_plan.raw_assignment_count;

    phase_start = n5_parallel_wall_time_seconds();
    if (!n5_dense_modular_contract_sliced_multi(
            network,
            &slice_plan,
            n5_crt_prime_basis(),
            N5_CRT_PRIME_COUNT,
            &dense_result
        ) ||
        dense_result.raw_assignment_count != slice_plan.raw_assignment_count)
    {
        return 0;
    }

    output->contraction_seconds +=
        n5_parallel_wall_time_seconds() - phase_start;
    for (prime_index = 0;
         prime_index < N5_CRT_PRIME_COUNT;
         ++prime_index)
    {
        output->residues[prime_index] = dense_result.values[prime_index];
    }
    output->zero_representatives = dense_result.zero_assignment_count;
    output->shared_dense_eliminations =
        dense_result.raw_assignment_count - dense_result.zero_assignment_count;
    output->dense_contractions =
        output->shared_dense_eliminations * N5_CRT_PRIME_COUNT;
    output->max_dense_table_entries = dense_result.max_table_entries;

    output->mode = N5_ALL_EVEN_MODE_GENERIC_SLICED;
    return 1;
}

/*
 * 对一个 symmetry-reduced flow class 执行完整 six-prime fixed-flow solve。
 *
 * representation policy 与 frozen generic production architecture 一致：
 *
 *   1. 构造 ordinary canonical network；
 *   2. 构造 exact k=2 second-Walsh rewrite；
 *   3. 只有 rewrite structural peak 严格更小时才采用 rewrite；
 *   4. peak <= 32768 时直接 dense contraction；
 *   5. rewrite hard case 优先 lambda-only symmetry slicing；
 *   6. 否则使用 generic raw slicing fallback。
 *
 * callback 不做 flow orbit weighting；主线程按 class index 固定顺序聚合。
 */
static void n5_all_even_solve_class_task(
    unsigned int class_index,
    void* raw_context
)
{
    N5AllEvenWorkerContext* context;
    const N5ReducedFlowClass* flow_class;
    N5AllEvenClassResult* output;

    N5SparseNetwork ordinary_network;
    N5SparseNetwork rewritten_network;
    N5SparseContractionPlan ordinary_plan;
    N5SparseContractionPlan rewritten_plan;

    N5SparseNetwork* selected_network;
    N5SparseContractionPlan* selected_plan;

    int ordinary_built;
    int rewritten_built;
    int use_rewrite;
    int pure_k2;
    int success;

    double start_time;

    context = (N5AllEvenWorkerContext*)raw_context;
    output = &context->class_results[class_index];
    flow_class = &context->flow_catalog->classes[class_index];

    memset(output, 0, sizeof(*output));
    start_time = n5_parallel_wall_time_seconds();

    ordinary_built = 0;
    rewritten_built = 0;
    success = 0;

    pure_k2 =
        n5_all_even_matrix_is_pure_k2(
            context->size_matrix
        );

    /*
     * pure-k2 时 second-Walsh identity 已经是 exact representation；
     * 无需为了 performance comparison 再 materialize 巨大的 ordinary
     * first-Walsh network。
     *
     * 这只是 representation selection optimization，
     * 不改变任何 coefficient。
     */
    if (!pure_k2)
    {
        double phase_start = n5_parallel_wall_time_seconds();
        if (!n5_contraction_build_fixed_flow_sparse_network(
                context->size_matrix,
                flow_class->representative_flow,
                &ordinary_network
            ))
        {
            goto cleanup;
        }

        ordinary_built = 1;
        output->network_build_seconds += n5_parallel_wall_time_seconds() - phase_start;

        phase_start = n5_parallel_wall_time_seconds();
        if (!n5_sparse_choose_plan(
                &ordinary_network,
                &ordinary_plan
            ) ||
            !n5_sparse_plan_is_valid(&ordinary_plan))
        {
            goto cleanup;
        }
        output->base_plan_seconds += n5_parallel_wall_time_seconds() - phase_start;
    }

    {
        double phase_start = n5_parallel_wall_time_seconds();
    if (!n5_contraction_build_second_walsh_sparse_network(
            context->size_matrix,
            flow_class->representative_flow,
            &rewritten_network
        ))
    {
        goto cleanup;
    }

    rewritten_built = 1;
    output->network_build_seconds += n5_parallel_wall_time_seconds() - phase_start;
    phase_start = n5_parallel_wall_time_seconds();

    if (!n5_sparse_choose_plan(
            &rewritten_network,
            &rewritten_plan
        ) ||
        !n5_sparse_plan_is_valid(&rewritten_plan))
    {
        goto cleanup;
    }

    output->base_plan_seconds += n5_parallel_wall_time_seconds() - phase_start;
    }

    use_rewrite =
        pure_k2 ||
        rewritten_plan.estimated_peak <
            ordinary_plan.estimated_peak;

    output->rewrite_selected = use_rewrite;

    selected_network = use_rewrite
        ? &rewritten_network
        : &ordinary_network;

    selected_plan = use_rewrite
        ? &rewritten_plan
        : &ordinary_plan;

    if (selected_plan->estimated_peak <= UINT64_C(32768))
    {
        success = n5_all_even_contract_direct_six_prime(
            selected_network,
            selected_plan,
            output
        );
    }
    else if (use_rewrite)
    {
        N5LambdaSlicePlan lambda_plan;

        {
        double phase_start = n5_parallel_wall_time_seconds();
        if (n5_production_choose_lambda_slices(
                selected_network,
                UINT64_C(32768),
                N5_PRODUCTION_MAX_SLICE_VARIABLES,
                &lambda_plan
            ) &&
            lambda_plan.target_reached &&
            lambda_plan.estimated_peak <= UINT64_C(32768) &&
            lambda_plan.variable_count > 0)
        {
            output->slice_plan_seconds += n5_parallel_wall_time_seconds() - phase_start;
            success = n5_all_even_contract_lambda_symmetry_six_prime(
                selected_network,
                flow_class->representative_flow,
                &lambda_plan,
                output
            );
        }
        else
        {
            output->slice_plan_seconds += n5_parallel_wall_time_seconds() - phase_start;
            success = n5_all_even_contract_generic_sliced_six_prime(
                selected_network,
                output
            );
        }
        }
    }
    else
    {
        success = n5_all_even_contract_generic_sliced_six_prime(
            selected_network,
            output
        );
    }

cleanup:
    if (rewritten_built)
    {
        n5_sparse_network_free(&rewritten_network);
    }

    if (ordinary_built)
    {
        n5_sparse_network_free(&ordinary_network);
    }

    output->wall_seconds =
        n5_parallel_wall_time_seconds() - start_time;
    output->success = success;

    if (context->report_progress)
    {
        printf(
            "  [flow %3u/%u] dim=%u weight=%" PRIu64
            " raw=%" PRIu64
            " reps=%" PRIu64
            " sym=%u time=%.2f s [build %.3f plan %.3f catalog %.3f contraction_total %.3f condition %.3f cplan %.3f dense %.3f] %s\n",
            class_index + 1u,
            context->flow_catalog->class_count,
            (unsigned int)flow_class->dimension,
            flow_class->flow_orbit_size,
            output->raw_lambda_slices,
            output->slice_representatives,
            output->slice_symmetry_count,
            output->wall_seconds,
            output->network_build_seconds,
            output->base_plan_seconds + output->slice_plan_seconds,
            output->slice_catalog_seconds,
            output->contraction_seconds,
            output->condition_seconds,
            output->conditioned_plan_seconds,
            output->dense_seconds,
            success ? "complete" : "failed"
        );

        fflush(stdout);
    }
}

/*
 * 返回 all-even production 的六个 CRT primes。
 */

/*
 * generic all-even exact solver。
 */
int n5_all_even_solve_ex(
    const N5SizeMatrix* size_matrix,
    unsigned int requested_threads,
    int report_progress,
    N5AllEvenSolveResult* result
)
{
    N5ReducedFlowCatalog flows;
    N5AllEvenClassResult* class_results;
    N5AllEvenWorkerContext worker_context;
    uint32_t primes[N5_CRT_PRIME_COUNT];

    unsigned int resolved_threads;
    unsigned int class_index;
    unsigned int prime_index;

    double start_time;

    if (size_matrix == NULL || result == NULL ||
        !n5_size_matrix_is_legal(size_matrix) ||
        !n5_size_matrix_is_all_even(size_matrix))
    {
        return 0;
    }

    memset(result, 0, sizeof(*result));
    memset(&flows, 0, sizeof(flows));
    class_results = NULL;

    if (!n5_crt_get_primes(primes, N5_CRT_PRIME_COUNT))
    {
        return 0;
    }

    start_time = n5_parallel_wall_time_seconds();

    if (!n5_contraction_reduce_v3_flows_subspace(
            size_matrix,
            &flows
        ))
    {
        return 0;
    }

    if (flows.class_count == 0)
    {
        n5_contraction_reduced_flow_catalog_free(&flows);
        return 0;
    }

    class_results = (N5AllEvenClassResult*)calloc(
        (size_t)flows.class_count, sizeof(*class_results));
    if (class_results == NULL)
    {
        n5_contraction_reduced_flow_catalog_free(&flows);
        return 0;
    }

    resolved_threads = n5_parallel_resolve_thread_count(
        requested_threads,
        flows.class_count
    );

    worker_context.size_matrix = size_matrix;
    worker_context.flow_catalog = &flows;
    worker_context.class_results = class_results;
    worker_context.report_progress = report_progress;

    if (!n5_parallel_for(
            flows.class_count,
            resolved_threads,
            n5_all_even_solve_class_task,
            &worker_context
        ))
    {
        free(class_results);
        n5_contraction_reduced_flow_catalog_free(&flows);
        return 0;
    }

    for (prime_index = 0;
         prime_index < N5_CRT_PRIME_COUNT;
         ++prime_index)
    {
        result->residues[prime_index] = 0;
    }

    for (class_index = 0;
         class_index < flows.class_count;
         ++class_index)
    {
        const N5ReducedFlowClass* flow_class;
        const N5AllEvenClassResult* class_result;

        flow_class = &flows.classes[class_index];
        class_result = &class_results[class_index];

        if (!class_result->success)
        {
            free(class_results);
            n5_contraction_reduced_flow_catalog_free(&flows);
            return 0;
        }

        for (prime_index = 0;
             prime_index < N5_CRT_PRIME_COUNT;
             ++prime_index)
        {
            uint32_t modulus;
            uint32_t weighted;

            modulus = primes[prime_index];

            weighted = (uint32_t)(
                ((uint64_t)class_result->residues[prime_index] *
                 (flow_class->flow_orbit_size % modulus)) %
                modulus
            );

            result->residues[prime_index] = (uint32_t)(
                ((uint64_t)result->residues[prime_index] + weighted) %
                modulus
            );
        }

        result->total_raw_lambda_slices +=
            class_result->raw_lambda_slices;
        result->total_slice_representatives +=
            class_result->slice_representatives;
        result->total_dense_contractions +=
            class_result->dense_contractions;
        result->total_shared_dense_eliminations +=
            class_result->shared_dense_eliminations;
        result->total_zero_representatives +=
            class_result->zero_representatives;
        result->summed_task_seconds +=
            class_result->wall_seconds;

        if (class_result->max_dense_table_entries >
            result->max_dense_table_entries)
        {
            result->max_dense_table_entries =
                class_result->max_dense_table_entries;
        }

        if (class_result->rewrite_selected)
        {
            ++result->rewrite_selected_classes;
        }

        if (class_result->affine_slice_quotient)
        {
            ++result->affine_lambda_sliced_classes;
        }

        if (class_result->slice_symmetry_count >
            result->max_slice_symmetry_count)
        {
            result->max_slice_symmetry_count =
                class_result->slice_symmetry_count;
        }

        switch (class_result->mode)
        {
        case N5_ALL_EVEN_MODE_DIRECT:
            ++result->direct_flow_classes;
            break;

        case N5_ALL_EVEN_MODE_LAMBDA_SYMMETRY_SLICED:
            ++result->lambda_symmetry_sliced_classes;
            break;

        case N5_ALL_EVEN_MODE_GENERIC_SLICED:
            ++result->generic_sliced_classes;
            break;

        default:
            free(class_results);
            n5_contraction_reduced_flow_catalog_free(&flows);
            return 0;
        }
    }

    if (!n5_crt_reconstruct(
            result->residues,
            n5_crt_prime_basis(),
            N5_CRT_PRIME_COUNT,
            &result->crt
        ) ||
        !n5_biguint_to_decimal(
            &result->crt.value,
            result->decimal_value,
            N5_CRT_DECIMAL_BUFFER
        ))
    {
        free(class_results);
        n5_contraction_reduced_flow_catalog_free(&flows);
        return 0;
    }

    result->reduced_flow_count = flows.class_count;
    result->raw_flow_count = flows.raw_flow_count;
    result->weighted_flow_count = flows.weighted_flow_count;
    result->worker_threads = resolved_threads;
    result->wall_seconds =
        n5_parallel_wall_time_seconds() - start_time;

    free(class_results);
    n5_contraction_reduced_flow_catalog_free(&flows);
    return 1;
}


/*
 * 默认 quiet wrapper。
 */
int n5_all_even_solve(
    const N5SizeMatrix* size_matrix,
    unsigned int requested_threads,
    N5AllEvenSolveResult* result
)
{
    return n5_all_even_solve_ex(
        size_matrix,
        requested_threads,
        0,
        result
    );
}
