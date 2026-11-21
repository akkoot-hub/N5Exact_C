#include "n5_odd_production.h"

#include <stdlib.h>
#include <string.h>

#include "n5_odd_charge.h"
#include "n5_parallel.h"

/*
 * 一个 reduced charge class 的独立 worker output。
 *
 * 每个 worker 只写 class_results[class_index]，因此 contraction 阶段不需要锁。
 */
typedef struct
{
    int success;
    N5OddAdaptiveClassResult adaptive;
} N5OddProductionClassResult;

/*
 * parallel-for 的共享 context。
 * size_matrix / library / charge_catalog 只读；class_results 分离写。
 */
typedef struct
{
    const N5SizeMatrix* size_matrix;
    const N5OddKernelLibrary* library;
    const N5OddChargeCatalog* charge_catalog;

    uint64_t target_peak;
    unsigned int max_slices;

    N5OddProductionClassResult* class_results;
} N5OddProductionParallelContext;

/*
 * 一个 reduced charge class 的 generic adaptive production worker。
 */
static void n5_odd_production_adaptive_worker(
    unsigned int class_index,
    void* raw_context
)
{
    N5OddProductionParallelContext* context;
    N5OddProductionClassResult* output;
    const N5OddChargeCase* charge_case;

    context = (N5OddProductionParallelContext*)raw_context;

    if (context == NULL ||
        context->size_matrix == NULL ||
        context->library == NULL ||
        context->charge_catalog == NULL ||
        context->class_results == NULL ||
        class_index >= context->charge_catalog->case_count)
    {
        return;
    }

    output = &context->class_results[class_index];
    memset(output, 0, sizeof(*output));

    charge_case = &context->charge_catalog->cases[class_index];

    if (!n5_odd_adaptive_contract_fixed_charge(
            context->size_matrix,
            charge_case->r_by_cell,
            context->library,
            context->target_peak,
            context->max_slices,
            &output->adaptive
        ))
    {
        return;
    }

    output->success = output->adaptive.success;
}


int n5_odd_production_solve_adaptive_modular(
    const N5SizeMatrix* size_matrix,
    const N5OddKernelLibrary* library,
    unsigned int requested_threads,
    uint64_t target_peak,
    unsigned int max_slices,
    N5OddProductionResult* result
)
{
    N5OddChargeCatalog charge_catalog;
    N5OddProductionClassResult* class_results;
    N5OddProductionParallelContext context;
    uint32_t primes[N5_CRT_PRIME_COUNT];
    unsigned int class_index;
    unsigned int prime_index;
    double start_time;

    if (size_matrix == NULL ||
        library == NULL ||
        result == NULL ||
        target_peak == 0u ||
        max_slices > N5_SLICING_MAX_VARIABLES)
    {
        return 0;
    }

    memset(result, 0, sizeof(*result));
    n5_odd_charge_catalog_init(&charge_catalog);

    if (!n5_crt_get_primes(primes, N5_CRT_PRIME_COUNT))
    {
        return 0;
    }
    class_results = NULL;

    start_time = n5_parallel_wall_time_seconds();

    if (!n5_odd_charge_reduce(size_matrix, &charge_catalog) ||
        charge_catalog.case_count == 0u)
    {
        n5_odd_charge_catalog_free(&charge_catalog);
        return 0;
    }

    class_results = (N5OddProductionClassResult*)calloc(
        charge_catalog.case_count,
        sizeof(N5OddProductionClassResult)
    );

    if (class_results == NULL)
    {
        n5_odd_charge_catalog_free(&charge_catalog);
        return 0;
    }

    memset(&context, 0, sizeof(context));
    context.size_matrix = size_matrix;
    context.library = library;
    context.charge_catalog = &charge_catalog;
    context.target_peak = target_peak;
    context.max_slices = max_slices;
    context.class_results = class_results;

    result->worker_threads = n5_parallel_resolve_thread_count(
        requested_threads,
        charge_catalog.case_count
    );

    if (!n5_parallel_for(
            charge_catalog.case_count,
            result->worker_threads,
            n5_odd_production_adaptive_worker,
            &context
        ))
    {
        free(class_results);
        n5_odd_charge_catalog_free(&charge_catalog);
        return 0;
    }

    result->class_count = charge_catalog.case_count;
    result->raw_charge_count = charge_catalog.raw_charge_count;

    /*
     * Deterministic exact aggregation。
     *
     * Worker completion order 完全不影响结果：这里只按 canonical class index
     * 从 0 到 case_count-1 聚合 weight * residue。
     */
    for (class_index = 0u;
         class_index < charge_catalog.case_count;
         ++class_index)
    {
        const N5OddProductionClassResult* class_result;
        const N5OddAdaptiveClassResult* adaptive;
        uint64_t weight;
        unsigned int mode_index;

        class_result = &class_results[class_index];

        if (!class_result->success)
        {
            free(class_results);
            n5_odd_charge_catalog_free(&charge_catalog);
            return 0;
        }

        adaptive = &class_result->adaptive;
        weight = charge_catalog.cases[class_index].weight;
        mode_index = (unsigned int)adaptive->mode;

        if (mode_index >= 5u)
        {
            free(class_results);
            n5_odd_charge_catalog_free(&charge_catalog);
            return 0;
        }

        ++result->mode_count[mode_index];

        if (adaptive->original_primitive_nnz >
            result->max_original_primitive_nnz)
        {
            result->max_original_primitive_nnz =
                adaptive->original_primitive_nnz;
        }

        if (adaptive->selected_primitive_nnz >
            result->max_selected_primitive_nnz)
        {
            result->max_selected_primitive_nnz =
                adaptive->selected_primitive_nnz;
        }

        if (adaptive->original_peak > result->max_original_peak)
        {
            result->max_original_peak = adaptive->original_peak;
        }

        if (adaptive->selected_peak > result->max_selected_peak)
        {
            result->max_selected_peak = adaptive->selected_peak;
        }

        if (adaptive->sliced_peak > result->max_sliced_peak)
        {
            result->max_sliced_peak = adaptive->sliced_peak;
        }

        result->total_raw_slice_assignments +=
            adaptive->raw_slice_assignments;

        result->total_zero_slice_assignments +=
            adaptive->zero_slice_assignments;

        if (adaptive->max_dense_table_entries >
            result->max_dense_table_entries)
        {
            result->max_dense_table_entries =
                adaptive->max_dense_table_entries;
        }

        for (prime_index = 0u;
             prime_index < N5_CRT_PRIME_COUNT;
             ++prime_index)
        {
            uint32_t modulus;
            uint64_t weighted;

            modulus = primes[prime_index];

            weighted =
                (weight % modulus) *
                (uint64_t)adaptive->residues[prime_index];
            weighted %= modulus;

            result->residues[prime_index] = (uint32_t)(
                ((uint64_t)result->residues[prime_index] + weighted) % modulus
            );
        }
    }

    if (!n5_crt_reconstruct(
            result->residues,
            n5_crt_prime_basis(),
            N5_CRT_PRIME_COUNT,
            &result->crt
        ))
    {
        free(class_results);
        n5_odd_charge_catalog_free(&charge_catalog);
        return 0;
    }

    result->wall_seconds =
        n5_parallel_wall_time_seconds() - start_time;

    free(class_results);
    n5_odd_charge_catalog_free(&charge_catalog);
    return 1;
}
