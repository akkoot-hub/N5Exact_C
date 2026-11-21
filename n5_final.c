#include <inttypes.h>
#include <stdio.h>
#include <string.h>
#include <inttypes.h>

#include "n5_all_even.h"
#include "n5_catalog.h"
#include "n5_final.h"
#include "n5_odd_adaptive.h"
#include "n5_odd_kernel.h"
#include "n5_odd_production.h"
#include "n5_orbit.h"
#include "n5_parallel.h"

/*
 * Print a six-lane CRT residue vector for independent review.
 */
static void n5_final_print_residue_vector(
    const char* label,
    const uint32_t residues[N5_CRT_PRIME_COUNT]
)
{
    unsigned int index;

    printf("      %s residues = [", label);
    for (index = 0u; index < N5_CRT_PRIME_COUNT; ++index)
    {
        if (index != 0u)
        {
            printf(", ");
        }
        printf("%" PRIu32, residues[index]);
    }
    printf("]\n");
}

static void n5_final_print_matrix(const N5SizeMatrix* size_matrix)
{
    unsigned int row;
    unsigned int column;

    for (row = 0u; row < 4u; ++row)
    {
        if (row != 0u)
        {
            printf("; ");
        }

        for (column = 0u; column < 4u; ++column)
        {
            if (column != 0u)
            {
                printf(" ");
            }

            printf("%u", (unsigned int)size_matrix->matrix[row][column]);
        }
    }
}

/* Obtain the six-prime CRT basis used by the production calculation. */
static int n5_final_get_common_primes(
    uint32_t primes[N5_CRT_PRIME_COUNT]
)
{
    if (primes == NULL)
    {
        return 0;
    }
    return n5_crt_get_primes(primes, N5_CRT_PRIME_COUNT);
}

/*
 * 将 fixed-K contribution residues 乘 top-level orbit size。
 *
 * weighted_residues 返回 |O| C(K) 的 six-prime residues；
 * weighted_crt 返回对应的非负 exact integer reconstruction。
 */
static int n5_final_weight_orbit(
    const uint32_t residues[N5_CRT_PRIME_COUNT],
    unsigned int orbit_size,
    const uint32_t primes[N5_CRT_PRIME_COUNT],
    uint32_t weighted_residues[N5_CRT_PRIME_COUNT],
    N5CRTResult* weighted_crt
)
{
    unsigned int lane;

    if (residues == NULL || primes == NULL ||
        weighted_residues == NULL || weighted_crt == NULL)
    {
        return 0;
    }

    for (lane = 0u; lane < N5_CRT_PRIME_COUNT; ++lane)
    {
        weighted_residues[lane] = (uint32_t)(
            ((uint64_t)residues[lane] *
             (uint64_t)(orbit_size % primes[lane])) %
            primes[lane]
        );
    }

    return n5_crt_reconstruct(
        weighted_residues,
        primes,
        N5_CRT_PRIME_COUNT,
        weighted_crt
    );
}

/*
 * 把一组 weighted six-prime residues 累加到 sector 与 total accumulators。
 */
static void n5_final_accumulate_weighted(
    const uint32_t weighted[N5_CRT_PRIME_COUNT],
    const uint32_t primes[N5_CRT_PRIME_COUNT],
    uint32_t sector[N5_CRT_PRIME_COUNT],
    uint32_t total[N5_CRT_PRIME_COUNT]
)
{
    unsigned int lane;

    for (lane = 0u; lane < N5_CRT_PRIME_COUNT; ++lane)
    {
        sector[lane] = (uint32_t)(
            ((uint64_t)sector[lane] + weighted[lane]) % primes[lane]
        );

        total[lane] = (uint32_t)(
            ((uint64_t)total[lane] + weighted[lane]) % primes[lane]
        );
    }
}

/*
 * 将一个 CRT value 写成十进制；用于 per-orbit 输出。
 */
static int n5_final_crt_decimal(
    const N5CRTResult* crt,
    char output[N5_CRT_DECIMAL_BUFFER]
)
{
    if (crt == NULL || output == NULL)
    {
        return 0;
    }

    return n5_biguint_to_decimal(
        &crt->value,
        output,
        N5_CRT_DECIMAL_BUFFER
    );
}

/*
 * 输出从定义生成的顶层轨道目录。
 *
 * 这里打印的 ID 就是 canonical K row-major lexicographic order 中的 i+1，
 * 不是读取 frozen orbit numbering table。
 */
static void n5_final_print_catalog(const N5Catalog* catalog)
{
    unsigned int index;

    printf("\nDefinition-derived top-level orbit catalog:\n");
    printf(" ID  family    |O|  Stab  canonical K\n");

    for (index = 0u; index < catalog->orbit_count; ++index)
    {
        const N5TopOrbit* orbit;

        orbit = &catalog->orbits[index];

        printf("%3u  %-8s %4u %5u  ",
               index + 1u,
               n5_orbit_family_name(orbit->family),
               orbit->member_count,
               orbit->stabilizer_size);
        n5_final_print_matrix(&orbit->canonical);
        printf("\n");
    }
}

/*
 * 对一个 all-even orbit 调用已经认证的 generic all-even production engine，
 * 并把 weighted contribution 累加到 sector / total residues。
 */
static int n5_final_solve_all_even_orbit(
    unsigned int global_id,
    const N5TopOrbit* orbit,
    const N5FinalRunOptions* options,
    const uint32_t primes[N5_CRT_PRIME_COUNT],
    N5FinalRunResult* result
)
{
    N5AllEvenSolveResult solved;
    uint32_t weighted[N5_CRT_PRIME_COUNT];
    N5CRTResult weighted_crt;
    char weighted_decimal[N5_CRT_DECIMAL_BUFFER];

    memset(&solved, 0, sizeof(solved));
    memset(&weighted_crt, 0, sizeof(weighted_crt));

    if (!n5_all_even_solve(
            &orbit->canonical,
            options->requested_threads,
            &solved
        ) ||
        !n5_final_weight_orbit(
            solved.residues,
            orbit->member_count,
            primes,
            weighted,
            &weighted_crt
        ) ||
        !n5_final_crt_decimal(
            &weighted_crt,
            weighted_decimal
        ))
    {
        return 0;
    }

    n5_final_accumulate_weighted(
        weighted,
        primes,
        result->all_even_residues,
        result->total_residues
    );

    ++result->processed_all_even_orbits;
    result->total_all_even_reduced_flows += solved.reduced_flow_count;
    result->all_even_orbit_seconds += solved.wall_seconds;

    if (solved.max_dense_table_entries > result->global_max_dense_table_entries)
    {
        result->global_max_dense_table_entries = solved.max_dense_table_entries;
    }

    if (options->print_orbit_results)
    {
        printf("      method=all-even  reduced_flows=%u raw=%" PRIu64
               " threads=%u time=%.3fs\n",
               solved.reduced_flow_count,
               solved.weighted_flow_count,
               solved.worker_threads,
               solved.wall_seconds);
        printf("      lambda raw/reps=%" PRIu64 "/%" PRIu64
               " max_dense=%" PRIu64 "\n",
               solved.total_raw_lambda_slices,
               solved.total_slice_representatives,
               solved.max_dense_table_entries);
        printf("      C(K) = %s\n", solved.decimal_value);
        n5_final_print_residue_vector("C(K)", solved.residues);
        printf("      |O|C(K) = %s\n", weighted_decimal);
        n5_final_print_residue_vector("|O|C(K)", weighted);
        printf("\n");
    }

    (void)global_id;
    return 1;
}

/*
 * 对一个 odd orbit 调用已经认证的 adaptive odd production engine，
 * 并把 weighted contribution 累加到 sector / total residues。
 */
static int n5_final_solve_odd_orbit(
    unsigned int global_id,
    const N5TopOrbit* orbit,
    const N5OddKernelLibrary* library,
    const N5FinalRunOptions* options,
    const uint32_t primes[N5_CRT_PRIME_COUNT],
    N5FinalRunResult* result
)
{
    N5OddProductionResult solved;
    uint32_t weighted[N5_CRT_PRIME_COUNT];
    N5CRTResult weighted_crt;
    char c_decimal[N5_CRT_DECIMAL_BUFFER];
    char weighted_decimal[N5_CRT_DECIMAL_BUFFER];
    unsigned int mode;

    memset(&solved, 0, sizeof(solved));
    memset(&weighted_crt, 0, sizeof(weighted_crt));

    if (!n5_odd_production_solve_adaptive_modular(
            &orbit->canonical,
            library,
            options->requested_threads,
            N5_ODD_PRODUCTION_DEFAULT_TARGET_PEAK,
            N5_ODD_PRODUCTION_DEFAULT_MAX_SLICES,
            &solved
        ) ||
        !n5_final_weight_orbit(
            solved.residues,
            orbit->member_count,
            primes,
            weighted,
            &weighted_crt
        ) ||
        !n5_final_crt_decimal(&solved.crt, c_decimal) ||
        !n5_final_crt_decimal(&weighted_crt, weighted_decimal))
    {
        return 0;
    }

    n5_final_accumulate_weighted(
        weighted,
        primes,
        result->odd_residues,
        result->total_residues
    );

    ++result->processed_odd_orbits;
    result->total_odd_reduced_classes += solved.class_count;
    result->total_odd_raw_charge_count += solved.raw_charge_count;
    result->odd_orbit_seconds += solved.wall_seconds;
    result->odd_total_raw_slice_assignments += solved.total_raw_slice_assignments;
    result->odd_total_zero_slice_assignments += solved.total_zero_slice_assignments;

    for (mode = 0u; mode < 5u; ++mode)
    {
        result->odd_mode_count[mode] += solved.mode_count[mode];
    }

    if (solved.max_original_peak > result->odd_max_original_peak)
    {
        result->odd_max_original_peak = solved.max_original_peak;
    }

    if (solved.max_selected_peak > result->odd_max_selected_peak)
    {
        result->odd_max_selected_peak = solved.max_selected_peak;
    }

    if (solved.max_sliced_peak > result->odd_max_sliced_peak)
    {
        result->odd_max_sliced_peak = solved.max_sliced_peak;
    }

    if (solved.max_dense_table_entries > result->global_max_dense_table_entries)
    {
        result->global_max_dense_table_entries = solved.max_dense_table_entries;
    }

    if (options->print_orbit_results)
    {
        printf("      method=odd-adaptive classes=%u raw=%" PRIu64
               " threads=%u time=%.3fs\n",
               solved.class_count,
               solved.raw_charge_count,
               solved.worker_threads,
               solved.wall_seconds);
        printf("      modes[z/d/dr/s/sr]=%u/%u/%u/%u/%u\n",
               solved.mode_count[N5_ODD_ADAPTIVE_MODE_ZERO],
               solved.mode_count[N5_ODD_ADAPTIVE_MODE_DIRECT],
               solved.mode_count[N5_ODD_ADAPTIVE_MODE_DIRECT_REWRITTEN],
               solved.mode_count[N5_ODD_ADAPTIVE_MODE_SLICED],
               solved.mode_count[N5_ODD_ADAPTIVE_MODE_SLICED_REWRITTEN]);
        printf("      peak original/selected/sliced=%" PRIu64 "/%" PRIu64
               "/%" PRIu64 " raw_slices=%" PRIu64 "\n",
               solved.max_original_peak,
               solved.max_selected_peak,
               solved.max_sliced_peak,
               solved.total_raw_slice_assignments);
        printf("      C(K) = %s\n", c_decimal);
        n5_final_print_residue_vector("C(K)", solved.residues);
        printf("      |O|C(K) = %s\n", weighted_decimal);
        n5_final_print_residue_vector("|O|C(K)", weighted);
        printf("\n");
    }

    (void)global_id;
    return 1;
}

/*
 * 完成三个 residue sums 的 CRT reconstruction 与 decimal serialization。
 */
static int n5_final_reconstruct_totals(
    const uint32_t primes[N5_CRT_PRIME_COUNT],
    N5FinalRunResult* result
)
{
    return
        n5_crt_reconstruct(
            result->all_even_residues,
            primes,
            N5_CRT_PRIME_COUNT,
            &result->all_even_crt
        ) &&
        n5_crt_reconstruct(
            result->odd_residues,
            primes,
            N5_CRT_PRIME_COUNT,
            &result->odd_crt
        ) &&
        n5_crt_reconstruct(
            result->total_residues,
            primes,
            N5_CRT_PRIME_COUNT,
            &result->total_crt
        ) &&
        n5_biguint_to_decimal(
            &result->all_even_crt.value,
            result->all_even_decimal,
            N5_CRT_DECIMAL_BUFFER
        ) &&
        n5_biguint_to_decimal(
            &result->odd_crt.value,
            result->odd_decimal,
            N5_CRT_DECIMAL_BUFFER
        ) &&
        n5_biguint_to_decimal(
            &result->total_crt.value,
            result->total_decimal,
            N5_CRT_DECIMAL_BUFFER
        ) &&
        n5_biguint_to_decimal(
            &result->total_crt.modulus_product,
            result->modulus_decimal,
            N5_CRT_DECIMAL_BUFFER
        );
}

/*
 * 从五位输入定义开始执行完整 one-shot production。
 */
int n5_final_run(
    const N5FinalRunOptions* options,
    N5FinalRunResult* result
)
{
    N5Catalog catalog;
    N5OddKernelLibrary odd_library;
    uint32_t primes[N5_CRT_PRIME_COUNT];
    double start_time;
    double catalog_start;
    unsigned int orbit_index;

    if (options == NULL || result == NULL)
    {
        return 0;
    }

    memset(result, 0, sizeof(*result));
    memset(&catalog, 0, sizeof(catalog));
    memset(&odd_library, 0, sizeof(odd_library));

    start_time = n5_parallel_wall_time_seconds();
    catalog_start = start_time;

    if (!n5_catalog_build_from_definition(&catalog) ||
        !n5_odd_kernel_build_library(&odd_library) ||
        !n5_final_get_common_primes(primes))
    {
        n5_catalog_free(&catalog);
        return 0;
    }

    result->catalog_seconds =
        n5_parallel_wall_time_seconds() - catalog_start;

    printf("Definition-derived prerequisites:\n");
    printf("six-parameter candidates       = %" PRIu64 "\n", catalog.candidate_count);
    printf("legal K                        = %u\n", catalog.legal_count);
    printf("top symmetry group             = %u\n", catalog.top_group_count);
    printf("canonical top-level orbits     = %u\n", catalog.orbit_count);
    printf("all-even / odd                 = %u / %u\n",
           catalog.family_counts[N5_FAMILY_ALL_EVEN],
           catalog.orbit_count - catalog.family_counts[N5_FAMILY_ALL_EVEN]);
    printf("P6 / P8a / P8b / P10 / P16    = %u / %u / %u / %u / %u\n",
           catalog.family_counts[N5_FAMILY_P6],
           catalog.family_counts[N5_FAMILY_P8A],
           catalog.family_counts[N5_FAMILY_P8B],
           catalog.family_counts[N5_FAMILY_P10],
           catalog.family_counts[N5_FAMILY_P16]);
    printf("catalog + KernelLibrary time   = %.3f s\n", result->catalog_seconds);
    printf("CRT primes                     = %u / %u / %u / %u / %u / %u\n",
           primes[0], primes[1], primes[2], primes[3], primes[4], primes[5]);

    if (options->print_catalog)
    {
        n5_final_print_catalog(&catalog);
    }

    printf("\nSequential exact production:\n\n");

    for (orbit_index = 0u; orbit_index < catalog.orbit_count; ++orbit_index)
    {
        const N5TopOrbit* orbit;
        unsigned int global_id;
        int success;

        orbit = &catalog.orbits[orbit_index];
        global_id = orbit_index + 1u;

        if (options->print_orbit_results)
        {
            printf("[orbit %u] %-8s |O|=%u Stab=%u\n",
                   global_id,
                   n5_orbit_family_name(orbit->family),
                   orbit->member_count,
                   orbit->stabilizer_size);
            printf("      K = ");
            n5_final_print_matrix(&orbit->canonical);
            printf("\n");
            fflush(stdout);
        }

        if (orbit->family == N5_FAMILY_ALL_EVEN)
        {
            success = n5_final_solve_all_even_orbit(
                global_id,
                orbit,
                options,
                primes,
                result
            );
        }
        else
        {
            success = n5_final_solve_odd_orbit(
                global_id,
                orbit,
                &odd_library,
                options,
                primes,
                result
            );
        }

        if (!success)
        {
            fprintf(stderr, "Production failed at generated orbit %u.\n", global_id);
            n5_catalog_free(&catalog);
            return 0;
        }

        ++result->processed_orbits;
    }

    if (!n5_final_reconstruct_totals(primes, result))
    {
        n5_catalog_free(&catalog);
        return 0;
    }

    printf("\nAggregate six-prime residues:\n");
    n5_final_print_residue_vector("all-even subtotal", result->all_even_residues);
    n5_final_print_residue_vector("odd subtotal", result->odd_residues);
    n5_final_print_residue_vector("N5 total", result->total_residues);

    result->total_wall_seconds =
        n5_parallel_wall_time_seconds() - start_time;

    n5_catalog_free(&catalog);
    return 1;
}
