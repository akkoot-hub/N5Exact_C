#include <errno.h>
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "n5_final.h"
#include "n5_parallel.h"

static int parse_thread_value(
    const char* text,
    unsigned int* requested
)
{
    char* end;
    unsigned long value;

    if (text == NULL || requested == NULL)
    {
        return 0;
    }

    if (strcmp(text, "auto") == 0)
    {
        *requested = 0u;
        return 1;
    }

    errno = 0;
    end = NULL;
    value = strtoul(text, &end, 10);

    if (text[0] == '\0' || end == text || *end != '\0' ||
        errno == ERANGE || value == 0ul || value > (unsigned long)UINT_MAX)
    {
        return 0;
    }

    *requested = (unsigned int)value;
    return 1;
}

static int parse_threads(
    int argc,
    char** argv,
    unsigned int hardware_threads,
    unsigned int* requested
)
{
    int i;
    int seen_threads;

    *requested = 0u;
    seen_threads = 0;

    for (i = 1; i < argc; ++i)
    {
        const char* value_text;

        value_text = NULL;

        if (strcmp(argv[i], "--threads") == 0 || strcmp(argv[i], "-j") == 0)
        {
            if (i + 1 >= argc || seen_threads)
            {
                return 0;
            }
            value_text = argv[++i];
        }
        else if (strncmp(argv[i], "--threads=", 10) == 0)
        {
            if (seen_threads)
            {
                return 0;
            }
            value_text = argv[i] + 10;
        }
        else if (strcmp(argv[i], "--help") == 0 || strcmp(argv[i], "-h") == 0)
        {
            printf("Usage: N5Exact [--threads auto|N] [-j N]\n");
            printf("  auto: use all available logical processors (default)\n");
            printf("  N:    use N total participating threads, from 1 to %u\n", hardware_threads);
            return 2;
        }
        else
        {
            return 0;
        }

        if (!parse_thread_value(value_text, requested))
        {
            return 0;
        }
        seen_threads = 1;
    }

    if (*requested > hardware_threads)
    {
        return 0;
    }

    return 1;
}

int main(int argc, char** argv)
{
    N5FinalRunOptions options;
    N5FinalRunResult result;
    unsigned int hardware_threads;
    int parse_result;
    double start;
    double elapsed;

    hardware_threads = n5_parallel_hardware_threads();
    parse_result = parse_threads(
        argc,
        argv,
        hardware_threads,
        &options.requested_threads
    );
    if (parse_result == 2)
    {
        return 0;
    }
    if (parse_result == 0)
    {
        fprintf(stderr,
                "Invalid thread setting. Use --threads auto or a value from 1 to %u.\n",
                hardware_threads);
        return 2;
    }

    options.print_catalog = 1;
    options.print_orbit_results = 1;
    memset(&result, 0, sizeof(result));

    printf("N5 exact production runner\n");
    printf("Input dimension: 5-bit\n");
    printf("Hardware logical processors: %u\n", hardware_threads);
    if (options.requested_threads == 0u)
    {
        printf("Participating threads: AUTO (up to %u available logical processors)\n",
               hardware_threads);
    }
    else
    {
        printf("Participating threads requested: %u\n", options.requested_threads);
        printf("Note: the main thread participates; a stage may use fewer threads when it has fewer independent tasks.\n");
    }

    start = n5_parallel_wall_time_seconds();
    if (!n5_final_run(&options, &result))
    {
        fprintf(stderr, "Production computation did not complete.\n");
        return 1;
    }
    elapsed = n5_parallel_wall_time_seconds() - start;

    printf("\n========================================\n");
    printf("COMPUTED RESULTS\n");
    printf("========================================\n");
    printf("Generated orbits: %u\n", result.processed_orbits);
    printf("All-even orbits: %u\n", result.processed_all_even_orbits);
    printf("Odd orbits: %u\n", result.processed_odd_orbits);
    printf("All-even subtotal: %s\n", result.all_even_decimal);
    printf("Odd subtotal: %s\n", result.odd_decimal);
    printf("N5 total: %s\n", result.total_decimal);
    printf("CRT modulus product: %s\n", result.modulus_decimal);
    printf("Total wall time: %.3f s\n", elapsed);
    printf("Total wall time: %.3f min\n", elapsed / 60.0);
    return 0;
}
