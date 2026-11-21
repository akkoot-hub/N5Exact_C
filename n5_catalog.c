#include <stdlib.h>
#include <string.h>

#include "n5_catalog.h"

static int n5_catalog_matrix_from_six_parameters(
    int p, int q, int r, int s, int t, int h,
    N5SizeMatrix* output
)
{
    int values[4][4];
    unsigned int row;
    unsigned int column;

    if (output == NULL)
    {
        return 0;
    }

    values[0][0] = p + r + s + t + h - 8;
    values[0][1] = h + p - q;
    values[0][2] = 8 + q - s - t - h - p;
    values[0][3] = 8 - r - h - p;
    values[1][0] = 8 + q - r - t - h - p;
    values[1][1] = 8 - s - h - p;
    values[1][2] = r + s + t + 2 * h + p - q - 8;
    values[1][3] = p;
    values[2][0] = t + h - q;
    values[2][1] = q;
    values[2][2] = 8 - r - t - h;
    values[2][3] = r;
    values[3][0] = 8 - s - t - h;
    values[3][1] = s;
    values[3][2] = t;
    values[3][3] = h;

    for (row = 0; row < 4; ++row)
    {
        for (column = 0; column < 4; ++column)
        {
            if (values[row][column] < 0 || values[row][column] > 8)
            {
                return 0;
            }
            output->matrix[row][column] = (uint8_t)values[row][column];
        }
    }

    return n5_size_matrix_is_legal(output);
}

static unsigned int n5_catalog_enumerate_legal_matrices(
    N5SizeMatrix** output,
    uint64_t* candidate_count
)
{
    N5SizeMatrix* matrices;
    unsigned int capacity;
    unsigned int count;
    uint64_t candidates;
    int p, q, r, s, t, h;

    if (output == NULL || candidate_count == NULL)
    {
        return 0;
    }

    *output = NULL;
    *candidate_count = 0;
    capacity = 256u;
    count = 0u;
    candidates = 0u;

    matrices = (N5SizeMatrix*)malloc((size_t)capacity * sizeof(*matrices));
    if (matrices == NULL)
    {
        return 0;
    }

    for (p = 0; p <= 8; ++p)
    for (q = 0; q <= 8; ++q)
    for (r = 0; r <= 8; ++r)
    for (s = 0; s <= 8; ++s)
    for (t = 0; t <= 8; ++t)
    for (h = 0; h <= 8; ++h)
    {
        N5SizeMatrix matrix;
        ++candidates;

        if (!n5_catalog_matrix_from_six_parameters(p, q, r, s, t, h, &matrix))
        {
            continue;
        }

        if (count == capacity)
        {
            N5SizeMatrix* resized;
            unsigned int new_capacity = capacity << 1u;
            resized = (N5SizeMatrix*)realloc(
                matrices, (size_t)new_capacity * sizeof(*matrices));
            if (resized == NULL)
            {
                free(matrices);
                return 0;
            }
            matrices = resized;
            capacity = new_capacity;
        }

        matrices[count++] = matrix;
    }

    *candidate_count = candidates;
    *output = matrices;
    return count;
}

void n5_catalog_free(N5Catalog* catalog)
{
    if (catalog != NULL)
    {
        free(catalog->orbits);
        memset(catalog, 0, sizeof(*catalog));
    }
}

int n5_catalog_build_from_definition(N5Catalog* catalog)
{
    N5SizeMatrix* legal_matrices = NULL;
    N5CellPermutation* group = NULL;
    N5TopOrbit* orbits = NULL;
    uint64_t candidate_count = 0u;
    unsigned int legal_count;
    unsigned int group_count;
    unsigned int orbit_count;
    unsigned int orbit_index;

    if (catalog == NULL)
    {
        return 0;
    }

    memset(catalog, 0, sizeof(*catalog));

    legal_count = n5_catalog_enumerate_legal_matrices(
        &legal_matrices, &candidate_count);
    if (legal_count == 0u || legal_matrices == NULL)
    {
        goto fail;
    }

    group = (N5CellPermutation*)malloc(
        (size_t)N5_TOP_GROUP_SIZE * sizeof(*group));
    orbits = (N5TopOrbit*)calloc((size_t)legal_count, sizeof(*orbits));
    if (group == NULL || orbits == NULL)
    {
        goto fail;
    }

    group_count = n5_generate_top_group(group);
    if (group_count == 0u)
    {
        goto fail;
    }

    orbit_count = n5_build_top_orbits(
        legal_matrices, legal_count, group, group_count, orbits, legal_count);
    if (orbit_count == 0u)
    {
        goto fail;
    }

    catalog->candidate_count = candidate_count;
    catalog->legal_count = legal_count;
    catalog->top_group_count = group_count;
    catalog->orbit_count = orbit_count;
    catalog->orbits = orbits;
    orbits = NULL;

    for (orbit_index = 0u; orbit_index < catalog->orbit_count; ++orbit_index)
    {
        unsigned int family = (unsigned int)catalog->orbits[orbit_index].family;
        if (family >= 7u)
        {
            goto fail;
        }
        ++catalog->family_counts[family];
    }

    free(group);
    free(legal_matrices);
    return 1;

fail:
    free(orbits);
    free(group);
    free(legal_matrices);
    n5_catalog_free(catalog);
    return 0;
}
