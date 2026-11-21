#ifndef N5_CATALOG_H
#define N5_CATALOG_H

#include <stdint.h>

#include "n5_core.h"
#include "n5_orbit.h"

typedef struct
{
    uint64_t candidate_count;
    unsigned int legal_count;
    unsigned int top_group_count;
    unsigned int orbit_count;
    unsigned int family_counts[7];
    N5TopOrbit* orbits;
} N5Catalog;

/* Rebuild the legal size-matrix set and its symmetry-orbit catalog. */
int n5_catalog_build_from_definition(N5Catalog* catalog);

/* Release storage owned by a catalog. */
void n5_catalog_free(N5Catalog* catalog);

#endif
