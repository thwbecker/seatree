#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

/* MMSC_INT_TYPE should be int32_t or int64_t depending on the size of an
   integer in the Fortran code that will be calling these functions. */
#ifndef MMSC_INT_TYPE
#define MMSC_INT_TYPE int32_t
#endif

/***** mallocf *****/

void *mallocf_(MMSC_INT_TYPE *nwords) {
    void *iptout;
    
    iptout = malloc(*nwords);
    if (iptout == NULL) {
        fprintf(stderr, "MALLOCF: Out of memory.\n");
        exit(1);
    }
    
#ifdef MMSC_DEBUG
    fprintf(stderr, "MALLOCF: Allocated %ld bytes at %p\n",
            (long int)*nwords, iptout);
#endif
    
    return iptout;
}

/***** reallocf *****/
 
void *reallocf_(void **iptr, MMSC_INT_TYPE *nwords) {
    void *iptout;
    
    iptout = realloc(*iptr, *nwords);
    if (iptout == NULL) {
        fprintf(stderr, "REALLOCF: Reallocation error - aborting.\n");
        exit(1);
    }

#ifdef MMSC_DEBUG
    fprintf(stderr, "REALLOCF: Reallocated from %p to %p with %ld bytes\n",
            *iptr, iptout, (long int)*nwords);
#endif
    
    return iptout;
}

/***** freef *****/

void *freef_(void **iptr) {
    free(*iptr);

#ifdef MMSC_DEBUG
    fprintf(stderr, "FREEF: Freed %p\n", *iptr);
#endif
    
    return 0;
}
