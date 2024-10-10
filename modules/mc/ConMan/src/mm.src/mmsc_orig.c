int MALLOCF(nwords)
    int *nwords;
{
 
    int iptout;
    int nbytes;
 
    nbytes = (*nwords);
    if(!(iptout= malloc(nbytes))) {
         printf("MALLOCF: Out of memory.\n");
         exit(1);
    }
 
    return iptout;
 
}
 
/* int REALLOCF(iptr,nwords)
    int *iptr;
    int *nwords;
{
 
    int iptout;
    int nbytes;
 
    nbytes = (*nwords);
    if(!(iptout= realloc((char *)*iptr,nbytes))) {
         printf("REALLOCF: Reallocation error - aborting\n");
         exit(1);
    }
 
    return iptout;
 
} */
 
int FREEF(iptr)
    int *iptr;
{
 
    free((char *)*iptr);
 
    return 0;
 
}
 
int mallocf_(nwords)
    int *nwords;
{
 
    int iptout;
    int nbytes;
 
    nbytes = (*nwords);
    if(!(iptout= malloc(nbytes))) {
         printf("MALLOCF: Out of memory.\n");
         exit(1);
    }
 
    return iptout;
 
}
 
int reallocf_(iptr,nwords)
    int *iptr;
    int *nwords;
{
 
    int iptout;
    int nbytes;
 
    nbytes = (*nwords);
    if(!(iptout= realloc((char *)*iptr,nbytes))) {
         printf("REALLOCF: Reallocation error - aborting\n");
         exit(1);
    }
 
    return iptout;
 
}
 
 
int freef_(iptr)
    int *iptr;
{
 
    free((char *)*iptr);
 
    return 0;
 
}
 
int mallocf(nwords)
    int *nwords;
{
 
    int iptout;
    int nbytes;
 
    nbytes = (*nwords);
    if(!(iptout= malloc(nbytes))) {
         printf("MALLOCF: Out of memory.\n");
         exit(1);
    }
 
    return iptout;
 
}
 
/* int reallocf(iptr,nwords)
    int *iptr;
    int *nwords;
{
 
    int iptout;
    int nbytes;
 
    nbytes = (*nwords);
    if(!(iptout= realloc((char *)*iptr,nbytes))) {
         printf("REALLOCF: Reallocation error - aborting\n");
         exit(1);
    }
 
    return iptout;
 
}*/
 
 
int freef(iptr)
    int *iptr;
{
 
    free((char *)*iptr);
 
    return 0;
 
}
 
 
 
 
