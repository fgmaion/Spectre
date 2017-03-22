int assign2DPkMemory(){
  twodim_pk           = (double **) realloc(twodim_pk, n0*n1*n2*sizeof(double*));

  for(j=0; j<n0*n1*n2; j++)  twodim_pk[j]    = (double *)  malloc(3*sizeof(double));

  d2_binnedpk         = (double **) realloc(d2_binnedpk, 50*sizeof(double* ));

  for(j=0; j<50; j++) d2_binnedpk[j] = (double *) malloc(50*sizeof(double));

  for(k=0; k<50; k++) for(j=0; j<50; j++)  d2_binnedpk[k][j] = 0.0;

  return 0;
}
