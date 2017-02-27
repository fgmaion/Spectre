int scale_Cov(int N){
  for(j=0; j<order; j++)  gsl_matrix_set(sigma_norm, j, j, sqrt((double) N)*gsl_matrix_get(sigma_norm, j, j));

  return 0;
}
