// -- GNU C Compiler allocates memory for VLAs on the stack.[5] VLAs, like all objects in C, are limited to SIZE_MAX bytes.
// -- A struct must have a fixed size known at compile time. If you want an array with a variable length, you have to dynamically allocate memory.

struct regress{
  int    modes_perbin[KBIN_NO];
  double mean_modk[KBIN_NO], detA[KBIN_NO];
  double Sum_Pi[KBIN_NO], Sum_Li[KBIN_NO], Sum_Li2[KBIN_NO], Sum_PiLi[KBIN_NO];
  double Monopole[KBIN_NO], Quadrupole[KBIN_NO];
  double logk_limits[KBIN_NO];
  
  int*   kind;
  
  double* kLi;
  double* kM2;
};


int regress_mem(struct regress* inst){
  inst->kind = calloc(n0*n1*nx, sizeof(int));
  
  inst->kLi  = calloc(n0*n1*nx, sizeof(double));
  inst->kM2  = calloc(n0*n1*nx, sizeof(double));
  
  return 0;
}
