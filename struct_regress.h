// -- GNU C Compiler allocates memory for VLAs on the stack.[5] VLAs, like all objects in C, are limited to SIZE_MAX bytes.
// -- A struct must have a fixed size known at compile time. If you want an array with a variable length, you have to dynamically allocate memory.

#ifndef REGRESS_GUARD
#define REGRESS_GUARD

struct Regress{
  int    modes_perbin[KBIN_NO];
  double mean_modk[KBIN_NO], detA[KBIN_NO];
  double Sum_Pi[KBIN_NO], Sum_Li[KBIN_NO], Sum_Li2[KBIN_NO], Sum_PiLi[KBIN_NO];
  double Monopole[KBIN_NO], Quadrupole[KBIN_NO];
  double logk_limits[KBIN_NO];
  
  int*   kind;
  
  double* kLi;
  double* kM2;

} unit, half;

typedef struct Regress regress;

#endif
