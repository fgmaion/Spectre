int fiberCollision_cat(){
  // int fiberSelected = Vipers_Num;

  double* fiber_flag;
  
  double dx, dy;

  fiber_flag = malloc(Vipers_Num*sizeof(double));

  for(j=0; j<Vipers_Num; j++)  fiber_flag[j] = 1.;

  printf("\n\nGenerating fiber collision catalogue");

  for(i=0; i<Vipers_Num; i++){
    printf("\n %d", i);

    Index             = (int) gsl_rng_uniform_int(gsl_ran_r, Vipers_Num);

    if(fiber_flag[Index] == 1.){
      for(j=0; j<Vipers_Num; j++){
	dx            = fabs(xCoor[Index] - xCoor[j]);

	dy            = fabs(yCoor[Index] - yCoor[j]);
  
	fiber_flag[j] = fiberCollision(dx, dy, 0.0254, 0.6087);
      }
    }
  }

  double remaining = 0.0;

  for(j=0; j<Vipers_Num; j++)  remaining += fiber_flag[j];

  remaining /= Vipers_Num;

  printf("\n\n%e", remaining);

  double observedVolume = 0.;

  return 0;
}
