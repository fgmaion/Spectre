/*
int IntegralConstraintCorrection(){    
    ConvPkZeroPoint     = SphConvZeroPoint();
    printf("\nConvolved P(k) zero point calculated to be: %e", ConvPkZeroPoint);
    
    // sprintf(filepath, "%s/Data/ConvolvedPk/IntegralConstraint_Uncorrected/midK_Pk_%ssplint_padded.dat", root_dir, surveyType);
    // output = fopen(filepath, "w");    
    // fprintf(output, "%f \t %f \n", 0.0, ConvolvedPkZeroPoint);  
    // for(j=0; j<kBinNumb-1; j++) fprintf(output, "%f \t %e \n", midKBin[j], ConvolvedPk[j]);
    // fclose(output);
    
    
    // Window fn. evaluated on the same regular grid in k on which the convolved P(k) is estimated, normalised such that W^2 = 1.0 at k=0.0
    // for(j=0; j<kBinNumb-1; j++) ConvolvedPk[j] -= splintWindowfunc(midKBin[j])*ConvolvedPkZeroPoint;

    // sprintf(filepath, "%s/Data/ConvolvedPk/IntegralConstraint_Corrected/midK_Pk_%ssplint_padded.dat", root_dir, surveyType);
    // output = fopen(filepath, "w");    
    // for(j=0; j<kBinNumb-1; j++) fprintf(output, "%f \t %e \n", midKBin[j], ConvolvedPk[j]);
    // fclose(output);

    return 0;
}
*/

int AnisoICC(){    
	// Integral constraint correction in the 3D anisotropic case. 
	
    ConvPkZeroPoint     = ConvolveCell(n2/2 + xwfKernelsize, n1/2 + ywfKernelsize, n0/2 + zwfKernelsize);
    
    printf("\nConvolved P(vec k) zero point calculated to be: %e", ConvPkZeroPoint);

    printf("\nWindow fn. zero point calculated to be:  %f", ZeroPointNorm());

    // printICC();

	for(kkshift= zminKernelshift; kkshift< zmaxKernelshift + 1; kkshift++){
	    for(jjshift= yminKernelshift; jjshift< ymaxKernelshift + 1; jjshift++){
	        for(iishift= xminKernelshift; iishift< xmaxKernelshift + 1; iishift++){
			    filterIndex   = (kkshift + zmaxKernelshift)*ywfKernelsize*xwfKernelsize + (jjshift + ymaxKernelshift)*xwfKernelsize + (iishift + xmaxKernelshift);
			    
			    convIndex     = (n0/2 + kkshift)*(n1-2*ywfKernelsize)*(n2-2*xwfKernelsize) + (n1/2 + jjshift)*(n2-2*xwfKernelsize) + (n2/2 + iishift);

			    convolvedPk3d[convIndex] -= ConvPkZeroPoint*windowFunc3D[filterIndex]/ZeroPointNorm();
            }
        }
	}

    return 0;
}


int printICC(){
  printf("\n\nIntegral constraint correction. %d \t %d \t %d \n", xmaxKernelshift, ymaxKernelshift, zmaxKernelshift);

  for(iishift= -1; iishift< 1 + 1; iishift++){
      for(jjshift= -1; jjshift< 1 + 1; jjshift++){
          for(kkshift= -1; kkshift< 1 + 1; kkshift++){
          
	          filterIndex   = (kkshift + zmaxKernelshift)*ywfKernelsize*xwfKernelsize + (jjshift + ymaxKernelshift)*xwfKernelsize + (iishift + xmaxKernelshift);
	          
	          convIndex     = (n0/2 + kkshift)*(n1-2*ywfKernelsize)*(n2-2*xwfKernelsize) + (n1/2 + jjshift)*(n2-2*xwfKernelsize) + (n2/2 + iishift);
	          
	          printf("%e \t", (ConvPkZeroPoint*windowFunc3D[filterIndex]/ZeroPointNorm())/convolvedPk3d[convIndex]);
	      }
      
          printf("\n");
      }
      
      printf("\n\n");
  }
  
  return 0;
}
