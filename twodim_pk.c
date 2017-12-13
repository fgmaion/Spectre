int assign2DPkMemory(){
  twodim_pk           = (double **) realloc(twodim_pk, n0*n1*n2*sizeof(double*));

  for(j=0; j<n0*n1*n2; j++)  twodim_pk[j]    = (double *)  malloc(3*sizeof(double));

  d2_binnedpk         = (double **) realloc(d2_binnedpk, 50*sizeof(double* ));

  for(j=0; j<50; j++) d2_binnedpk[j] = (double *) malloc(50*sizeof(double));

  for(k=0; k<50; k++) for(j=0; j<50; j++)  d2_binnedpk[k][j] = 0.0;

  return 0;
}

int Cartesian2Dpk(int modeCount){
  DualBinning(modeCount, twodim_pk, d2_binnedpk);

  printf("\n\nWriting 2D pk.");
  /*
  sprintf(filepath,"%s/W1_Spectro_V7_2/data_v1.7/pk_2d/pk_2d_W%d_%.1lf_%.1lf_Jf_%d.dat", root_dir, fieldFlag, lo_zlim, hi_zlim, (int) foldfactor);

  output = fopen(filepath, "w");

  // corresponds to k=0.5 for **current** binning.
  for(k=0; k<kbin_no; k++){
    for(j=0; j<kbin_no; j++){
      // fprintf(output, "%e \t", d2_binnedpk[k][j]);
    }

    fprintf(output, "\n");
  }

  fclose(output);
  */
  return 0;
}

int DualBinning(int NumberModes, double** DualParamArray, double** BinnedDualParamArray){
  int m;

  int bin_no = 50;

  double firstBinLimits[bin_no];
  double secndBinLimits[bin_no];

  int modesPerBin[bin_no][bin_no];

  int firstColLowerBinIndex     = 0;
  int firstColUpperBinIndex     = 0;

  int secndColLowerBinIndex     = 0;
  int secndColUpperBinIndex     = 0;


  qsort(DualParamArray, NumberModes, sizeof(DualParamArray[0]), FirstColumnCompare);

  // Order by first then second column.
  printf("\nDual param array sorted.");

  for(j=0; j<bin_no; j++) firstBinLimits[j] = 0.01 + j*1./bin_no;
  for(j=0; j<bin_no; j++) secndBinLimits[j] = 0.01 + j*1./bin_no;

  // for(j=0; j<20; j++)  printf("\n%.3lf \t %.3lf \t %.3lf", DualParamArray[j][0], DualParamArray[j][1], DualParamArray[j][2]);

  for(j=0; j<bin_no; j++){
    for(k=0; k<bin_no; k++){
      BinnedDualParamArray[j][k] = 0.0;
      // mean_firstCol[j][k]        = 0.0;
      // mean_secndCol[j][k]        = 0.0;
      modesPerBin[j][k]          =   0;
    }
  }

  for(j=0; j<NumberModes; j++){
    if(DualParamArray[j][0]   >= firstBinLimits[0]){
      firstColLowerBinIndex = j;

      break;
    }
  }

  for(j=0; j<bin_no; j++){
    for(i=firstColLowerBinIndex; i<NumberModes; i++){
      if(DualParamArray[i][0] > firstBinLimits[j+1]){
        firstColUpperBinIndex = i;
        break;
      }
    }

    secndColLowerBinIndex = firstColLowerBinIndex;

    qsort(&DualParamArray[firstColLowerBinIndex], firstColUpperBinIndex - firstColLowerBinIndex, sizeof(DualParamArray[0]), SecondColumnCompare);

    // meets the lower bin limit
    for(i=secndColLowerBinIndex; i<NumberModes; i++){
      if(DualParamArray[i][1]   >= secndBinLimits[0]){
        secndColLowerBinIndex = i;

        break;
      }
    }

    // meets the upper bin limit.
    for(k=0; k<bin_no; k++){
      for(m=secndColLowerBinIndex; m<firstColUpperBinIndex; m++){
        if(DualParamArray[m][1] > secndBinLimits[k+1]){
          secndColUpperBinIndex = m;
          break;
        }
      }

      for(m=secndColLowerBinIndex; m<secndColUpperBinIndex; m++){
        BinnedDualParamArray[j][k]    += DualParamArray[m][2];
        // mean_firstCol[j][k]           += DualParamArray[m][0];
        // mean_secndCol[j][k]           += DualParamArray[m][1];
        modesPerBin[j][k]             += 1;
      }

      if(modesPerBin[j][k] != 0)  BinnedDualParamArray[j][k] /= modesPerBin[j][k];
      // if(modesPerBin[j][k] != 0)  mean_firstCol[j][k]        /= modesPerBin[j][k];
      // if(modesPerBin[j][k] != 0)  mean_secndCol[j][k]        /= modesPerBin[j][k];

      secndColLowerBinIndex = secndColUpperBinIndex;
    }

    firstColLowerBinIndex = firstColUpperBinIndex;
  }

  return 0;
}

