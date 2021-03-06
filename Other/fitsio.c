int list_head(char filepath[]){
  fitsfile *fptr;
  char      card[FLEN_CARD];
  int       status = 0, nkeys, ii, hdupos;           // MUST initialize status.
  
  fits_open_file(&fptr, filepath, READONLY, &status);

  fits_get_hdu_num(fptr, &hdupos);                   // Get the current HDU position
  
  for(; !status; hdupos++){
    fits_get_hdrspace(fptr, &nkeys, NULL, &status);  // get number of keywords

    printf("\n\nHeader listing for HDU number: %d", hdupos);
    
    for(int ii = 1; ii <= nkeys; ii++){
      fits_read_record(fptr, ii, card, &status);     // read and print each keyword.

      if (fits_read_record(fptr, ii, card, &status))  break;
      
      printf("\n %s", card);
    }
    
    fits_movrel_hdu(fptr, 1, NULL, &status);         // try to move to next HDU
  }

  if(status == END_OF_FILE)  status = 0; // reset after normal error.
  
  fits_close_file(fptr, &status);

  if(status){          // print any error message.
    fits_report_error(stderr, status);
  }
  
  return(status);
}

int list_body(char filepath[]){
  fitsfile *fptr;         // FITS file pointer, defined in fitsio.h */
  char      keyname[FLEN_KEYWORD], colname[FLEN_VALUE], coltype[FLEN_VALUE];
  int       status = 0;   // CFITSIO status value MUST be initialized to zero! */
  int       hdupos, hdutype, bitpix, naxis, ncols, ii;
  long      naxes[10], nrows;

  fits_open_file(&fptr, filepath, READONLY, &status);
  
  fits_get_hdu_num(fptr, &hdupos);  // Get the current HDU position

  for(; !status; hdupos++){   // Main loop for each HDU
    fits_get_hdu_type(fptr, &hdutype, &status);  // Get the HDU type

    printf("\nHDU #%d  ", hdupos);

    if(hdutype == IMAGE_HDU){   // primary array or image HDU
      fits_get_img_param(fptr, 10, &bitpix, &naxis, naxes, &status);

      printf("Array:  NAXIS = %d,  BITPIX = %d\n", naxis, bitpix);

      for(ii = 0; ii < naxis; ii++){
        printf("   NAXIS%d = %ld\n",ii+1, naxes[ii]);
      }
    }

    else{  // a table HDU 
      fits_get_num_rows(fptr, &nrows, &status);
      fits_get_num_cols(fptr, &ncols, &status);

      if(hdutype == ASCII_TBL){
        printf("ASCII Table:  ");
      }
        
      else{
        printf("Binary Table:  ");
      }

      printf("%d columns x %ld rows\n", ncols, nrows);
      printf(" COL NAME             FORMAT\n");

      for(ii = 1; ii <= ncols; ii++){
        fits_make_keyn("TTYPE", ii, keyname, &status); /* make keyword */
        fits_read_key(fptr, TSTRING, keyname, colname, NULL, &status);
        fits_make_keyn("TFORM", ii, keyname, &status); /* make keyword */
        fits_read_key(fptr, TSTRING, keyname, coltype, NULL, &status);

        printf(" %3d %-16s %-16s\n", ii, colname, coltype);
      }
    }

    fits_movrel_hdu(fptr, 1, NULL, &status);  // try move to next ext.
  }

  if(status == END_OF_FILE) status = 0;

  fits_close_file(fptr, &status);

  if(status){
    fits_report_error(stderr, status); /* print any error message */
  }

  return(status);
}

int list_table(){
  fitsfile *fptr;      

  char     *val, value[1000], nullstr[]="*";
  char      keyword[FLEN_KEYWORD], colname[FLEN_VALUE];

  int       status = 0;   /*  CFITSIO status value MUST be initialized to zero!  */
  int       hdunum, hdutype, ncols, ii, anynul, dispwidth[1000];
  int       firstcol, lastcol = 0, linewidth;
  
  long      jj, nrows;
  long      frow, felem, nelem, longnull;
  
  float     floatnull = 0.;
  
  
  // see also: https://heasarc.gsfc.nasa.gov/docs/software/fitsio/cexamples/cookbook.c
  fits_open_file(&fptr, "/home/mjw/venice-4.0.2/nagoya_v6_samhain_W1W4.fits", READONLY, &status);

  if(fits_get_hdu_num(fptr, &hdunum) == 1){
    fits_movabs_hdu(fptr, 2, &hdutype, &status);  // This is the primary array;  try to move to the first extension and see if it is a table
  }
  
  else{
    fits_get_hdu_type(fptr, &hdutype, &status); // Get the HDU type
  }
  
  if(hdutype == IMAGE_HDU)  printf("Error: this program only displays tables, not images\n");

  else{
    fits_get_num_rows(fptr, &nrows, &status);
    fits_get_num_cols(fptr, &ncols, &status);
   
    while(lastcol < ncols){
      linewidth = 0;                                            // find the number of columns that will fit within 80 characters
      firstcol = lastcol + 1;
      
      for(lastcol = firstcol; lastcol <= ncols; lastcol++){
        fits_get_col_display_width(fptr, lastcol, &dispwidth[lastcol], &status);

        linewidth += dispwidth[lastcol] + 1;

        if(linewidth > 80)  break;
      }

      if(lastcol > firstcol)  lastcol--;  // the last col didn't fit

      // print column names as column headers
      printf("\n    ");

      for(ii=firstcol; ii<=lastcol; ii++){
        fits_make_keyn("TTYPE", ii, keyword, &status);

        fits_read_key(fptr, TSTRING, keyword, colname, NULL, &status);

        colname[dispwidth[ii]] = '\0';  /* truncate long names */

        printf("%*s ",dispwidth[ii], colname);
      }

      printf("\n");  /* terminate header line */

      // print each column, row by row (there are faster ways to do this).
      val = value;
      
      for(jj = 1; jj <= 20 && !status; jj++){
          printf("%4d ", jj);
          
          for (ii = firstcol; ii <= lastcol; ii++){
            /* read value as a string, regardless of intrinsic datatype */
            if(fits_read_col_str(fptr, ii, jj, 1, 1, nullstr, &val, &anynul, &status)){
              break;  /* jump out of loop on error */
            }
            
            printf("%-*s ", dispwidth[ii], value);
          }

          printf("\n");
      }
    }
  }
  
  fits_close_file(fptr, &status);

  if(status) fits_report_error(stderr, status); /* print any error message */

  return(status);
}

int read_maskfits(char filepath[], double sampling){
  // e.g. filepath = "/home/mjw/venice-4.0.2/nagoya_v6_samhain_W1W4.fits"
  fitsfile *fptr;
  
  char      nullstr[]="*";

  int       status = 0;   //  CFITSIO status value MUST be initialized to zero!
  int       hdunum, hdutype, ncols, anynul;

  long      nrows;
  
  double    *ra,  *dec;
  double    doublenull = 0.;

  list_head(filepath);  
  
  // see also: https://heasarc.gsfc.nasa.gov/docs/software/fitsio/cexamples/cookbook.c
  fits_open_file(&fptr, filepath, READONLY, &status);

  if(fits_get_hdu_num(fptr, &hdunum) == 1){
    fits_movabs_hdu(fptr, 2, &hdutype, &status);  // This is the primary array;  try to move to the first extension and see if it is a table
  }

  else{
    fits_get_hdu_type(fptr, &hdutype, &status); // Get the HDU type
  }

  if(hdutype == IMAGE_HDU)  printf("Error: this program only displays tables, not images\n");

  else{
    fits_get_num_rows(fptr, &nrows, &status);
    fits_get_num_cols(fptr, &ncols, &status);

    rand_number = (int) floor(nrows*sampling);
        
    rand_ra     = calloc(rand_number, sizeof(*rand_ra)); // sizeof(double).
    rand_dec    = calloc(rand_number, sizeof(*rand_dec));
    
    fits_read_col(fptr, TDOUBLE, 1, 1, 1, rand_number, &doublenull, rand_ra,  &anynul, &status);
    fits_read_col(fptr, TDOUBLE, 2, 1, 1, rand_number, &doublenull, rand_dec, &anynul, &status);
  }

  fits_close_file(fptr, &status);

  if(status) fits_report_error(stderr, status); /* print any error message */
    
  // for(j=0; j<20; j++)  printf("\n%.13lf \t %.13lf", ra[j], dec[j]);
    
  return(status);
}
