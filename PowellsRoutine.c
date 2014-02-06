#define ITMAX 200  // Max number of iterations


void powell(float p[], float **xi, int n, float ftol, int *iter, float *fret, float (*func)(float [])){
    /* 
      Minimization of a function func of n variables.  Input consists of an initial starting point p[1..n];
      an initial matrix xi[1..n][1..n], whose columns contain the initial set of directions (usually the n 
      unit vectors); and ftol, the fractional tolerance in the function value such that failure to decrease
      by more than this amount in one iteration signals completion.  On output, p is set to the best point 
      found, xi is the then-current direction set, fret is the returned function value at p, and iter is the
      number of iterations taken.  The routine linmin is used. 
    */  

    void linmin(float p[], float xi[], int n, float *fret, float (*func)(float []));
    
    int   i, ibig, j;
    
    float del, fp, fptt, t, *pt, *ptt, *xit;
    
    pt    = vector(1, n);
    ptt   = vector(1, n); 
    xit   = vector(1, n);
    
    *fret = (*func)(p);
    
    for(j=1; j<=n; j++) pt[j] = p[j]; 
    // Save the initial point. 
    
    for(*iter=1;; ++(*iter)){
        fp   = (*fret);
        ibig =   0;
        del  = 0.0;  // Will be the biggest function decrease, marking the direction to be removed to minimise
                     // the accumulation of linear dependence in the basis vector set. 

        for(i=1; i<=n; i++){
            for(j=1; j<=n; j++) xit[j] = xi[j][i];  // Copy the direction. 
        
            fptt = (*fret);
        
            linmin(p, xit, n, fret, func);          // Minimise along it. 
    
            if(fabs(fptt - (*fret)) > del){
                del  = fabs(fptt - (*fret));        // Record if its the largest decrease so far. 
                ibig = i;                                                
        
            }   
        }

        if(2.0*fabs(fp - (*fret)) <= ftol*(fabs(fp) + fabs(*fret))){
            free_vector(xit, 1, n);
            free_vector(ptt, 1, n);
            free_vector(pt,  1, n);   
            return;
        }

        if(*iter == ITMAX) nrerror("Powell exceeding maximum iterations.");
    
        for(j=1; j<=n; j++){
            ptt[j] = 2.0*p[j] - pt[j];      // Construct the extrapolated point
            xit[j] =     p[j] - pt[j];      //               average direction moved
            pt[j] =     p[j];              // Save the old starting point. 
        }
    
        fptt = (*func)(ptt);                // Function evaluated at extrapolated point. 
    
        if(fptt < fp){
            t = 2.0*(fp - 2.0*(*fret) + fptt)*SQR(fp - (*fret) - del) - del*SQR(fp - fptt);
        
            if(t<0.0){
                linmin(p, xit, n, fret, func);  // Move to the minimum of the new direction
            
                for(j=1; j<= n; j++){
                    xi[j][ibig] =  xi[j][n];    // Save the new diction 
                    xi[j][n]    = xit[j];
                }
            }
        }
    
    // Back for another iteration. 
    }
}
