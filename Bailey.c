/* result */
/* http://journals.cambridge.org/action/displayAbstract?fromPage=online&aid=1736752 */

/* proof */
/* http://www.jstor.org/stable/113644?seq=3&uid=16784632&uid=3738032&uid=2&uid=3&uid=16756968&uid=67&uid=62&sid=21106479401483#page_scan_tab_contents */

/* links */
/* http://www.mscand.dk/article/viewFile/10471/8492 */
/* http://mathworld.wolfram.com/LegendrePolynomial.html */


double coeff_A(int r){
    double start = (double) 2.*r - 1.;
    
    for(j=r-1; j>=1; j--)  start *= (2.*j - 1.);
    
    if(r<=0) start = 1.;
    
    return start/gsl_sf_fact(r);
}


int productLegendrePolys_coefficients(int p, int q){
    // expression legitimate for p >= q.  If fed in wrong order, switch to begin with. 
    int swap, r; 

    if(p<q){
        swap =    p;
        
        p    =    q;
        q    = swap;
    }

    printf("\n\nL_i coefficient of the product L_%d L_%d", p, q);

    for(r=q; r>=0; r--){
        Interim = coeff_A(p-r)*coeff_A(r)*coeff_A(q-r)/coeff_A(p+q-r)*(2.*p + 2.*q - 4.*r + 1.)/(2.*p + 2.*q  - 2.*r + 1.);
    
        printf("\n%d \t %lf", p+q-2*r, Interim);
    }  

    return 0;
}


double productLegendrePolys(int p, int q, double mu){
    // expression legitimate for p >= q.  If fed in wrong order, switch to begin with. 
    int swap, r; 

    if(p<q){
        swap =    p;
        
        p    =    q;
        q    = swap;
    }

    Interim = 0.0;

    for(r=0; r<=q; r++)  Interim += coeff_A(p-r)*coeff_A(r)*coeff_A(q-r)/coeff_A(p+q-r)*(2.*p + 2.*q - 4.*r + 1.)/(2.*p + 2.*q  - 2.*r + 1.)*gsl_sf_legendre_Pl(p+q-2*r, mu);

    return Interim;
}


int clippingCoefficient_test(double r, double mu){
  printf("\n\nClipping coefficient test");
  
  double A, B, C;
  
  double xi_0, xi_2, xi_4;
  
  xi_0 = xi_2 = xi_4 = Pk_powerlaw_truncated_xi(r);
  
  A    = xi_0*gsl_sf_legendre_Pl(0, mu) + xi_2*gsl_sf_legendre_Pl(2, mu) + xi_4*gsl_sf_legendre_Pl(4, mu);
  
  // 'control result'
  A    = A*A;
    
  // 'new' coefficients: derived from Bailey.c
  B    = (xi_0*xi_0 + xi_2*xi_2/5. + xi_4*xi_4/9.)*gsl_sf_legendre_Pl(0, mu); 
 
  // 'original' coefficients, mathematica.
  C    = (xi_0*xi_0 + xi_2*xi_2/4. + 11.*xi_4*xi_4/180.)*gsl_sf_legendre_Pl(0, mu); 
    
  printf("\n\nClipping coefficient test.");
    
  printf("\n%e \t %e \t %e", A, B, C);
    
  return 0;    
}


int Bailey_test(){
  printf("\n\nBailey test");
  
  productLegendrePolys_coefficients(2, 2);    

  // printf("\n%lf \t %lf", gsl_sf_legendre_Pl(4, 0.13)*gsl_sf_legendre_Pl(4, 0.13),  productLegendrePolys(4, 4, 0.13));

  // printf("\n\n");
    
  // convergence of quadrupole coefficient of xi_4 W2_l.   
  // for(kk=0; kk<20; kk=kk+2)  productLegendrePolys_coefficients(2, kk);
  
  // for(kk=0; kk<20; kk=kk+2)  printf("\n%.6f \t %.6f", gsl_sf_legendre_Pl(2, 0.13)*gsl_sf_legendre_Pl(kk, 0.13),  productLegendrePolys(2, kk, 0.13));
    
  // clippingCoefficient_test(25., 0.13);
    
  return 0;    
}
