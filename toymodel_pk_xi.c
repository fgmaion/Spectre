// Pk model, power law with high k truncation.
double Pk_powerlaw_truncated(double k){
    double n  = 0.96;
    double k0 =  0.09;
    double kc =   3.;

    return 2.*pow(pi, 2.)*pow(k0, -3.)*pow(k/k0, n)*exp(-1.*k/kc);
}
 
 
double Pk_powerlaw_truncated_xi(double r){
    double dn = 0.96;
    double k0 =  0.09;
    double kc =   3.;
    
    double y  = kc*r;
    
    return pow(kc/k0, dn + 3.)*sin((2.+dn)*atan(y))*tgamma(2.+dn)*pow(y*pow(1.+y*y, 1.+dn/2.), -1.);
} 


int inputTruncPowerLawPk(){
    pt2Pk   = &Pk_powerlaw_truncated;
    
    pt2Xi   = &Pk_powerlaw_truncated_xi;

    return 0;
}


double Pk_powerlaw(double k, double r0, double gamma){
    return 4.*pi*pow(k, -3.)*pow(k*r0, gamma)*tgamma(2.-gamma)*sin((pi/2.)*(2.-gamma)); 
} 


double Pk_powerlaw_xi(double r, double r0, double gamma){
    return pow(r0/r, gamma);
}

double toyRSD_OnePlusOneHalfMuSq(int order){
    switch(order){
        case 0:
            return (7./6.);
        
        case 2:
            return (1./3.); 
    }
}

