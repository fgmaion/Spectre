double bispecJ(double modk1, double modk2, double mu){
    return 1. - (0.5 - (3./14.)*pow(Om_m, -1./143.)) + 0.5*mu*(k1/k2 + k2/k1) + pow(mu, 2.)*(0.5 - (3./14.)*pow(Om_m, -1./143.);
}

double bispecJPcombo(double modk1, double modk2, double mu){
    return (bispecJ(modk1, modk2, mu)/b1 + b2*pow(b1, -2.))*(*pt2Pk)(modk1)(*pt2Pk)(modk2);
}

int theory(int a, int b, int c, double k){
    // a, b, c specify the integer ratio of the sides. 

    double k1, k2, k3;
    double mu_alpha, mu_beta, mu_gamma;
    
    k1 = a*k;
    k2 = b*k;
    k3 = c*k;
    
    mu_alpha = (pow(k2, 2.) + pow(k3, 2.) - pow(k1, 2.))/(2.*k2*k3);
    
    mu_beta  = (pow(k1, 2.) + pow(k2, 2.) - pow(k3, 2.))/(2.*k1*k2);
    
    mu_gamma = (pow(k3, 2.) + pow(k1, 2.) - pow(k2, 2.))/(2.*k1*k3);
    

    return bispecJPcombo(k1, k2, mu12) + bispecJPcombo(k2, k3, mu23) + bispecJPcombo(k3, k1, mu31);
}

int calc_bispectra(){

    theory(1, 1, 2);

    return 0;
}
