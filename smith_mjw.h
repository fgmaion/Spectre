int halofit();

// treatment of amp is pretty horrible. 
double amp;

double omega_m(double aa, double om_m0, double om_v0);
double omega_v(double aa, double om_m0, double om_v0);

double CarollPressTurner(double om_m, double om_v);

double BondEfstathiou(double rk, double gams, double p_index, double sig8);

double sigint(double r);

double BunnWhite(double rk);
double EisensteinHu(double yy, double om_m, double om_b);
