def covariance_corrections(Ns, Nb, Np):
    # correction factor for psi.
    P = (1 - (Nb + 1)/(Ns -1))
    
    print "P:%.4le" % P
    
    A = 2/((Ns-Nb-1.)*(Ns-Nb-4.))
    B = (Ns-Nb-2.)/((Ns-Nb-1.)*(Ns-Nb-4.))

    m = 1. + B*(Nb-Np)/(1. + A + B*(Np + 1.))

    print "A:%.4le, B:%.4le, m:%.4le" % (A, B, m)


# Pezzotta et al.
# covariance_corrections(153., 20., 4.)

# Wilson et al.
covariance_corrections(153., 40., 3.)  # k_max = 0.2
covariance_corrections(153., 50., 3.)  # k_max = 0.4
covariance_corrections(153., 56., 3.)  # k_max = 0.6 
covariance_corrections(153., 60., 3.)  # k_max = 0.8 

