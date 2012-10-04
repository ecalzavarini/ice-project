kt = 0.0333333
NY = 100

set xl 'dimensionless time k_{T} t / L^{2}'
set yl 'melting front position x_{m}/L'

p 'melt.dat' u (kt*$1/NY**2.):($2) t 'k_T = 0.033' w lp

