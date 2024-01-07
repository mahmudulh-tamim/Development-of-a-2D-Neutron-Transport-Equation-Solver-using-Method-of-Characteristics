%Leonard Optimized Quadrature

mu=[0.099812, 0.395534, 0.891439];
w=[0.017620, 0.188561, 0.793819];

%gauss legendre
N_p=16;
[mu, w] =lgwt(N_p,-1,1);
polar_discretization_number=size(mu,1);