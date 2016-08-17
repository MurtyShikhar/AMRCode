function[g_n_chi, d_n] = UpdateContrast_TV(N, IfWeighted, g_prev, d_prev)

global chi;
global w;
global side;

%Calculate gn_TV, also finds delta, reciprocal size, volume
[g_n_TV, grad_chi_x, grad_chi_y, delta, V] = Calc_g_n_TV(chi, IfWeighted);	%N*N x 1 matrix

%Calculates g_n_chi

g_n_chi = grad_TV(w, chi, g_n_TV);

%Calculates direction d_n

d_n = direction_TV(w, chi, g_n_chi, g_prev, d_prev);

%Calculate alpha

alpha = Calc_alpha_TV(w, chi, d_n, grad_chi_x, grad_chi_y, delta, V);

