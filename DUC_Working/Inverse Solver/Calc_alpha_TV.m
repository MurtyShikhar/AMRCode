function[alpha] = Calc_alpha_TV(w, chi, d_n, grad_chi_x, grad_chi_y)

global u_inc;

cost = cost_function();
%This gives F_s + F_D,n

b = (sqrt(V)*(

