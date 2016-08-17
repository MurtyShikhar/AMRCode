function[f_n_TV] = calc_f_n_TV(chi_prev);
global N;
global chi;
global w;

delta = sqrt(data_error(w, chi)*reciprocal*reciprocal);

%calculate fn_TV
V = side*side

%reshape chi_D
chi_2D = reshape(chi, N, N);
chi_prev_2d = reshape(chi_prev,N,N);
%take derivatives
[grad_chi_prev_x, grad_chi_prev_y] = gradient(chi_2D_prev);
[grad_chi_x, grad_chi_y] = gradient(chi_2D);
dv = (side/N)*(side/N)

f_n_TV = 0;

for i = 1:N
    for j = 1:N
        n1 = abs(grad_chi_x(i,j))^2 + abs(grad_chi_y(i,j))^2;
        n2 = delta^2;
        d1 = abs(grad_chi_prev_x(i,j))^2 + abs(grad_chi_prev_y(i, ...
                                                          j))^2;
        f_n_TV = (n1 + n2)/(d1 + n2)*dv + f_n_TV;
    end
end
