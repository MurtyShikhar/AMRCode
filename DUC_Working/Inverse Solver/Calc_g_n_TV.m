function [g_n_TV, grad_chi_x, grad_chi_y, delta, V] = Calc_g_n_TV(chi_prev, IfWeighted);

reciprocal = N/size;

delta = sqrt(data_error(w, chi_prev)*reciprocal*reciprocal);

%Calculate gn_TV
V = side*side;

chi_2D = reshape(chi_prev, N, N);
[grad_chi_x, grad_chi_y] = gradient(chi_2D);

g_n_TV = zeros(N,N);
temp = zeros(N,N);
for i = 1:N
	for j = 1:N
		numer = [grad_chi_x(i,j), grad_chi_y(i,j)];
		denom = delta^2 + abs(grad_chi_x(i,j))^2 + abs(grad_chi_y(i,j))^2;
		temp(i,j) = (1/V).*numer./denom;	
	end
end

g_n_TV = divergence(temp);

g_n_TV = reshape(g_n_TV, N^2, 1);