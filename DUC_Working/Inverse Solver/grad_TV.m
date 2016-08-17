function[grad]  = grad_TV(w,chi, g_n_TV)
global eta_D;
global eta_S;
global N;
global M;
global u_inc;

aux_matrix = zeros(N*N,M);
grad = zeros(N*N,1);
for j = 1:M
   aux_matrix(:,j) = chi.*u_inc(:,j);
end
cost = cost_function();
for i = 1:N*N
    denom = 0;
    numer = 0;
    for j = 1:M
        numer = numer + (w(i,j) - aux_matrix(i,j))*(u_inc(i,j))
        denom = denom + norm(u(i,j))^2;
    end
    
    grad(i) = (numer + cost*(g_n_TV(i)))/denom;
    
end



