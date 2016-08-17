function[cost] = data_error(w)
global f;
global M;
global aux_GS


cost = 0;
aux_matrix =aux_GS*w;
for j=1:M
    cost = cost + norm(f(:,j) - aux_matrix(:,j) )^2;
end

