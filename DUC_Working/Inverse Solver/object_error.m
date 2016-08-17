function[cost] = object_error(w,chi)

global f;
global aux_GS;
global aux_GD;
global u_inc;
global M;

cost = 0;
aux_matrix = aux_GD*w;
for j = 1:M
    cost = cost + norm(chi.*u_inc(:,j) - w(:,j) + chi.*(aux_matrix(:,j)) ...
                       )^2;
end


