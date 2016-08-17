function [u_inc] = IncidentField(theta, N)

global Dom_Coords; % is an N^2 x 2 matrix 
global M; % number of sources/ receivers
global k_b; % k_b is a scalar

for i = 1:M
    for j = 1:N*N
        u_inc(j,i) = exp(1i*k_b*(Dom_Coords(j,1)*cos(theta(i)) + Dom_Coords(j,2)*sin(theta(i))));
    end
end


        