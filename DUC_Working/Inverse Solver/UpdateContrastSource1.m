function [g_n, v_n] = UpdateContrastSource(g_prev, v_prev , N)

% all globals
global f;   %M x M matrix
global chi; % N^2 x 1 matrix
global w;
global aux_GD; % auxilliary GD N^2 x N^2
global aux_GS; % auxilliary GS M x N^2
global aux_GD_star; 
global aux_GS_star;
global eta_S; % is a scalar
global eta_D; % is a scalar that depends on the current iteration
global u_inc; % is an N^2 x M matrix
global M; % number of sources/ receivers

%w is N^2 x M matrix, where each column is a w_j for the jth incident source
%chi is N^2 x 1 matrix
%eta_S is scalar
%eta_D is scalar
%f is MxM matrix. f_ij is reading at ith receiver and incident field from jth transmitter
%u_inc is N^2 x M matrix
%M is number of receivers, N^2 is number of cells
%side of object domain
%diameter of circle where receivers are placed
%k_b 
%g_prev is N^2 x M matrix 
%v_prev is N^2 x M matrix 

%Dom_Coords is N^2 x 2 matrix of (x,y) pairs
%Source_Coords is M x 2 matrix of (x,y) pairs

%BEGIN------------------------------------------------------------------------

%Calculate DATA ERROR

%dat_error(i,j) is the data error for the jth source, and the ith receiver.
dat_error = zeros(M,M);
%aux_matrix is the matrix product of aux_GS and w. it is an MxM matrix
aux_matrix = (aux_GS)*w;
for i = 1:M
    for j = 1:M
        dat_error(i,j) = f(i,j) - aux_matrix(i,j);
    end
end
   
% object error is an N^2 xM matrix, one column for each source
obj_error = zeros(N*N, M);
aux_matrix = (aux_GD)*w;

for i = 1:N*N
    for j = 1:M
        obj_error(i,j) = chi(i)*u_inc(i,j) - w(i,j) +  chi(i)*aux_matrix(i,j);
    end
end

%Calculate G_S*

%change
G_s_star = conj(aux_GS_star)* conj(dat_error);

%Calculate G_D*
chi_obj_mult = zeros(N*N, M);
for i=1:N*N
    %change
    chi_obj_mult(i,:) = chi(i).*(conj(obj_error(i,:)));
end

chi_obj_mult_2 = zeros(N*N, M);
for i=1:N*N
    %change
    chi_obj_mult_2(i,:) = conj(chi(i)).*obj_error(i,:);
end

%change
G_d_star = conj(aux_GD_star)*chi_obj_mult;

%Calculate g_n which is N^2 x M mtrix
g_n = (-1*eta_S) .* G_s_star - eta_D .* (obj_error - G_d_star);

alpha_g_n = (-1*eta_S) .* conj(aux_GS_star) * dat_error - eta_D .* (obj_error - conj(aux_GD_star)*chi_obj_mult_2);

%Calculate v_n which is a N^2 x M matrix

sum_numer = 0;
sum_denom = 0;

for i = 1:M
    sum_numer = sum_numer + real(dot(g_n(:,i), g_n(:,i) - g_prev(:,i)));
    sum_denom = sum_denom + dot(g_prev(:,i), g_prev(:,i));
end
%%
% *BOLD TEXT*

%disp(sum_numer)
%disp(sum_denom)
%disp(v_prev)

v_n = g_n + (sum_numer/sum_denom).*v_prev;

%Calculate alpha

%line search for calculating alpha

% alpha_n = linesearch(g_n, v_n, w, chi);
% alpha_n = alpha_n*(-1)
% pause;

sum_numer = 0;
alpha_denom1 = 0;
alpha_denom2 = 0;


mat_1 = aux_GS*v_n;
mat_2 = aux_GD*v_n;

for i = 1:N*N
    mat_2(i,:) = chi(i)*mat_2(i,:);
end

for i = 1:M
    sum_numer = sum_numer - real(dot(alpha_g_n(:,i), v_n(:,i))); %change
    alpha_denom1 = alpha_denom1 + norm(mat_1(:,i))^2;
    alpha_denom2 =alpha_denom2 + norm(v_n(:,i) - mat_2(:,i))^2;
end

alpha_denom1 = alpha_denom1*(eta_S);
alpha_denom2 = alpha_denom2*(eta_D);
sum_denom = alpha_denom1 + alpha_denom2;

alpha_n = sum_numer/sum_denom

%Update contrast source
w = w  + alpha_n.*v_n;   %%% move in the cojugate direction
%w = w - alpha_n.*v_n;  %%%%%% move in the opposite direction
%DONE-------------------------------

%Calculate data error (Equation 26)
    %Calculate G_D
    %Calculate G_S

%Calculate object error (Equation 27)

%Update w (Equation 28)
    %G_D*,  (Equation 32)
    %G_S*   (Equation 31)
    %Calculate g, (Equation 30)
        %Calculate alpha (Equation 34)
        %Calculate v (Equation 29)
