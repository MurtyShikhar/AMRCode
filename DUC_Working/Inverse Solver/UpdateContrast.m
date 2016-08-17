function [] = UpdateContrast(N, flag)
% u is an N^2 x M matrix ,for M sources

% all globals
global chi; 
global w;
global aux_GD; % auxilliary GD N^2 x N^2
global eta_D; % is a scalar that depends on the current iteration
global u_inc; % is an N^2 x M matrix

global M; % number of sources/ receivers

% u_inc is an N^2 x M matrix

if(flag ~= 1)
	aux_matrix = (aux_GD)*w;

	%total_field is an N^2 x M matrix, one column per source,and one row for each of the N^2 squares
	%Update u (Equation 36)
	total_field = aux_matrix + u_inc;

	%Update contrast (Equation 38)

	for i = 1:N*N
		numer = 0;
		denom = 0;
		for j=1:M
			numer = numer + (w(i,j))*(conj(total_field(i,j)));
			denom = denom + norm(total_field(i,j))^2;
		end
		chi(i) = numer/denom;
	end
elseif (flag == 1)
end;
%Update eta_D
aux_eta_D=0;
for j=1:M
    aux_eta_D=aux_eta_D+(norm(chi.*u_inc(:,j)))^2;
end

eta_D=1/aux_eta_D