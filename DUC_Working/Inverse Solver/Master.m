%function [chi] = Master(configfile, N)
%Input is the config file
%Input N, number of discretizations of source
%Output is the value of chi in each cell
% all globals
clear all;
clc;
close all;

global theta;   %M x 1 matrix
global f;   %M x M matrix
global status; % for debug purposes
global side;  % side of square domain
global diameter; %diameter of measuring boundary
global aux_GD; % auxilliary GD N^2 x N^2
global aux_GS; % auxilliary GS M x N^2
global a;     % dimension of each of the N x N squares
global k_b; %k_b is a scalar
global aux_GD_star; 
global aux_GS_star;
global w_0;  % initialize the contrast source matrix N^2 x M
global chi_0; % initialize the contrast matrix N^2 x 1
global eta_D; % is a scalar that depends on the current iteration
global eta_S;
global u_inc; % is an N^2 x M matrix
global Dom_Coords; % is an N^2 x 2 matrix 
global Source_Coords; % is an M x2 matrix
global M; % number of sources/ receivers

global w;
global chi;
global g_0;



N = 64;
configfile = 'square.txt';
[status, theta, lambda, diameter, tr_x, tr_y, bl_x, bl_y, f] = read_synthetic_from_fwd(configfile);


M = length(theta);
if (tr_x == tr_y)
	side = 2*tr_x;
else
	error('expecting a square domain, length = %s, breadth = %s', 2*tr_x, 2*tr_y);
end;

%% 


%Create Coordinate Matrices
Source_Coords = zeros(M,2);

Dom_Coords = Dom_Coord(side, N);

for i = 1:M
    [x,y] = Source_Coord(i, diameter, M);
    Source_Coords(i,1) = x;
    Source_Coords(i,2) = y;
end ;


%Radius of equivalent area circle of square cell
a = side/(N*sqrt(pi));
k_b = (2*pi)/(lambda);

%Calculate aux N^2 x N^2 integral matrix
%Calculate aux N^2 x N^2 integral matrix
aux_GD = aux_GD_compute(N);  %% CHANGE IT BACK
aux_GD_star = (aux_GD)';   %% (aux_GD)' is the conjugate transpose of the matrix aux_GD

pause;
aux_GD = (k_b*k_b).*(aux_GD);
aux_GD_star = (conj(k_b*k_b)).*(aux_GD_star);  %% k_b is a real number

%Calculate aux M x N^2 integral matrix
%Calcualte aux N^2 x M integral matrix
aux_GS = aux_GS_compute(N);  % CHANGE IT BACK!
aux_GS_star = (aux_GS)';   %% (aux_GS)' is the conjugate transpose of the matrix aux_GS

aux_GS = (k_b*k_b).*(aux_GS);
aux_GS_star = (conj(k_b*k_b)).*(aux_GS_star); %% k_b is a real number



%Populate u_inc, either manually or from FEM Code [u_inc] = IncidentField(thetas)
[u_inc] = IncidentField(theta, N);

[v_0] = Initialize(N);

chi_temp = reshape(chi_0, [N,N]);
x_temp = reshape(Dom_Coords(:,1), [N,N]);
y_temp = reshape(Dom_Coords(:,2), [N,N]);

surf(x_temp, y_temp, real(chi_temp));
colormap(copper)

%set(hSurf, 'FaceColor', [1 0 0], 'FaceAlpha', 0.5);
figure;

      for i = 1:N*N
           if( real(chi_0(i))<0)
			chi_0(i)= 1i*imag(chi_0(i));
		end;
       
         if(imag(chi_0(i))<0)
           chi_0(i)= real(chi_0(i));    
         end
      end;
	   


pause;

surf(x_temp, y_temp, imag(chi_temp));
 pause;

g_0 = ones(N*N,M);
w = w_0;
chi = chi_0;

display('first iteration');

cost = cost_function()
[g_n, v_n] = UpdateContrastSource(g_0, v_0, N);
cost = cost_function()

UpdateContrast(N,0);

display('CSI loop: ');

while (cost > 0.001)
    [g_n, v_n] = UpdateContrastSource( g_n, v_n, N);
    UpdateContrast(N,0);

     %uncomment this to run CSI with positive chi (both real and imag)
     

      for i = 1:N*N
           if( real(chi(i))<0)
			chi(i)= 1i*imag(chi(i));
		end;
       
         if(imag(chi(i))<0)
           chi(i)= real(chi(i));    
         end
      end;

    chi_temp = reshape(chi, [N,N]);    
    figure(1),surf(x_temp, y_temp, real(chi_temp)); 
    figure(2),surf(x_temp, y_temp, imag(chi_temp));
    cost = cost_function()
    %Graphically see chi on the mesh
    pause;
    
    
end;
