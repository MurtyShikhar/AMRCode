%Calculate aux M x N^2 GS integral matrix
%%%% Equations --  %%%%%
function[Aux_GS] = aux_GS_compute(N)

global Dom_Coords; % is an N^2 x 2 matrix 
global Source_Coords; % is an M x2 matrix
global M; % number of sources/ receivers
global a;     % dimension of each of the N x N squares
global k_b; % k_b is a scalar

Aux_GS = zeros(M,N*N);

    % Using Richmond's Formulae
for i = 1:M
    for l = 1:N*N
        p_il = sqrt((Source_Coords(i,1) - Dom_Coords(l,1))^2 + (Source_Coords(i,2) - Dom_Coords(l,2))^2);
       
        Aux_GS(i,l) = (1i/4)*(2*pi*a/k_b)*(conj(besselj(1, k_b*a))*besselh(0,1,k_b*p_il)); 
    end;
end;

   % Using Numerical Integration
% for i = 1:M
%     for l = 1:N*N
%         Num=101; %No. of divisions each cell is divided(let it be odd)
%         leng=side/(N*Num);% Length of finest cell divided
%             for i1=1:Num
%                 for j1=1:Num
%                     intcoordx=((-1-Num)/2 +i1)*leng;
%                     intcoordy=((-1-Num)/2 +j1)*leng;
%                     p_il = sqrt((Source_Coords(i,1) - (Dom_Coords(l,1)+intcoordx))^2 + (Source_Coords(i,2) - (Dom_Coords(l,2)+intcoordy))^2);
%                     Aux_GD(i,l)=Aux_GD(i,l)+(leng*leng)*(1i/4)*besselh(0,1,k_b*p_il);
%                 end;
%             end;
%         end;
%     end;
% end;
