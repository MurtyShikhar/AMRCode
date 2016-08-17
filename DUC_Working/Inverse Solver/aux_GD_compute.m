%Calculate aux N^2 x N^2 GD integral matrix

function[Aux_GD] = aux_GD_compute(N)

global Dom_Coords; % is an N^2 x 2 matrix
global a;     % dimension of each of the N x N squares
global k_b; % k_b is a scalar
global side;

Aux_GD = zeros(N*N,N*N);

        % Using Richmond's formulae
for i = 1:N*N
    for l = 1:N*N
        p_il = sqrt((Dom_Coords(i,1) - Dom_Coords(l,1))^2 + (Dom_Coords(i,2) - Dom_Coords(l,2))^2);
        
        if(i == l)
            Aux_GD(i,l) = (1i/4)*(2/(k_b*k_b))*(pi*k_b*a*besselh(1, 1, k_b*a) + 2*1i);
        else
            Aux_GD(i,l) = (1i/4)*(2*pi*a/k_b)*(conj(besselj(1, k_b*a))*besselh(0,1,k_b*p_il));
        end
    end
end


       % Using Numerical Integration
% Num=101; %No. of divisions each cell is divided(let it be odd)
% leng=side/(N*Num);% Length of finest cell divided
% for i1=1:Num
%     intcoordx(i1)=((-1-Num)/2 +i1)*leng;
% end;
% for j1=1:Num
%     intcoordy(j1)=((-1-Num)/2 +j1)*leng;
% end;
% for i = 1:N*N
%     for l = 1:N*N
%         aa=leng/sqrt(pi);
%         if(i==l)
%             for i1=1:Num
%                 for j1=1:Num
%                     p_il = sqrt((Dom_Coords(i,1) - (Dom_Coords(l,1)+intcoordx(i1)))^2 + (Dom_Coords(i,2) - (Dom_Coords(l,2)+intcoordy(j1)))^2);
%                     if(intcoordx(i1)==0 && intcoordy(j1)==0)
%                         Aux_GD(i,l)=Aux_GD(i,l)+ (1i/4)*(2/(k_b*k_b))*(pi*k_b*aa*besselh(1, 1, k_b*aa) + 2*1i);
%                     else
%                         Aux_GD(i,l)=Aux_GD(i,l)+(leng*leng)*(1i/4)*besselh(0,1,k_b*p_il);
%                     end;
%                 end;
%             end;
%         else
%             for i1=1:Num
%                 for j1=1:Num
%                     p_il = sqrt((Dom_Coords(i,1) - (Dom_Coords(l,1)+intcoordx(i1)))^2 + (Dom_Coords(i,2) - (Dom_Coords(l,2)+intcoordy(j1)))^2);
%                     Aux_GD(i,l)=Aux_GD(i,l)+(leng*leng)*(1i/4)*besselh(0,1,k_b*p_il);
%                 end;
%             end;
%         end;
%     end;
% end;
