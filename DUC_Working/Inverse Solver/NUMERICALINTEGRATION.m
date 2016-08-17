clc;close all;
clear all;
N=2;
k_b=2*pi;
side=1/5;
Dom_Coords=Dom_Coord(side,N);
a = side/(N*sqrt(pi));
M=24;



Aux_GD = zeros(N*N,N*N);
Aux_GD1 = zeros(N*N,N*N);
        % Using Richmond's formulae
for i = 1:N*N
    for l = 1:N*N
        p_il = sqrt((Dom_Coords(i,1) - Dom_Coords(l,1))^2 + (Dom_Coords(i,2) - Dom_Coords(l,2))^2);
        
        if(i == l)
            Aux_GD1(i,l) = (1i/4)*(2/(k_b*k_b))*(pi*k_b*a*besselh(1, 1, k_b*a) + 2*1i);
        else
            Aux_GD1(i,l) = (1i/4)*(2*pi*a/k_b)*(conj(besselj(1, k_b*a))*besselh(0,1,k_b*p_il));
        end;
    end;
end;


       % Using Numerical Integration
for i = 1:N*N
    for l = 1:N*N
        Num=1001; %No. of divisions each cell is divided(let it be odd)
        leng=side/(N*Num);% Length of finest cell divided
        aa=leng/sqrt(pi);
        if(i==l)
            for i1=1:Num
                for j1=1:Num
                    intcoordx=((-1-Num)/2 +i1)*leng;
                    intcoordy=((-1-Num)/2 +j1)*leng;
                    p_il = sqrt((Dom_Coords(i,1) - (Dom_Coords(l,1)+intcoordx))^2 + (Dom_Coords(i,2) - (Dom_Coords(l,2)+intcoordy))^2);
                    if(intcoordx==0 && intcoordy==0)
                        Aux_GD(i,l)=Aux_GD(i,l)+ (1i/4)*(2/(k_b*k_b))*(pi*k_b*aa*besselh(1, 1, k_b*aa) + 2*1i);
                    else
                        Aux_GD(i,l)=Aux_GD(i,l)+(leng*leng)*(1i/4)*besselh(0,1,k_b*p_il);
                    end;
                end;
            end;
        else
            for i1=1:Num
                for j1=1:Num
                    intcoordx=((-1-Num)/2 +i1)*leng;
                    intcoordy=((-1-Num)/2 +j1)*leng;
                    p_il = sqrt((Dom_Coords(i,1) - (Dom_Coords(l,1)+intcoordx))^2 + (Dom_Coords(i,2) - (Dom_Coords(l,2)+intcoordy))^2);
                    Aux_GD(i,l)=Aux_GD(i,l)+(leng*leng)*(1i/4)*besselh(0,1,k_b*p_il);
                end;
            end;
        end;
    end;
end;
error=100.*abs(Aux_GD-Aux_GD1)./abs(Aux_GD)