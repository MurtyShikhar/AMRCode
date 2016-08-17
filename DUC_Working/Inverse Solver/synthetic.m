%Synthetic data generation using MoM. Re: paper by Richmond (1965) IEEE TAP.
clc;
%Problem geometry description
a = 0.5; %side of the square

%synthetic problem specs
dl = 1.0/10; %discretization of cell
dr = dl/sqrt(pi); %radius of equivalent circle with same area
n = (ceil(a/dl))^2; %number of pixels
xc = (-a/2+dl/2):dl:(a/2-dl/2); %x-centre of pixels
yc = (-a/2+dl/2):dl:(a/2-dl/2); %y-centre of pixels
[XC,YC] = meshgrid(xc,yc);

eps_true = 3; %homogeneous true epsilon
eps_mat = ones(sqrt(n),sqrt(n)) * eps_true;
x = eps_mat(:) - eb; %contrast variable
eb = 1; %epsilon of background medium
k = 2*pi*sqrt(eb); %wavevector of background medium
angles = 10; %number of incident waves
thetas = 0:2*pi/angles:(2*pi*(1-1/angles)); %list of incident angles

%measurement points
m = 15; %no of measurements for each inc angle
m_thetas = 0:2*pi/m:(2*pi*(1-1/m)); %list of measured angles
r_meas = a + 0.5;
xm = r_meas * cos(m_thetas);
ym = r_meas * sin(m_thetas);
s = zeros(m*angles,1); %vector holding the scattered field
    
%Variable for synthetic or inverse problem
gen_synthetic = true;

%The matrix eqn we are solving is as follows. Step 1, solve Fu=b:
%F(p,q) = 1+(j/2)[(pi k dr) H_1^(2)(k dr) - 2j]x(q), when p=q
% = (j pi k dr/2) J_1(k dr)H_0^(2)(k rho_pq)x(q), when p\neq q
%where p,q both belong to D and x = \epsilon_r - 1.
%Now that u is solved for, step 2, obtain s = M u:
%M(p,q) = -(j pi k dr/2) J_1(k dr)H_0^(2)(k rho_pq)x(q), where p belongs
%to S and q belongs to D, where D represents domain (state eqn) and S is
%the measurement domain (data eqn). Repeat for each incidence angle.

%some often recurring terms
diag_el = (1i/2)*(pi*k*dr*besselh(1,2,k*dr)-2*1i);
off_diag_el = (1i*pi*k*dr/2)*besselj(1,k*dr);

F = zeros(n); %initialize with zeros
b = zeros(n,1); %initialize rhs with zeros
%build fwd matrix
for p=1:n
    for q=1:n
        if(p==q)
            F(p,q) = 1+ (eps_mat(p)-1)*diag_el;
        else
            rho = sqrt((XC(p)-XC(q))^2+(YC(p)-YC(q))^2);
            F(p,q) = off_diag_el*(eps_mat(q)-1)*besselh(0,2,k*rho);
        end
    end
end
%solve and compute for scattered fields
for j=1:angles
    theta = thetas(j);
    %field solution
    for q=1:n
        b(q) = exp(1i*(k*cos(theta)*XC(q)+k*sin(theta)*YC(q)));
    end
    u = F\b;
    %compute and store scattered field
    for p=1:m
        for q=1:n
            rho = sqrt((XC(q)-xm(p))^2+(YC(q)-ym(p))^2);
            s((j-1)*m+p) = s((j-1)*m+p) + rho*(eps_mat(q)-1)*u(q);
        end
    end
end
s = s * (-1i*pi*besselj(1,k*dr)*k/2); %multiply the constant
save('scatfield.mat','s');