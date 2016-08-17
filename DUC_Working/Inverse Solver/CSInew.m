% Data input;
% Paper followed is CSI_State of Heart

green_measureXcells= aux_GS; 	% cell refers to the portion in which the contrast is assumed to be piecewise constant; i,j entry of this variable corresponds to the integral of green's function on j^th cell by assuming the source to be at i^th measurement point  
green_cellsXcells= aux_GD;        	% i,j entry of this variable corresponds to the integral of green's function on j^th cell by assuming the source to be at i^th cell
E_measured=conj(f); 		%  i,j entry of this variable corresponds to measured field on i^th measurement point corresponding to the the source located at j^th measurement angle

E_incident=conj(u_inc); 				%  (i,j) entry of this variable denotes the incident field at i^th cell corresponding to the source located at the j^th measurement angle
N_incid_angles=size(E_measured,2);		% no. of incident angles
N_cells=size(green_cellsXcells,1);		% no. of cells in which the scattere is divided

% The scheme described above is followed consistently in representation of various physical quantities like Electric field. The i,j entry of actual electric feld u would denote the same thing as E_incident above

% initializing w(contrast sources) through back propagation;

w=(1j)*ones(N_cells,N_incid_angles); 		
for j=1:N_incid_angles
    k1=green_measureXcells'*E_measured(:,j);
    n1=norm(k1)^2;
    d1=norm(green_measureXcells*k1)^2;
    w(:,j)=(n1/d1)*k1;						% eq (35) of the CSI paper
end

% initializing u(actual fields) inside the cells
u=(1j)*ones(N_cells,N_incid_angles);
for j=1:N_incid_angles
    u(:,j)=E_incident(:,j)+green_cellsXcells*w(:,j);				%simple integral source eq or eq 36 of the CSI paper 
end

% initializing chi (contrast) inside the cell. Note that this would be a vector as it is assumed to be source independent.
n1=zeros(N_cells,1);
d1=0;
for j=1:N_incid_angles
    n1=n1+w(:,j).*conj(u(:,j));		
    d1=d1+abs(u(:,j)).^2;
end
chi=n1./d1;			% eq 38 of CSI paper

for i = 1 : N_cells
    if(real(chi(i))<0)
        				% since real part of chi can't be -ve, you can enforce this condition here
    end
   
end

%initializing the constant etaS which is the weight of the data error in CSI cost function. It is defined in eq 25 of the paper

d1=0;
for j=1:N_incid_angles
    d1=d1+norm(E_measured(:,j))^2;		% eq 25
end
etaS=1/d1;	


% From now on each symbol used in the CSI paper has the same variable name here. So, I would just refer which equation of the paper is the line referring to 

iter=1; % iteration count

rho=E_measured-green_measureXcells*w;   						%eq 26
r = C_Mat_element(chi,E_incident) - w + C_Mat_element(chi,green_cellsXcells*w);		%eq27
v=zeros(N_cells,N_incid_angles);	%initialize variable v
g_prev=v; 			% initialize the direction g_(n-1) corresponding to previous iteration
cost = inf;			% initializing cost
disp('Entering Iterations')
while (iter<=10000)
    
    d1=0;
    for j=1:N_incid_angles
        d1=d1+norm(chi.*E_incident(:,j))^2;	
    end
    etaD=1/d1;			%eq25
    
    
    temp=C_Mat_element(conj(chi),r);
    g_curr=-etaS*green_measureXcells'*rho - etaD*(r-green_cellsXcells'*temp);		%eq30
    
    temp=g_curr-g_prev;
    n1=0;d1=0;
    for j=1:N_incid_angles
        n1=n1+real(g_curr(:,j)'*temp(:,j));      % big doubt here: inn_prod = conj trans or just trans??; doubt resolved should be conjugate transposed
        d1=d1+g_prev(:,j)'*g_prev(:,j);   
    end
    if (d1==0)
        d1=1;     % only for 1st iteration since d and v are both 0
    end
    v=g_curr+(n1/d1)*v;			%eq29
    
    n1=0;
    d1=0;
    temp=v-C_Mat_element(chi,green_cellsXcells*v);
    for j=1:N_incid_angles
        n1=n1-real(g_curr(:,j)'*v(:,j));      % big doubt here: inn_prod = conj trans or just trans??
        d1=d1+etaS*(norm(green_measureXcells*v(:,j))^2)  +  etaD*(norm(temp(:,j))^2);   
    end
    alpha=n1/d1;		%eq34
    
    
    w=w+1*alpha*v;			%eq 28
    
    % Now update u and chi
    u=E_incident+green_cellsXcells*w;			%eq36
    
    n1=0;d1=0;
    for j=1:N_incid_angles
        n1=n1+w(:,j).*conj(u(:,j));
        d1=d1+abs(u(:,j)).^2;
    end
    chiold = chi;
    chi=n1./d1;		% eq 38
    
for i = 1 : N_cells
        if(real(chi(i))<0)
        			% again decide what to do here
        end
        chi(i)=real(chi(i));	% enforcing the constraint that chi is real
    end
    
    
    g_prev=g_curr;		% updation current direction to be the previous direction (For the next iteration current dirction should be the previous direction)
 

% from now on the code is not neccessary....
   
    temp=C_Mat_element(chi,E_incident)- w + C_Mat_element(chi,green_cellsXcells*w);
    costinit = cost;
    cost=0;
    for j=1:N_incid_angles
        cost=cost+ etaS*(norm(E_measured(:,j)-green_measureXcells*w(:,j))^2)+   etaD*(norm(temp(:,j))^2);		% calculating the cost
    end
    disp(cost);
    if(cost > costinit)
        chi = chiold;				% artificial measure adopted by us so as not to go in the direction increasing cost. Would not be needed if thCSI always gives correct directions. So, try to remove this part and just try to carry on a large no. of iterations of CSI without checking anything
       break;
    end

% end of unnecessary part
    chi_temp = reshape(chi, [N,N]);    
    surf(x_temp, y_temp, real(chi_temp)); 
    
    %Graphically see chi on the mesh
    pause;



    iter=iter+1;


end