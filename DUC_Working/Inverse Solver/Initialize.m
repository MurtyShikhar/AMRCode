function [v_0] = Initialize(N)

% all globals
global f;   %M x M matrix
global aux_GD; % auxilliary GD N^2 x N^2
global aux_GS; % auxilliary GS M x N^2
global aux_GS_star;
global w_0;  % initialize the contrast source matrix N^2 x M
global chi_0; % initialize the contrast matrix N^2 x 1
global eta_S; % is a scalar
global eta_D; % is a scalar that depends on the current iteration
global u_inc; % is an N^2 x M matrix
global M; % number of sources/ receivers

    numerator=(aux_GS_star*f);
    denominator=(aux_GS*numerator);
    w_0=zeros(N^2,M);
    for j=1:M
        %norm_2(numerator(:,j))
        %pause;
        w_0(:,j)=((norm(numerator(:,j))^2)/norm(denominator(:,j))^2).*numerator(:,j);
    end;
    
    %%%%%---eqn (25) first part----%%%%%
    
    eta_S=0;
  for j=1:M
      eta_S=eta_S+(norm(f(:,j)))^2;
  end

  eta_S=1/eta_S;
  
  %%%%%%----this part calculates v_0----%%%%%
  u_0=u_inc+aux_GD*w_0;

for i=1:N*N
  numer = 0;
  denom = 0;
  for j = 1:M
    numer = numer + w_0(i,j)*conj(u_0(i,j));
    denom = denom + norm(u_0(i,j))^2;
  end
  chi_0(i,1) = numer/denom;
end
  v_0=zeros(N^2,M);
  
  %%%%%---eqn (25) second part----%%%%%
  
  aux_eta_D=0;
  for j=1:M
      aux_eta_D=aux_eta_D+(norm(chi_0.*u_inc(:,j)))^2;
  end
  eta_D=1/aux_eta_D;   