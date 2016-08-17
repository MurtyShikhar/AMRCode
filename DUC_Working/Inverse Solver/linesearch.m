% does line search to calculate alpha_n in equation (33)
% alpha_n is of the form (gamma)^s where s is (0,1,2..)
function[alpha] = linesearch(grad_fn, v_n, a,b)

global M;

alpha = 1;
rho = 0.5;
c = 1e-4;


fval = cost_fn(a,b);
gval  = 0;
grad_fn = -1.*(grad_fn);
for i = 1:M
	gval = gval + dot(grad_fn(:,i), v_n(:,i));
end

while (cost_fn(a + alpha.*v_n, b)  < fval + c*alpha*gval)
	alpha = alpha*rho;
end


