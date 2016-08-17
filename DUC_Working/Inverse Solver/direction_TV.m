function [d_n] = direction_TV(w, chi, g_n_chi, g_prev, d_prev)

d_n = zeros(N^2,1);
for i = 1:N^2
	numer = dot(g_n_chi, g_n_chi - g_prev);
	denom = dot(g_prev, g_prev);

	d_n(i) = g_n_chi(i) + (Re(numer)./denom).*d_prev;
end