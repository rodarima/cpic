N = 32;
Nt = N * N;
Np = 1;
q = -1;
e0 = 1;
L = N;
xmin = -L/2;
xmax = L/2;
H = 1/N;
H2 = H*H;
tol = 1e-10;
maxit = 100;
use_contour = 1;
scale_arrows = 2.0;
wx = 1500;
wy = 500;

range = linspace(xmin, xmax, N);
[X, Y] = meshgrid(range, range);

% Add one charge exactly in the middle of the grid
rho = zeros(Nt, 1);
%rho = rho - q/(Nt-1);
rho(Nt*1/4 + N*12/16) = -q/(e0*H2);
rho(Nt*2/4 + N*4/16) = +q/(e0*H2);
%rho(Nt*3/4:Nt*3/4+N/2) = + q/(e0*H2*3);
rho(Nt*3/4 + N*8/16) = -q/(e0*H2);

sum_rho = sum(rho)

rho(rho == 0) = -sum_rho / (Nt - 3);

sum(rho)

% Big matrix of coefficients
A = zeros(Nt, Nt);

for i = 1:Nt
	j = i - 1;
	r = 1 + mod(j+1, Nt);
	l = 1 + mod(j-1, Nt);
	t = 1 + mod(j-N, Nt);
	b = 1 + mod(j+N, Nt);

	A(i,i) = -4;
	A(i,r) = 1;
	A(i,l) = 1;
	A(i,t) = 1;
	A(i,b) = 1;
endfor

A(1,1) = -3;

%A = [1,zeros(1, Nt-1);A];
%rho = [0;rho];
%A(1,:) = 0;
%A(1,1) = -r;

%phi = A\rho;
%phi = pcr(-A, -rho, tol, maxit);

[L,U] = lu(A);

y = linsolve(L, rho);
phi2 = linsolve(U, y);

phi2 = phi;

phi = reshape(phi, N, N);

[Ex, Ey] = gradient(-phi, H);

graphics_toolkit qt;
f = figure('visible','off');

subplot(1, 2, 1);
if(use_contour)
	contourf(X, Y, phi);
	colorbar()
else
	surf(X, Y, phi);
endif
title("Electric potential (\\phi)");

subplot(1, 2, 2);
%contourf(X, Y, E);
%colorbar()
q=quiver(X, Y, Ex, Ey, scale_arrows);
%set(q, "showarrowhead", false);
set(q, "maxheadsize", 0.1);

title("Electric field (E)");

size_str = sprintf("-S%d,%d", wx, wy);
print('contour2d.png', '-dpng', size_str);
