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
maxit = 1000;
use_contour = 1;

range = linspace(xmin, xmax, N);
[X, Y] = meshgrid(range, range);

% Add one charge exactly in the middle of the grid
rho = zeros(Nt, 1);
%rho = rho - q/(Nt-1);
rho(Nt/2 + N*2/8) = -q/(e0*H2);
rho(Nt/2 + N*6/8) = -q/(e0*H2);

A = gallery("tridiag", N, 1, -4, 1);

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

A = [1,zeros(1, Nt-1);A];
rho = [0;rho];
%A(1,:) = 0;
%A(1,1) = -r;

phi = A\rho;
%phi = pcg(-A, -rho, tol, maxit);

phi = reshape(phi, N, N);

E = -gradient(phi);

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
if(use_contour)
	contourf(X, Y, E);
	colorbar()
else
	surf(X, Y, E);
endif
title("Electric field (E)");

print('contour2d.png', '-dpng', '-S1024,300');
