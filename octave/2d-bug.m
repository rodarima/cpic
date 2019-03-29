N = 5;
Nt = N * N;
Np = 1;
q = -1e2;
e0 = 1;
L = N;
xmin = -L/2;
xmax = L/2;
H = L/N;
H2 = H*H;
tol = 1e-10;
maxit = 100;
use_contour = 1;
scale_arrows = 2.0;
wx = 1500;
wy = 500;

range = linspace(xmin, xmax, N);
[X, Y] = meshgrid(range, range);

rho = zeros(N, N);
%rho(1, 3) = -q/(e0*H2);
rho(5, 3) = -q/(e0*H2);
%rho(3, 1) = +4.26e-02;

rho


rho = reshape(rho, Nt, 1);

sum_rho = sum(rho);

rho = rho - sum_rho / Nt;

rho_sq = reshape(rho, N, N);

printf('Sum rho = %e\n', sum(rho));
rho_sq

% Big matrix of coefficients
A = zeros(Nt, Nt);

function k = delta_index(i, j, di, dj, n, m)
	ni = mod(i+di-1, n);
	nj = mod(j+dj-1, m);

	k = ni * n + nj + 1;
endfunction

for i = 1:N
	for j = 1:N
		k = (i - 1) * N + j;
%		r = 1 + mod(j+1, Nt);
%		l = 1 + mod(j-1, Nt);
%		t = 1 + mod(j-N, Nt);
%		b = 1 + mod(j+N, Nt);

		r = delta_index(i, j, 0,+1, N, N);
		l = delta_index(i, j, 0,-1, N, N);
		t = delta_index(i, j,-1, 0, N, N);
		b = delta_index(i, j,+1, 0, N, N);

		printf("For ij=(%d,%d) k=%d, r %d, l %d, t %d b %d\n", i, j, k, r, l, t, b);

		A(k,k) = -4;
		A(k,r) = 1;
		A(k,l) = 1;
		A(k,t) = 1;
		A(k,b) = 1;
	endfor
endfor

A(1,1) = -3;

A

%A = [1,zeros(1, Nt-1);A];
%rho = [0;rho];
%A(1,:) = 0;
%A(1,1) = -r;

%phi = A\rho;
%phi = pcr(-A, -rho, tol, maxit);

[L,U,P] = lu(A);

L
U
P

y = linsolve(L, rho);
phi2 = linsolve(U, y);

phi = phi2;

phi = reshape(phi, N, N);

phi

[Ex, Ey] = gradient(-phi, H);

graphics_toolkit qt;
%f = figure('visible','off');

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
%print('contour2d.png', '-dpng', size_str);
