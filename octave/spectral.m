N = 8;
Nt = N * N;
Np = 1;
q = 1;
e0 = 1;
L = N;
xmin = -L/2;
xmax = L/2;
H = L/N;
H2 = H*H;
tol = 1e-10;
maxit = 100;
use_contour = 1;
d_term = 0.0; % Hockey book uses this so 1.0, but seems wrong
scale_arrows = 2.0;
wx = 1500;
wy = 500;
dx2 = L/N * 2

range = linspace(xmin, xmax, N);
[X, Y] = meshgrid(range, range);

% Add one charge exactly in the middle of the grid
rho = zeros(Nt, 1);
%rho = rho - q/(Nt-1);
%rho(2) = -q/(e0*H2);
%rho(12) = -q/(e0*H2);
%rho(Nt*1/4 + N*12/16) = -q/(e0*H2);
%rho(Nt*2/4 + N*4/16) = +q/(e0*H2);
%rho(Nt*3/4:Nt*3/4+N/2) = + q/(e0*H2*3);
%rho(Nt*3/4 + N*8/16) = -q/(e0*H2);

sum_rho = sum(rho)

%rho(rho == 0) = -sum_rho / (Nt - 3);

%rho = rho - sum_rho / Nt;

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



% ------------------------------------- MFT -------------------------------- %

rho_2d = reshape(rho, N, N);

rho_2d(1,2) = -q;
rho_2d(3,3) = -q;

g = fft2(rho_2d, N, N);

G = zeros(N, N);

for k = 0:N-1
	for l = 0:N-1
		G(k+1,l+1) = 1/(  2*(cos(2*pi*k/N) + cos(2*pi*l/N)) - 4.0 +
		d_term);
	endfor
endfor

% Correct division by 0
G(1,1)=0;

g
rho_2d

g = g .* G;
g

phi_fft = real(ifft2(g, N, N));

phi_fft

Ex = (shift(phi_fft, +1, 2) - shift(phi_fft, -1, 2)) ./ dx2
Ey = (shift(phi_fft, +1, 1) - shift(phi_fft, -1, 1)) ./ dx2


%return
% -------------------------------------------------------------------------- %


%A = [1,zeros(1, Nt-1);A];
%rho = [0;rho];
%A(1,:) = 0;
%A(1,1) = -r;

%phi = A\rho;
%phi = pcr(-A, -rho, tol, maxit);

[L,U] = lu(A);

y = linsolve(L, rho);
phi_linsolve = linsolve(U, y);




phi = phi_linsolve;

phi = reshape(phi, N, N);

%[Ex, Ey] = gradient(-phi, H);
Ex = (shift(phi, +1, 2) - shift(phi, -1, 2)) ./ dx2
Ey = (shift(phi, +1, 1) - shift(phi, -1, 1)) ./ dx2

subplot(2, 2, 1);
if(use_contour)
	contourf(X, Y, phi);
	colorbar()
else
	surf(X, Y, phi);
endif
title("Electric potential (\\phi) using linsolve");

subplot(2, 2, 2);
%contourf(X, Y, E);
%colorbar()
q=quiver(X, Y, Ex, Ey, scale_arrows);
%set(q, "showarrowhead", false);
set(q, "maxheadsize", 0.1);

title("Electric field (E)");




phi = phi_fft;

phi = reshape(phi, N, N);

%[Ex, Ey] = gradient(-phi, H);
Ex = (shift(phi, +1, 2) - shift(phi, -1, 2)) ./ dx2
Ey = (shift(phi, +1, 1) - shift(phi, -1, 1)) ./ dx2

subplot(2, 2, 3);
if(use_contour)
	contourf(X, Y, phi);
	colorbar()
else
	surf(X, Y, phi);
endif
title("Electric potential (\\phi) using MFT");

subplot(2, 2, 4);
%contourf(X, Y, E);
%colorbar()
q=quiver(X, Y, Ex, Ey, scale_arrows);
%set(q, "showarrowhead", false);
set(q, "maxheadsize", 0.1);

title("Electric field (E)");
