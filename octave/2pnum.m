N = 201; % Numero de puntos en la malla
L = 1; % Tamaño en metros del dominio
q = 1;
q1 = q;
q2 = -q;
H = L/N;
NC = 2;
xf = 2/10
x_ = L * xf;
periodic = 1;
extra = 1;
invertida = 0;
reemplazo = 0;
dirichlet = 0;
gaussian = 0;
magic = 1;

% Cargo una matriz tridiagonal con coeficientes 1 -2 y 1
A = full(gallery("tridiag", N, 1, -2, 1));

rho = zeros(N,1);

x = linspace(-L/2, L/2, N);

rho(floor(N/2 - xf*N)) = q1/H * H^2;
rho(floor(N/2 + xf*N)) = q2/H * H^2;


if(periodic)
	% Añado los coeficientes en la frontera, para que sea un dominio periódico
	A(1,N) = 1;
	A(N,1) = 1;
endif

if(extra)
	% Añado una equación para fijar el valor de $\phi_1 = 0$
	A = [[1, zeros(1,N-1)]; A];
	%A = [[1, zeros(1,N-1)]; A];
	%A = [A; [zeros(1,N-1), 1]];

	% Añado el lado derecho de la ecuación extra $\phi_1 = 0$
	rho = [0; rho];
	%rho = [0; rho; 0];
endif

% Resuelvo el potencial electrico
phi = A\rho;

E = (shift(phi, -1) - shift(phi, +1)) / H;

% Compruebo que la solucion funciona
disp("Error:")
disp(norm(A*phi - rho))





function r = th_phi(x, x_, L, q)

	r = zeros(size(x));

	i0 = abs(x) <= x_;
	i1 = x <= -x_;
	i2 = x >= x_;

	r(i0) = q ./ L .* (L - 2.*x_) .* x(i0);
	r(i1) = q ./ L .* (-L - 2.*x(i1)) .* x_;
	r(i2) = q ./ L .* (L - 2.*x(i2)) .* x_;

	return;
endfunction

function r = th_E(x, x_, L, q)

	o = ones(size(x));
	r = -q/(2*L) * (L - 4*x_ * o);

	return;
endfunction

phit = th_phi(x, x_, L, q);
Et = th_E(x, x_, L, q);







plot(
	x, phi,	"-;phi numeric;",
%	x, E,	"o;E numeric;",
	x, phit,"--;phi theoric;"
%	x, Et,	"-;E theoric;"
)

%axis([-L/2 L/2 -.25 .25]);
%axis([-L/2 L/2 -2 2]);
grid()

print -dpng -r200 field.png
