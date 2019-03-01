N = 10; % Numero de puntos en la malla
L = 1e-10; % Tamaño en metros del dominio
q = -1;
dx = L/N;
periodic = 1;

% Cargo una matriz tridiagonal con coeficientes 1 -2 y 1
A = full(gallery("tridiag", N, 1, -2, 1));

% El rango de A es N-1, así que necesito una ecuación adicional
disp("Rank of A:")
disp(rank(A))

b = zeros(N,1);

% Coloco un electron en el medio
% b(N/2) = 1;
b(1/2 * N) = -q*dx^2;
%b(3/8 * N) = -q*dx^2;
%b(5/8 * N) = -q*dx^2;

% La suma de cargas debe ser cero
%b = b - 1/(N-1);


if(periodic)
	disp("PERIODIC")
	% Añado los coeficientes en la frontera, para que sea un dominio periódico
	A(1,N) = 1;
	A(N,1) = 1;
	%A(1,1) = -1;

	% Añado una equación para fijar el valor de $\phi_1 = 0$
	A = [[1, zeros(1,N-1)]; A];

	%A(1,:) = 0;
	%A(1,1) = 1;

	disp("Rank of A after adding the extra equation:")
	disp(rank(A))

	% Añado el lado derecho de la ecuación extra $\phi_1 = 0$
	b = [0; b];
endif

% Resuelvo el potencial electrico
phi = A\b;

% Compruebo que la solucion funciona
disp("Error:")
disp(norm(A*phi - b))

% Y aproximo el campo eléctrico $E = - \nabla \phi$ como $\phi_{i-1} - \phi_{i+1}$
E = shift(phi, -1) - shift(phi, +1);

x = linspace(0, L, N);

deltax = (N/2)*dx - x;

% El valor teórico debería ser este:
Et = q ./ deltax.^2;
phit = q ./ deltax;


plot(
	x, phi/10,	"-;Computed phi;"
	,x, E,		"-;Computed E;"
%	,x, phit/1e12,	";Theorical phi;"
%	,x, Et/1e20,		";Theorical E;"
)

grid()

%print -dpng -r200 field.png
