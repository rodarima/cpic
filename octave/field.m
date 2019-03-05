N = 22; % Numero de puntos en la malla
L = 1; % Tamaño en metros del dominio
q = -1;
H = L/N;
NC = 1;
x0 = L/2;
periodic = 0;
extra = 1;
invertida = 0;
reemplazo = 0;
dirichlet = 0;
gaussian = 0;
magic = 1;

% Cargo una matriz tridiagonal con coeficientes 1 -2 y 1
A = full(gallery("tridiag", N, 1, -2, 1));

% El rango de A es N-1, así que necesito una ecuación adicional
disp("Rank of A:")
disp(rank(A))

b = zeros(N,1) - NC/(N-1)*abs(q)*H^2;


x = linspace(0, L, N);
r = abs(x0 - x);

% Coloco un electron en el medio
% b(N/2) = 1;
b(floor(N/2)) = -(q) * H^2;
%b(floor(N/4)) = -(q) * H^2;
%b(1/2 * N + 1) = -q/2*H^2;
%b(1/2 * N+1) = -q*H^2/2;
%b(3/8 * N) = -q*H^2;
%b(5/8 * N) = -q*H^2;

% La suma de cargas debe ser cero
%b = b - 1/(N-1);

%b(1) = 10;
%b(N) = b(1);

if(gaussian)

	b = exp(-((r./(L./30)).^2)./2).';
endif

if(magic)
	b = 2*q./(r.^2 .* abs(r));
	b = b.';
endif

if(periodic)
	disp("PERIODIC")
	% Añado los coeficientes en la frontera, para que sea un dominio periódico
	%A(1,N) = 1;
	%A(N,1) = 1;

	%A(1,1) = -1;

	A(1,:) = 0;
	A(1,1) = -2;
	A(1,2) = 1;

	A(N-1,:) = 0;
	A(N-1,N-2) = 1;
	A(N-1,N-1) = -2;

	A(N,:) = 0;
	A(N,1) = 1;
	A(N,N-1) = 1;

	% Añado una equación para fijar el valor de $\phi_1 = 0$
	%A = [[1, zeros(1,N-1)]; A];

	%A(1,:) = 0;
	%A(1,1) = 1;

	disp("Rank of A after adding the extra equation:")
	disp(rank(A))

	% Añado el lado derecho de la ecuación extra $\phi_1 = 0$
	%b = [0; b];
endif

if(extra)
	% Añado una equación para fijar el valor de $\phi_1 = 0$
	A = [[1, zeros(1,N-1)]; A];
	%A = [A; [zeros(1,N-1), 1]];

	% Añado el lado derecho de la ecuación extra $\phi_1 = 0$
	b = [0; b];
	%b = [0; b; 0];
endif

if(reemplazo)
	% Añado una equación para fijar el valor de $\phi_1 = 0$
	A(1,:) = [1, zeros(1,N-1)];
	A(N,:) = [zeros(1,N-1), 1];
	A(2,1) = 0;
	A(N-1,N) = 0;

	% Añado el lado derecho de la ecuación extra $\phi_1 = 0$
	b = zeros(N,1) - 1/(N-3)*abs(q)*H^2;
	b([1,N]) = 0;
	b(floor(N/2)) = -q*H^2;
endif


if(dirichlet)
	A(1,1) = -1;
	b(1) = b(1) / 2;


endif


%A = A./H^2;

alpha = 0.01
%A = A.*alpha
%b = b.*alpha

% Resuelvo el potencial electrico
phi = A\b;

% Compruebo que la solucion funciona
disp("Error:")
disp(norm(A*phi - b))

% Y aproximo el campo eléctrico $E = - \nabla \phi$ como $\phi_{i-1} - \phi_{i+1}$
E = (shift(phi, -1) - shift(phi, +1))/H;

% El valor teórico debería ser este:
Et = q ./ (abs(r) .* r);
phit = q ./ abs(r);


plot(
	x, b(2:N+1)/max(abs(b(2:N+1))),	"-;rho;",
	x, phi/max(abs(phi)),		"-;Computed phi;"
	,x, E/max(abs(E)),	"-;Computed E;"
	,x, phit/max(abs(phit)),	";Theorical phi;"
	,x, Et/max(abs(Et)),	";Theorical E;"
)

grid()

%print -dpng -r200 field.png
