L = 1;
N = 1001;
x = linspace(-L/2, L/2, N);
x_ = 1/8;
q = -1;

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

E = (shift(phit, -1) - shift(phit, +1))./(L/N);

plot(
	x, phit
	, x, Et
	,x, E
	);
axis([-L/2 L/2 -1 1]);
grid();
%set(gca, 'xtick', linspace(-0.5, 0.5, 9));
