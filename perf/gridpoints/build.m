% Number of maximum grid points
%Nmax=1073741824;
Nmax=8192**2;

% Intervals
k = 16;

% Number of plasma chunks
b = 128;

i = 1:k;

x = sqrt(Nmax .* i / k);

r = round(x / b) * b;

printf('%d\n', r);
