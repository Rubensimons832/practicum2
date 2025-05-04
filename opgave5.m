clear;
N = 50;
X = randn(N,N);
D = diag(linspace(N, 1, N));
A = X*D*inv(X);
disp(implicitqr(A));