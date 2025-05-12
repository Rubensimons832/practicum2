clear;
N = 100;
it = 30;

A = hilb(N);
x = randn(N,1);
b = A*x;

fout = zeros(it, 1);
for i = 1 : it
    fout(i) = norm(b - A*my_gmres1(A,b,i));
end

figure;
semilogy(1:it, fout); xlabel('iterations'); ylabel('||b - Ax^{(k)}||_2');