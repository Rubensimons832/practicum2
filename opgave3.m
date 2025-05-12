n = 1000;
A = sprand(n, n, 0.01);
it_max = 100;

ritz_vals = zeros(it_max, it_max);

for k = 1:it_max
    e = arnoldi(A, k);
    ritz_vals(1:length(e), k) = real(e);  % Bewaar reëel deel
end

% Exacte eigenwaarden
eigA = eig(full(A));
eigA_real = real(eigA);

% Plot in lineaire schaal
figure;
hold on;

% Dummy markers voor legend
h1 = plot(NaN, NaN, 'k.', 'MarkerSize', 8);
h2 = plot(NaN, NaN, 'r.', 'MarkerSize', 10);

% Plot Ritz-waarden
for i = 1:it_max
    vals = ritz_vals(1:i, i);
    plot(i * ones(size(vals)), vals, 'k.');
end

% Plot exacte eigenwaarden op stap 101
plot(101 * ones(size(eigA_real)), eigA_real, 'r.', 'MarkerSize', 10);

xlabel('Iteratiestap');
ylabel('Ritz-waarden (reëel deel)');
title('Convergentie van Ritz-waarden');
xlim([0 105]);
grid on;
legend([h1 h2], {'Ritz-waarden', 'Exacte eigenwaarden'}, 'Location', 'best');
