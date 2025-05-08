% Opgave 3 (uitgebreid): Visualiseer Ritz-waarden én exacte eigenwaarden
n = 1000;
A = sprand(n, n, 0.01);
it_max = 100;

ritz_vals = zeros(it_max, it_max);

for k = 1:it_max
    e = arnoldi(A, k);
    ritz_vals(1:length(e), k) = sort(real(e));
end

% Plot convergentie Ritz-waarden
figure;
hold on;
for i = 1:it_max
    vals = ritz_vals(1:i, i);
    plot(i * ones(size(vals)), vals, 'k.');
end

% Bereken exacte eigenwaarden (volledige matrix)
eigA = eig(full(A));
eigA_real = sort(real(eigA));  % Enkel het reële deel

% Voeg exacte eigenwaarden toe als rode puntjes op stap 101
plot(101 * ones(size(eigA_real)), eigA_real, 'r.', 'MarkerSize', 10);

xlabel('Iteratiestap');
ylabel('Ritz-waarden (reëel deel)');
title('Convergentie van Ritz-waarden en exacte eigenwaarden');
xlim([0 105]);
grid on;
legend('Ritz-waarden', 'Exacte eigenwaarden', 'Location', 'best');
