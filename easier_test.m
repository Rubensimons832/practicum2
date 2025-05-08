% MATRIX A met rang 4 en dominante eigenwaarde
n = 100;
D = zeros(n,1);
D(1:4) = [10, 5, 1, -1];  % Slechts 4 niet-nul eigenwaarden

V = orth(randn(n));  % Orthogonale basis
A = V * diag(D) * inv(V);  % Diagonaliseerbare matrix met rang 4

% ARNOLDI met slechts 2 iteraties
it = 3;
ritz = arnoldi(A, it);

% Exacte eigenwaarden
eigA = eig(A);
eigA_real = sort(real(eigA));
ritz_real = sort(real(ritz));

% PLOT
figure;
hold on;

% Exacte eigenwaarden in rood (alle 4 niet-nul, rest nul)
plot(1.1 * ones(size(eigA_real)), eigA_real, 'r.', 'MarkerSize', 15);

% Ritz-waarden in zwart (rechts van x=1)
plot(1.0 * ones(size(ritz_real)), ritz_real, 'k.', 'MarkerSize', 15);

xlim([0.9 1.2]);
ylim([-2 12]);
xlabel('Index');
ylabel('Waarde');
title('Vergelijking Ritz-waarden (zwart) en eigenwaarden (rood)');
legend('Exacte eigenwaarden', 'Ritz-waarden (2 iteraties)', 'Location', 'best');
grid on;
