A = randn(10);   % Kleine matrix
ritz = arnoldi(A, 10); 
eigA = eig(A);   % Exacte eigenwaarden

% Vergelijk met echte eigenwaarden (visueel)
disp([sort(real(ritz)), sort(real(eigA))]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

A = randn(1000); %+ 1i * randn(1000);   % Complexe matrix
ritz = arnoldi(A, 100);              % Benaderde eigenwaarden
eigA = eig(A);                      % Exacte eigenwaarden

% Plot in het complexe vlak
figure;
plot(real(eigA), imag(eigA), 'ro', 'DisplayName', 'Exact'); hold on;
plot(real(ritz), imag(ritz), 'kx', 'MarkerSize', 8, 'LineWidth', 2, 'DisplayName', 'Ritz');
legend('Location','best');
xlabel('Reëel deel');
ylabel('Imaginair deel');
title('Complexe vlak: exacte vs. Ritz-waarden');
axis equal; grid on;
