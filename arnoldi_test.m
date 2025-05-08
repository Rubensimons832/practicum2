A = randn(10);   % Kleine matrix
ritz = arnoldi(A, 10); 
eigA = eig(A);   % Exacte eigenwaarden

% Vergelijk met echte eigenwaarden (visueel)
disp([sort(real(ritz)), sort(real(eigA))]);
