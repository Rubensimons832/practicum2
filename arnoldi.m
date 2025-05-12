

function ritzwaarden = arnoldi(A,k)
    % ARNOLDI: Voert 'k' Arnoldi-iteraties uit op matrix A.
    % Geeft de eigenwaarden (Ritz-waarden) van de Hessenbergmatrix terug.
    
    % INITIALISATIE: k<<m (of k<m)
    n = size(A, 1); % A € R^(nxm)
    b = randn(n,1); % Startvector: willekeurig
    Q = zeros(n,k); % Orthonormale basis [q1, q2, ..., qk]
    H = zeros(k,k); % Hessenbergmatrix
    
    Q(:,1) = b/norm(b);                  % Algoritme 6: regel 1
    for j = 1:k                          % Algoritme 6: regel 2
        v = A * Q(:,j);                  % Algoritme 6: regel 3
        % Gram-Schmidt orthogonalisatie
        for i = 1:j                      % Algoritme 6: regel 4
            H(i,j) = Q(:,i)' * v;        % Algoritme 6: regel 5
            v = v - H(i,j) * Q(:,i);     % Algoritme 6: regel 6
        end                              % Algoritme 6: regel 7
        
        if j<k % Als het niet de laatste iteratie is.
            H(j+1,j) = norm(v);          % Algoritme 6: regel 8
            if H(j+1,j) == 0 % dan wordt Krylov-
                break;       % deelruimte niet meer
            end              % uitgebreid.(detail)
            Q(:,j+1) = v / H(j+1,j);     % Algoritme 6: regel 9
        end                              % Algoritme 6: regel 10
    end
    
    % Bereken Ritz-waarden (benaderde eigenwaarden)
    % mag volgens opgave gewoon met eig().
    ritzwaarden = eig(H(1:k,1:k));
end
