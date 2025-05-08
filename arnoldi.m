function e = arnoldi(A, it)
    % ARNOLDI: Voert 'it' Arnoldi-iteraties uit op matrix A.
    % Geeft de eigenwaarden (Ritz-waarden) van de Hessenbergmatrix terug.
    
    n = size(A, 1);
    
    % Startvector: willekeurig maar genormaliseerd
    b = randn(n,1);
    b = b / norm(b);
    
    Q = zeros(n, it);        % Orthonormale basis
    H = zeros(it, it);       % Hessenbergmatrix
    
    Q(:,1) = b;
    
    for k = 1:it
        v = A * Q(:,k);
        
        % Gram-Schmidt orthogonalisatie
        for j = 1:k
            H(j,k) = Q(:,j)' * v;
            v = v - H(j,k) * Q(:,j);
        end
        
        if k < it
            H(k+1,k) = norm(v);
            if H(k+1,k) == 0
                break; % Krylov-ruimte is opgesloten
            end
            Q(:,k+1) = v / H(k+1,k);
        end
    end
    
    % Bereken Ritz-waarden (benaderde eigenwaarden)
    e = eig(H(1:it,1:it));
end
