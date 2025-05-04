function x = my_gmres1(A, b, it)
    n = length(b);
    Q = zeros(n, it+1);
    H = zeros(it+1, it);

    Q(:,1) = b / norm(b);

    for j=1:it
        v = A*Q(:,j);
        for i=1:j
            H(i,j) = Q(:,i)' * v;
            v = v - H(i,j) * Q(:,i);
        end
        H(j+1,j) = norm(v);
        if H(j+1,j) == 0
            break;
        end
        Q(:,j+1) = v / H(j+1,j);
    end
    y = H\(norm(b)*eye(it+1, 1));
    x = Q(:,1:it)*y; 
end

