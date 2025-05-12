function c = implicitqr(A)
    [n,m] = size(A);
    assert(n==m);
    c = [];
    if n == 1
        c = [A(1, 1)];
    elseif n == 2
        trace = A(1,1) + A(2,2);
        det = A(1,1)*A(2,2) - A(1,2)*A(2,1);
        discrim = trace^2 - 4*det;
        c = [(trace + sqrt(discrim)) / 2, (trace - sqrt(discrim)) / 2];
    else
        H = hess(A);
        while true
            rho = H(n,n); %verschuiving
            H_alt = H-rho*eye(n);
            Q = eye(n);
            
            %Q via householder
            for k=1:n-1
                x = H_alt(k:end, k);
                e1 = zeros(length(x), 1); e1(1) = 1;
                v = sign(x(1))*norm(x)*e1 + x;
                v = v/norm(v);
                Hk = eye(n);
                Hk(k:end, k:end) = Hk(k:end, k:end)-2*(v*v');
                H_alt = Hk*H_alt;
                Q = Q*Hk';
            end

            H = Q'*H*Q;
            idx = find(abs(diag(H, -1)) < 1e-12);
            if ~isempty(idx)
                i = idx(1);
                c = [c, implicitqr(H(1:i, 1:i)), implicitqr(H(i+1:n, i+1:n))];
                break
            end
        end
    end
end