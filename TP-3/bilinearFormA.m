function A = bilinearFormA(Th, aFunc, bFunc, Q)

    x = Th.x(:);
    n = numel(x) - 1;
    A_full = sparse(n+1, n+1);

    a0 = Q.a0; b0 = Q.b0;
    xi = Q.x(:); w = Q.w(:);
    scale = 1/(b0 - a0);

    for k = 1:n
        h  = x(k+1) - x(k);
        J  = h * scale;
        dphi1 = -1/h; % dphi2 = -dphi1

        Ke11 = 0.0; Ke12 = 0.0; Ke22 = 0.0;

        for q = 1:numel(xi)
            xq  = x(k) + (xi(q) - a0) * J;
            wqJ = w(q) * J;

            alpha = aFunc(xq);
            beta  = bFunc(xq);

            phi1 = (x(k+1) - xq)/h;
            phi2 = (xq - x(k))/h;

            Ke11 = Ke11 + wqJ * (alpha*dphi1*dphi1 + beta*phi1*phi1);
            Ke12 = Ke12 + wqJ * (-alpha*dphi1*dphi1 + beta*phi1*phi2);
            Ke22 = Ke22 + wqJ * (alpha*dphi1*dphi1 + beta*phi2*phi2);
        end
        Ke21 = Ke12;

        A_full(k,  k)   = A_full(k,  k)   + Ke11;
        A_full(k,  k+1) = A_full(k,  k+1) + Ke12;
        A_full(k+1,k)   = A_full(k+1,k)   + Ke21;
        A_full(k+1,k+1) = A_full(k+1,k+1) + Ke22;
    end
    A = A_full(2:end-1, 2:end-1);
end
