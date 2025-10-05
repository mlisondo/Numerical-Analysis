function b = linearformL(Th, fFunc, Q)

    x = Th.x(:);
    n = numel(x) - 1;
    b_full = zeros(n+1,1);

    a0 = Q.a0; b0 = Q.b0;
    xi = Q.x(:); w = Q.w(:);
    scale = 1/(b0 - a0);

    for k = 1:n
        h  = x(k+1) - x(k);
        J  = h * scale;

        be1 = 0.0; 
        be2 = 0.0;

        for q = 1:numel(xi)
            xq  = x(k) + (xi(q) - a0) * J;
            wqJ = w(q) * J;

            fq  = fFunc(xq);

            phi1 = (x(k+1) - xq)/h;
            phi2 = (xq - x(k))/h;

            be1 = be1 + wqJ * fq * phi1;
            be2 = be2 + wqJ * fq * phi2;
        end

        b_full(k)   = b_full(k)   + be1;
        b_full(k+1) = b_full(k+1) + be2;
    end

    b = b_full(2:end-1);
end
