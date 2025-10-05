function errL2 = projectionError(Th, f, c, Q)
    x = Th.x(:);
    n = numel(x) - 1;

    a0 = Q.a0   ; b0 = Q.b0;
    xi = Q.x(:) ; w  = Q.w(:);

    scale = 1/(b0 - a0);

    acc = 0.0;

    for k = 1:n
        h = x(k+1) - x(k);
        J = h * scale;

        for q = 1:numel(xi)
            xq = x(k) + (xi(q) - a0) * (h * scale);
            phi1 = (x(k+1) - xq)/h;
            phi2 = (xq - x(k))/h;
            uh   = c(k)*phi1 + c(k+1)*phi2; % u_h en xq
            diff = f(xq) - uh;

            acc = acc + w(q)*J*(diff^2);
        end
    end

    errL2 = sqrt(acc);
end
