function b = projectionVector(Th, f, Q)
    x = Th.x(:);
    n = numel(x) - 1;
    b = zeros(n+1,1);

    % intervalo de referencia
    a0 = Q.a0; b0 = Q.b0;

    % puntos y pesos de cuadratura
    xi = Q.x(:); w = Q.w(:);

    scale = 1/(b0 - a0); % para el Jacobiano

    for k = 1:n
        h = x(k+1) - x(k);
        J = h * scale; % dx/dxi

        for q = 1:numel(xi)
            xq = x(k) + (xi(q) - a0) * (h * scale); % punto físico
            phi1 = (x(k+1) - xq)/h;   % asociada a nodo k
            phi2 = (xq - x(k))/h;     % asociada a nodo k+1
            fq   = f(xq);

            b(k)   = b(k)   + w(q)*J*fq*phi1;
            b(k+1) = b(k+1) + w(q)*J*fq*phi2;
        end
    end
end
