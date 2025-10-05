function A = projectionMatrix(Th)

    x = Th.x(:); % vector columna de nodos [x0; x1; ...; xn]
    n = numel(x) - 1; % n + 1 nodos = n elementos 
    A = sparse(n+1, n+1);

    for k = 1:n
        h = x(k+1) - x(k);

        A(k,  k)   = A(k,  k)   + 2*h/6; % diagonal
        A(k,  k+1) = A(k,  k+1) +   h/6; % fuera diagonal
        A(k+1,k)   = A(k+1,k)   +   h/6; % fuera diagonal
        A(k+1,k+1) = A(k+1,k+1) + 2*h/6; % diagonal
    end
end