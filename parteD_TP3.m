function parteD_TP3()
    eps_list = [1e-2, 1e-3, 1e-10];
    h_list   = [1e-1, 1e-2, 1e-3];

    [xi, wi] = gl10();
    Q.a0 = -1; Q.b0 = +1; Q.x = xi(:); Q.w = wi(:);

    for m = 1:numel(eps_list)
        eps  = eps_list(m);
        aFunc = @(x) eps;
        bFunc = @(x) 1.0;
        fFunc = @(x) 1.0;

        uex = @(x) 1 - cosh(x./sqrt(eps)) ./ cosh(1./sqrt(eps));

        figure; hold on; grid on; box on;

        xfine = linspace(-1,1,4001).';
        plot(xfine, uex(xfine), 'k-', 'LineWidth', 1.8, 'DisplayName', '"Exacta"');

        errs = zeros(numel(h_list),1);

        for j = 1:numel(h_list)
            h  = h_list(j);
            n  = round(2/h);

            gamma = 3;                              % ajusta: 2–6 va bien
            Th.x  = mesh_tanh(n, gamma);
            
            % Th.x = linspace(-1,1,n+1).';

            A  = bilinearFormA(Th, aFunc, bFunc, Q);
            b  = linearformL(Th, fFunc, Q);
            uh = [0; A \ b; 0]; % Dirichlet

            uex = @(x) 1 - cosh(x./sqrt(eps)) ./ cosh(1./sqrt(eps));
            errs(j) = L2error_GL10(Th, uh, uex, Q);

            % Curva
            plot(Th.x, uh, '-', 'DisplayName', sprintf('h=%.3g', h), ...
                 'LineWidth', 1.2, 'MarkerSize', 2);
        end

        xlabel('x'); ylabel('u_h');
        legend('Location','best'); box on;

        fprintf('epsilon = %.3g\n', eps);
        for j = 1:numel(h_list)
            fprintf('  h = %-8.3g   ||u - u_h||_{L2} = %.6e\n', h_list(j), errs(j));
        end
    end
end

function x = mesh_tanh(n, gamma)
    xi = linspace(-1, 1, n+1);
    x  = tanh(gamma*xi) ./ tanh(gamma);
    x(1) = -1; x(end) = 1;
    x = x(:);
end


function [x,w] = gl10()
    x = [ -0.9739065285171717
          -0.8650633666889845
          -0.6794095682990244
          -0.4338837391175581
          -0.1488743389816312
           0.1488743389816312
           0.4338837391175581
           0.6794095682990244
           0.8650633666889845
           0.9739065285171717 ];
    w = [ 0.0666713443086881
          0.1494513491505806
          0.2190863625159820
          0.2692667193099963
          0.2955242247147529
          0.2955242247147529
          0.2692667193099963
          0.2190863625159820
          0.1494513491505806
          0.0666713443086881 ];
end

function E = L2error_GL10(Th, uh, u_exact, Q)
    x = Th.x; n = numel(x)-1;
    a0 = Q.a0; b0 = Q.b0; xi = Q.x; w = Q.w; scale = 1/(b0-a0);
    err2 = 0.0;
    for e = 1:n
        h  = x(e+1) - x(e);
        J  = h * scale;
        for q = 1:numel(xi)
            xq   = x(e) + (xi(q) - a0) * J;
            phi1 = (x(e+1) - xq)/h;
            phi2 = (xq - x(e))/h;
            uhq  = uh(e)*phi1 + uh(e+1)*phi2;
            diff = u_exact(xq) - uhq;
            err2 = err2 + w(q) * J * (diff.^2);
        end
    end
    E = sqrt(err2);
end