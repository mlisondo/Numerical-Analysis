function errors = calc_error()
    format short e
    aFunc   = @(x) 1.0;
    bFunc   = @(x) 1.0;
    kmin = 4;
    kmax = 19;

    u_exact = @(x) 1e-4 .* x .* (1 - x) .* exp(15.*x);
    fFunc   = @(x) 1e-4 .* exp(15.*x) .* (224.*x.^2 - 164.*x - 28);

    [xi, wi] = gl10();
    Q.a0 = -1; Q.b0 = +1; Q.x = xi(:); Q.w = wi(:);

    kList  = (kmin:kmax).';
    errors = [kList, zeros(numel(kList),1)];

    for t = 1:numel(kList)
        k    = kList(t);
        n    = 2^k;
        Th.x = linspace(0,1,n+1).';

        A = bilinearFormA(Th, aFunc, bFunc, Q);
        b = linearformL(Th, fFunc, Q);

        uh = [0; A \ b; 0];       % Dirichlet: u(0)=u(1)=0

        errL2 = L2error_GL10(Th, uh, u_exact, Q);
        errors(t,2) = errL2;
    end

    alphas = NaN(size(kList));
    ratios = errors(1:end-1,2) ./ errors(2:end,2);
    alphas(1:end-1) = log(ratios) / log(2);

    for t = 1:numel(kList)
        if isnan(alphas(t))
            fprintf('%4d  % .5e    %s\n', kList(t), errors(t,2), '---');
        else
            fprintf('%4d  % .5e  % .6f\n', kList(t), errors(t,2), alphas(t));
        end
    end


    h = 2.^(-kList);
    figure;
    loglog(h, errors(:,2), '-o', 'LineWidth', 1.5, 'MarkerSize', 5, ...
           'DisplayName','error'); grid on; hold on;
    xlabel('$h$','Interpreter','latex','FontSize',14);
    ylabel('$\|u-u_h\|_{L^2}$','Interpreter','latex','FontSize',14);

    p = polyfit(log(h), log(errors(:,2)), 1);
    p_est = p(1);

    title(sprintf('Convergencia $L^2$; orden estimado $\\approx$ %.2f', p_est), ...
          'Interpreter','latex','FontSize',14);

    set(gca, 'XDir','reverse');
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
