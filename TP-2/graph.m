% ------------------------ run_errors_simple.m ------------------------
clear; clc;

Q.a0 = -1; Q.b0 = 1;
Q.x  = [-0.9739065285171717; -0.8650633666889845; -0.6794095682990244; ...
        -0.4333953941292472; -0.1488743389816312;  0.1488743389816312; ...
         0.4333953941292472;  0.6794095682990244;  0.8650633666889845; ...
         0.9739065285171717];
Q.w  = [0.06667134430868814; 0.1494513491505806; 0.2190863625159820; ...
        0.2692667193099964;  0.2955242247147529; 0.2955242247147529; ...
        0.2692667193099964;  0.2190863625159820; 0.1494513491505806; ...
        0.06667134430868814];

alphas = [0.5, 1.0, 1.5, 2.0, 2.5];
K  = 10;
hs = 2.^(-(1:K));
E  = zeros(K, numel(alphas));


for k = 1:K
    h  = 2^(-k);
    N  = round(2 / h);
    Th.x = linspace(-1, 1, N+1).';

    A = projectionMatrix(Th);

    for j = 1:numel(alphas)
        alpha = alphas(j);
        f = @(x) abs(x).^alpha;

        vec = projectionVector(Th, f, Q);
        c   = A \ vec;

        E(k, j) = projectionError(Th, f, c, Q); % ||f - Pi_1(f)||_varphi
    end
end

for j = 1:numel(alphas)
    fprintf('alpha = %.1f: ', alphas(j));
    fprintf('%.6e ', E(:, j).');
    fprintf('\n');
end
fprintf('\n');


mask = abs(alphas - 1.0) > 1e-12;
figure; hold on; grid on; box on;
for j = find(mask)
    plot(hs, E(:, j), 'o-', 'LineWidth', 1.25, 'MarkerSize', 5);
end
set(gca, 'XScale', 'log', 'YScale', 'log');   % fuerza ejes log

xlabel('$h$','Interpreter','latex','FontSize',16);
ylabel('$\|\,f - \Pi_{1}(f)\,\|_{\varphi}$','Interpreter','latex','FontSize',16);
leg = arrayfun(@(a) sprintf('\\alpha=%.1f', a), alphas(mask), 'UniformOutput', false);
legend(leg, 'Location', 'south east');
