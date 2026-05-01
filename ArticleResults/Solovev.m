%% Section 4.1: Comparison with rotating Solov'ev equilibria
% This script can be run as is to reproduce Fig. 2 of the article
% Van Parys et al., 
% "Investigation of finite aspect ratio effects in axisymmetric 
% magnetic equilibria with toroidal rotation and pressure anisotropy"

shifts_files = {'shifts_R10_M0_E1p5.mat', 'shifts_R10_M0p5_E1p5.mat', 'shifts_R10_M1_E1p5.mat'};
qp_files     = {'qP_R10_M0_E1p5.mat',     'qP_R10_M0p5_E1p5.mat',     'qP_R10_M1_E1p5.mat'    };

% Parameters of the Solov'ev equilibria
R0     = 10;
mach2s = [0, 0.25, 1];
F0s    = sqrt([450, 570, 930].^2 - [-874.624, -994.913, -1355.8]);
S3bcs  = [16, 12.5, 8] * 1e-3;
betas  = [1.145, 1.10266, 0.825];

colors      = [0.12, 0.47, 0.87; 0.60, 0.20, 0.80; 0.96, 0.30, 0.75];
case_labels = {'$\mathcal M^2=0$', '$\mathcal M^2=0.25$', '$\mathcal M^2=1$'};

figure; hold on;
for k = 1:numel(qp_files)
    shifts = load(shifts_files{k});
    qp     = load(qp_files{k});

    rvals  = qp.Expression1(1,:);
    aminor = rvals(end);
    eps_val = aminor / R0;
    qvals  = qp.Expression1(2,:);

    rhat     = rvals / aminor;
    deltahat = shifts.Expression1(1,:) / (aminor * eps_val);
    q0 = qvals(1);
    s0 = 2 * (qvals(end) - q0) / q0;
    beta = betas(k);

    [L, LX] = equilSol('equation_of_state', @isotropic_rotating, ...
        'mach20', mach2s(k), 'debug', 4, 'q0', q0, 's0', s0, 'beta', beta);
    LX.eps_val = eps_val;
    LX.Sbc(1)  = -0.25 / (eps_val * aminor);
    LX.Sbc(2)  = S3bcs(k) / (eps_val * aminor);
    LY = equilY(L, LX);

    fprintf('bp = %f, wrot = %f, bp+wrot = %f, bp+wrot/2 = %f\n', ...
        LY.bp, LY.Wkrot, LY.bp + LY.Wkrot, LY.bp + 0.5*LY.Wkrot);

    p1 = plot(LY.r_plt, LY.delta, '.', 'Color', colors(k,:));
    idx = round(linspace(1, numel(rhat), 15));
    p2 = plot(rhat(idx), deltahat(idx), 's', 'Color', colors(k,:), ...
        'MarkerSize', 5, 'MarkerFaceColor', 'none', 'LineStyle', 'none');
    p3 = plot(L.r_q, LY.delta_ana, '-.', 'Color', colors(k,:));

    if k == 1
        hCode = p1;  hSolovev = p2;  hExpansion = p3;
    end
end
hold off;

xlabel('$\hat r$', 'Interpreter', 'latex');
ylabel('$\hat \Delta$', 'Interpreter', 'latex');
grid on;

% Line-style legend (upper left)
main_ax = gca;
legend(main_ax, [hCode, hExpansion, hSolovev], {'code', 'expansion', "Solov'ev"}, ...
    'Location', 'NorthWest', 'Box', 'off', 'Interpreter', 'latex');

% Color legend (lower right)
ah = axes('position', main_ax.Position, 'visible', 'off');
hColors = arrayfun(@(k) line(NaN, NaN, 'Color', colors(k,:), 'LineWidth', 2, 'Parent', ah), ...
    1:numel(case_labels));
legend(ah, hColors, case_labels, 'Location', 'SouthEast', 'Box', 'off', 'Interpreter', 'latex');