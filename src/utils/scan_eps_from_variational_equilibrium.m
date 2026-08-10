function scan = scan_eps_from_variational_equilibrium(L,LX0,LY0,varargin)
%SCAN_EPS_FROM_VARIATIONAL_EQUILIBRIUM Fixed-physics aspect-ratio scan.
%   Each smaller-epsilon equilibrium starts from the preceding solution.
%   Epsilon is never used as a continuation path for the initial solve.

    if isempty(varargin)
        nsim = 10;
    else
        nsim = varargin{1};
    end
    values = flip(logspace(log10(1e-3),log10(LX0.eps_val),nsim));
    nr = numel(L.r_q);
    scan.epsilon = values;
    scan.residual = zeros(1,nsim);
    scan.t2 = zeros(nr,nsim);
    scan.delta = zeros(nr,nsim);
    scan.alpha = zeros(nr,nsim);
    scan.S2 = zeros(nr,nsim);
    scan.solutions = cell(1,nsim);

    LY = LY0;
    for k = 1:nsim
        if k > 1
            LX = LX0;
            LX.eps_val = values(k);
            LX.x = LY.x;
            LY = equilVariationalY(L,LX);
        end
        if ~LY.isconverged
            error('Variational epsilon scan failed at epsilon %.4g.',values(k));
        end
        scan.residual(k) = norm(LY.residual);
        scan.t2(:,k) = LY.t2(2:end-1);
        scan.delta(:,k) = LY.delta(2:end-1);
        scan.alpha(:,k) = values(k)*LY.deltap(2:end-1);
        scan.S2(:,k) = LY.S(2:end-1,1,1);
        scan.solutions{k} = LY;
        fprintf('Variational epsilon scan %d/%d: eps=%.3e, |R|=%.3e\n', ...
                k,nsim,values(k),scan.residual(k));
    end

    colors = turbo(nsim);
    figure;
    loglog(values,scan.residual,'o-');
    grid on
    xlabel('$\epsilon$','Interpreter','latex');
    ylabel('$\|R\|_2$','Interpreter','latex');

    figure;
    tiledlayout(2,2,'TileSpacing','compact','Padding','compact');
    names = {'$t_2$','$\Delta$','$\epsilon\Delta''$','$S_2$'};
    profiles = {scan.t2,scan.delta,scan.alpha,scan.S2};
    for panel = 1:4
        nexttile; hold on; grid on
        for k = 1:nsim
            plot(L.r_q,profiles{panel}(:,k),'Color',colors(k,:));
        end
        xlabel('$r/a$','Interpreter','latex');
        ylabel(names{panel},'Interpreter','latex');
    end
end
