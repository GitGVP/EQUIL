%% High-shift equilibria at matched edge displacement
% The three cases share epsilon, q(r), zero shaping BCs and resolution.
% Their beta values differ because the de Blank firehose and mirror bounds
% prevent strongly anisotropic states at a common beta.

epsilon = 0.4;
Ns = 3;
common_args = {'Ns',Ns,'q0',1,'m',4,'nq',12,'om_pts',64, ...
    'NLtol',1e-10,'nk',70,'debug',0};

% High-beta isotropic reference.
[Liso,Xiso] = equilVariationalSol(common_args{:}, ...
    'equation_of_state',@isotropic,'beta',7.24);
Xiso.eps_val = epsilon;
Xiso.Sbc = zeros(Ns,1);
Yiso = equilVariationalY(Liso,Xiso);
assert(Yiso.isconverged);

% Strong parallel pressure: p_perp/p_parallel=1/6 on axis. The adaptive
% beta controller reaches the target directly from its beta=0 seed.
parallel_beta_stages = 2.91;
[Lparallel,Xparallel,Yparallel] = continue_beta(common_args,epsilon, ...
    -5,parallel_beta_stages);

% Perpendicular branch: solve isotropic beta=20, then jump directly to the
% requested anisotropy while retaining the condensed local-B state.
perpendicular_theta_stages = [0,0.093];
[Lperp,Xperp,Yperp] = continue_theta(common_args,epsilon, ...
    20,perpendicular_theta_stages);

labels = {'isotropic','strong parallel','near-mirror perpendicular'};
betas = [7.24,parallel_beta_stages(end),20].';
theta0 = [0,-5,perpendicular_theta_stages(end)].';
solutions = {Yiso,Yparallel,Yperp};
edge_shift = cellfun(@(Y)Y.delta(end),solutions).';
shift_shear = cellfun(@(Y)max(abs(epsilon*Y.deltap)),solutions).';
min_scaled_J = cellfun(@(Y)Y.min_J_over_eps2r,solutions).';
max_elongation = cellfun(@(Y)max(Y.kappa),solutions).';
firehose_margin = cellfun(@(Y)Y.firehose_margin,solutions).';
mirror_margin = cellfun(@(Y)Y.mirror_margin,solutions).';
axis_pressure_ratio = cellfun(@(Y) ...
    Y.Pi_perp(1,1)/Y.Pi_parallel(1,1),solutions).';
residual = cellfun(@(Y)norm(Y.residual),solutions).';
dominant_block = strings(3,1);
for k = 1:3
    [~,iblock] = max(solutions{k}.residual_block_norms);
    dominant_block(k) = solutions{k}.residual_block_names{iblock};
end
summary = table(labels(:),betas,theta0,edge_shift,shift_shear, ...
    min_scaled_J,max_elongation,firehose_margin,mirror_margin, ...
    axis_pressure_ratio,residual,dominant_block,'VariableNames', ...
    {'case_name','beta','Theta0','Delta_edge','max_epsilon_Deltap', ...
     'min_J_over_epsilon2r','max_kappa','firehose_margin', ...
     'mirror_margin','axis_pperp_over_ppar','residual','dominant_block'});
disp(summary);
assert(max(edge_shift)-min(edge_shift) < 0.015, ...
    'The selected edge shifts are no longer matched.');
assert(max(max_elongation)-min(max_elongation) > 0.5, ...
    'The selected cases no longer resolve a strong elongation contrast.');

colors = [0.10,0.10,0.10; 0.20,0.45,0.85; 0.90,0.25,0.35];
Rlimits = [min(cellfun(@(Y)min(Y.RR,[],'all'),solutions)), ...
           max(cellfun(@(Y)max(Y.RR,[],'all'),solutions))];
Zmax = max(cellfun(@(Y)max(abs(Y.ZZ),[],'all'),solutions));

%% Flux-surface comparison on common axes
figure;
tiledlayout(1,3,'TileSpacing','compact','Padding','compact');
for k = 1:3
    nexttile; hold on; axis equal; grid on
    Y = solutions{k};
    contour(Y.RR,Y.ZZ,Y.psiN.*ones(size(Y.RR)), ...
        linspace(0,1,11),'LineWidth',1.5);
    contour(Y.RR,Y.ZZ,(Y.betaperp-min(Y.betaperp(:)))/(max(Y.betaperp(:))-min(Y.betaperp(:))), ...
        linspace(0,1,11),'r--','LineWidth',1.5);
    plot(Y.RR(end,:), Y.ZZ(end,:), 'k')
    xlim(Rlimits); ylim([-Zmax,Zmax]);
    title(labels{k},'Interpreter','none');
    xlabel('$R/R_0$','Interpreter','latex');
    if k == 1
        ylabel('$Z/R_0$','Interpreter','latex');
    end
end


figure;hold on; axis equal; grid on
for k = 1:3
    Y = solutions{k};
    plot(Y.RR(end,:), Y.ZZ(end,:))
end


figure;hold on;grid on
for k = 1:3
    Y = solutions{k};
    plot(Y.r_plt, Y.P)
end

%% Shift, elongation harmonic and geometrical elongation
figure;
tiledlayout(1,3,'TileSpacing','compact','Padding','compact');
nexttile; hold on; grid on
for k = 1:3
    plot(solutions{k}.r_plt,solutions{k}.delta,'LineWidth',2, ...
        'Color',colors(k,:),'DisplayName',labels{k});
end
xlabel('$\hat r$','Interpreter','latex');
ylabel('$\hat\Delta$','Interpreter','latex');
legend('Location','northwest','Box','off');

nexttile; hold on; grid on
for k = 1:3
    plot(solutions{k}.r_plt,squeeze(solutions{k}.S(:,:,1)), ...
        'LineWidth',2,'Color',colors(k,:),'DisplayName',labels{k});
end
xlabel('$\hat r$','Interpreter','latex');
ylabel('$\hat S_2$','Interpreter','latex');

nexttile; hold on; grid on
for k = 1:3
    plot(solutions{k}.r_plt,solutions{k}.kappa,'LineWidth',2, ...
        'Color',colors(k,:),'DisplayName',labels{k});
end
xlabel('$\hat r$','Interpreter','latex');
ylabel('$\kappa$','Interpreter','latex');

%% Radial packing and stability margins
figure;
tiledlayout(1,3,'TileSpacing','compact','Padding','compact');
nexttile; hold on; grid on
for k = 1:3
    Y = solutions{k};
    scaled_J = Y.J(2:end,:)./(epsilon^2*Y.r_plt(2:end));
    plot(Y.r_plt(2:end),min(scaled_J,[],2),'LineWidth',2, ...
        'Color',colors(k,:),'DisplayName',labels{k});
end
xlabel('$\hat r$','Interpreter','latex');
ylabel('$\min_\omega J/(\epsilon^2\hat r)$','Interpreter','latex');
legend('Location','best','Box','off');

nexttile; hold on; grid on
for k = 1:3
    plot(solutions{k}.r_plt,min(1-solutions{k}.sigma,[],2), ...
        'LineWidth',2,'Color',colors(k,:),'DisplayName',labels{k});
end
yline(0,'k:');
xlabel('$\hat r$','Interpreter','latex');
ylabel('$\min_\omega(1-\sigma)$','Interpreter','latex');

nexttile; hold on; grid on
for k = 1:3
    plot(solutions{k}.r_plt,min(solutions{k}.GB,[],2), ...
        'LineWidth',2,'Color',colors(k,:),'DisplayName',labels{k});
end
yline(0,'k:');
xlabel('$\hat r$','Interpreter','latex');
ylabel('$\min_\omega G_B$','Interpreter','latex');


function [L,X,Y] = continue_beta(common_args,epsilon,theta,beta_stages)
    Y = [];
    for beta = beta_stages
        [L,X] = equilVariationalSol(common_args{:}, ...
            'equation_of_state',@de_Blank,'Theta0',theta,'beta',beta);
        X.eps_val = epsilon;
        X.Sbc = zeros(L.P.Ns,1);
        if ~isempty(Y)
            X.x = Y.x;
            X.local_B_guess = Y.local_B_quadrature;
        end
        Y = equilVariationalY(L,X);
        assert(Y.isconverged,'Parallel beta continuation failed at %g.',beta);
    end
end

function [L,X,Y] = continue_theta(common_args,epsilon,beta,theta_stages)
    Y = [];
    for theta = theta_stages
        [L,X] = equilVariationalSol(common_args{:}, ...
            'equation_of_state',@de_Blank,'Theta0',theta,'beta',beta);
        X.eps_val = epsilon;
        X.Sbc = zeros(L.P.Ns,1);
        if ~isempty(Y)
            X.x = Y.x;
            X.local_B_guess = Y.local_B_quadrature;
        end
        Y = equilVariationalY(L,X);
        assert(Y.isconverged, ...
            'Perpendicular Theta continuation failed at %g.',theta);
    end
end
