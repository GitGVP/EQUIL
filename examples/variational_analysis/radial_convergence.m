%% h-refinement of the axis-regular B-spline discretization

if ~exist('use_pedestal','var')
    use_pedestal = false;
end
m_values = [3,4,5,6,8,10];
reference_m = 14;
spline_p = 4;
beta0 = 0.2;
epsilon = 0.2;

if use_pedestal
    center = 0.65;
    width = 0.10;
    scale = atan(center/width)-atan((center-1)/width);
    beta_shape = @(r)(atan((center-r.^2)/width) ...
        -atan((center-1)/width))/scale;
    betap_shape = @(r)(-2*r/width) ...
        ./(1+((center-r.^2)/width).^2)/scale;
end

all_m = [m_values,reference_m];
solutions = cell(size(all_m));
previous_L = [];
previous_Y = [];
for k = 1:numel(all_m)
    [L,X] = equilVariationalSol('beta',beta0,'Ns',1, ...
        'm',all_m(k),'nq',6,'spline_p',spline_p,'om_pts',64, ...
        'NLtol',1e-11,'nk',35,'debug',0);
    X.eps_val = epsilon;
    X.Sbc = 0;
    if use_pedestal
        X.beta_shape = beta_shape;
        X.betap_shape = betap_shape;
        X.kinetic_profiles.beta = @(r)beta0*beta_shape(r);
        X.kinetic_profiles.betap = @(r)beta0*betap_shape(r);
    end
    if ~isempty(previous_Y)
        X.x = transfer_state(previous_L,L,previous_Y.x);
    end
    Y = equilVariationalY(L,X);
    assert(Y.isconverged);
    solutions{k} = Y;
    previous_L = L;
    previous_Y = Y;
end

r = linspace(0,1,2001).';
reference = solutions{end};
reference_profiles = { ...
    interp1(reference.r_plt,reference.t2,r,'pchip'), ...
    interp1(reference.r_plt,reference.t2p,r,'pchip'), ...
    interp1(reference.r_plt,reference.delta,r,'pchip'), ...
    interp1(reference.r_plt,reference.deltap,r,'pchip')};
errors = zeros(numel(m_values),4);
for k = 1:numel(m_values)
    Y = solutions{k};
    profiles = {Y.t2,Y.t2p,Y.delta,Y.deltap};
    for field = 1:4
        value = interp1(Y.r_plt,profiles{field},r,'pchip');
        errors(k,field) = sqrt(trapz(r,(value-reference_profiles{field}).^2));
    end
end

h = 1./m_values;
orders = zeros(1,4);
for field = 1:4
    fit = polyfit(log(h),log(errors(:,field).'),1);
    orders(field) = fit(1);
end
fprintf(['Expected smooth-solution orders: L2 about %d, H1 about %d.\n', ...
    'Measured t2: L2 %.2f, H1 %.2f.\n', ...
    'Delta is pre-asymptotic on this mesh range: fitted slopes %.2f and %.2f.\n'], ...
    spline_p+1,spline_p,orders);
disp(table(m_values(:),h(:),errors(:,1),errors(:,2), ...
    errors(:,3),errors(:,4),'VariableNames', ...
    {'m','h','t2_L2','t2_H1','Delta_L2','Delta_H1'}));

figure;
tiledlayout(2,2,'TileSpacing','compact','Padding','compact');
titles = {'$\|t_2-t_{2,h}\|_{L^2}$', ...
          '$|t_2-t_{2,h}|_{H^1}$', ...
          '$\|\Delta-\Delta_h\|_{L^2}$', ...
          '$|\Delta-\Delta_h|_{H^1}$'};
for field = 1:4
    nexttile;
    loglog(h,errors(:,field),'ok-','LineWidth',1.5);
    grid on
    set(gca,'XDir','reverse');
    xlabel('$h=1/m$','Interpreter','latex');
    title(titles{field},'Interpreter','latex');
    legend(sprintf('slope %.2f',orders(field)), ...
        'Location','southeast','Box','off');
end

figure; hold on; grid on
show = [1,3,5,numel(m_values)];
colors = turbo(numel(show));
for j = 1:numel(show)
    k = show(j);
    plot(solutions{k}.r_plt,solutions{k}.delta,'LineWidth',1.5, ...
        'Color',colors(j,:),'DisplayName',sprintf('m=%d',m_values(k)));
end
plot(reference.r_plt,reference.delta,'k--','LineWidth',2, ...
    'DisplayName',sprintf('reference m=%d',reference_m));
xlabel('$\hat r$','Interpreter','latex');
ylabel('$\hat\Delta$','Interpreter','latex');
legend('Location','northwest','Box','off');


function xnew = transfer_state(Lold,Lnew,xold)
    xnew = zeros(Lnew.total_dofs,1);
    nprofiles = 3+Lnew.P.Ns+Lnew.P.Nh;
    for profile = 1:nprofiles
        old_rows = Lold.profile_starts(profile) ...
            +(0:Lold.profile_lengths(profile)-1);
        new_rows = Lnew.profile_starts(profile) ...
            +(0:Lnew.profile_lengths(profile)-1);
        old_values = Lold.B0{profile}*xold(old_rows);
        new_values = interp1(Lold.r_q,old_values,Lnew.r_q,'pchip');
        xnew(new_rows) = Lnew.B0{profile}\new_values;
    end
end
