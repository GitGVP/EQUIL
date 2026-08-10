%% Variational and standard counterparts of examples/simple_case.m

beta_final = 30;
epsilon = 0.4;
Ns = 3;
Sbc = zeros(Ns,1);

%% Variational solve: fixed-epsilon beta continuation is internal
[L,LX] = equilVariationalSol('beta',beta_final,'debug',4);
LX.eps_val = epsilon;
LX.Sbc = Sbc;
LY = equilVariationalY(L,LX);

assert(LY.isconverged);
assert(LY.local_B_residual < 1e-10);
fprintf(['Variational: beta=%g, |R|=%.3e, ', ...
    'max|epsilon Delta''|=%.4f\n'],beta_final,norm(LY.residual), ...
    max(abs(epsilon*LY.deltap)));
%% Flux surfaces and constant-omega curves, including an axis zoom
figure;
plot_equilibrium(LY,'variational',false);

figure;
plot(LY.r_plt, LY.S(:,:,1))


function plot_equilibrium(Y,title_text,axis_zoom)
    hold on; axis equal; grid on
    omega_indices = unique(floor(linspace(1,numel(Y.omega_plt),11)));
    if axis_zoom
        radial = Y.r_plt <= 0.2;
        for k = omega_indices
            plot(Y.RR(radial,k),Y.ZZ(radial,k),'--', ...
                'Color',[0,0,0,0.5],'LineWidth',0.7);
        end
        plot(Y.RR(radial,:).',Y.ZZ(radial,:).','Color',[0.2,0.45,0.8,0.35]);
    else
        contour(Y.RR,Y.ZZ,Y.psiN.*ones(size(Y.RR)), ...
            linspace(0,1,11),'LineWidth',1.5);
        for k = omega_indices
            plot(Y.RR(:,k),Y.ZZ(:,k),'--', ...
                'Color',[0,0,0,0.5],'LineWidth',0.5);
        end
    end
    xlabel('$R/R_0$','Interpreter','latex','FontSize',12);
    ylabel('$Z/R_0$','Interpreter','latex','FontSize',12);
    title(title_text,'Interpreter','none');
end
