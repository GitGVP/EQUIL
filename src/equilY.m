function LY = equilY(L, LX)
  % Analytical profiles if requested
  if L.P.do_ana
     LY = equil_ana(L, LX);
  end

  % Default initial guess if not provided in LX
  if ~isfield(LX, 'x') || isempty(LX.x)
     LX.x = equil_x0(L, LX);
  end

  % Initial state
  x = LX.x;
  if L.P.debug > 10
    [t2, ~, delta, ~, ~, ~, ~, ~, ~, ~, S, ~, ~] = ...
          equil_unpack_state(x, L.dof_count, L.P.Nb, LX.Sbc, L.P0, L.P1, L.P2);  
    figure;
    tiledlayout(3,1,'TileSpacing','compact','Padding','compact');
    ax1 = nexttile(1);
    hold on;grid on;
    ylabel('$t_2$', 'Interpreter','latex', 'FontSize',14)
    ax2 = nexttile(2);
    hold on;grid on;
    ylabel('$\hat \Delta$', 'Interpreter','latex', 'FontSize',14)
    ax3 = nexttile(3);
    hold on;grid on;
    xlabel('$\hat r$', 'Interpreter','latex', 'FontSize',14)
    ylabel('$\hat S_2$', 'Interpreter','latex', 'FontSize',14)
    plot(ax1, L.r_q, t2,'k--','LineWidth',3,...
        'DisplayName', sprintf('$k=%i$',0))
    plot(ax2, L.r_q, delta,'k--','LineWidth',3,...
        'DisplayName', sprintf('$k=%i$',0))
    plot(ax3, L.r_q, S(:,:,1),'k--','LineWidth',3,...
        'DisplayName', sprintf('$k=%i$',0))
    drawnow;
  end
  res_norms = [];LY.isconverged = false;
  % Newton iterations
  for k = 1:L.P.nk
    %Compute residuals and Jacobian at quadrature points
    residuals = L.P.residuals_fun(L.r_q,x,LX.eps_val,L.omega,...
        LX.kinetic_profiles,LX.qfun(L.r_q),LX.qpfun(L.r_q), L.P0,L.P1,L.P2,L.M_profiles,L.dof_count,L.P.Nb,LX.Sbc,L.P.equation_of_state);
    J = L.P.jacobian_fun(L.r_q,x,LX.eps_val,L.omega,...
        LX.kinetic_profiles,LX.qfun(L.r_q),LX.qpfun(L.r_q), L.P0, L.P1, L.P2,L.dof_count,L.P.Nb,LX.Sbc,L.P_templates, L.M_extended,L.P.equation_of_state);
    update = J \ residuals;
    % try with damping
    x = x - L.P.damping*update;
    res_norm = norm(residuals);
    res_norms(k) = res_norm;
    if L.P.debug > 1
      fprintf('Iter %d, eps %.1e, |res| = %.4e, |Δx| = %.3e\n', ...
              k, LX.eps_val, res_norm, norm(update));
      if L.P.debug > 10
          [t2, ~, delta, ~, ~, ~, ~, ~, ~, ~, S, ~, ~] = ...
          equil_unpack_state(x, L.dof_count, L.P.Nb, LX.Sbc, L.P0, L.P1, L.P2);  
            plot(ax1, L.r_q, t2,'.', 'DisplayName', sprintf('$k=%i$',k))
            plot(ax2, L.r_q, delta,'.', 'DisplayName', sprintf('$k=%i$',k))
            plot(ax3, L.r_q, S(:,:,1),'.','DisplayName', sprintf('$k=%i$',k))
            drawnow;
      end
    end
    % Break if converged
    if res_norm < L.P.NLtol
        LY.isconverged = true;
      break;
    end
  end
  if L.P.debug > 10
    legend('Interpreter','latex','Box','off');
  end
  LY.x = x;
  LY.res_norms = res_norms;
  LY = equilPP(L,LX,LY);
end