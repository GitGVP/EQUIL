function [residual, fields] = equil_variational_assemble(L, LX, state, B)
%EQUIL_VARIATIONAL_ASSEMBLE Direct weak-form residual assembly.
%   Only profile values and first radial derivatives enter this routine.

    fixed = equil_variational_B_constraint(L, LX, state, B);
    edge = equil_variational_local_B(L, LX, state.edge, []);
    pressure = fixed.pressure;
    sigma = fixed.sigma;
    Bphi = fixed.Bphi;
    omega = state.omega;
    w = L.w_r;

    residual = zeros(L.total_dofs, 1);

    % Weak t2 equation: the poloidal term carries T through psir, while
    % the explicit T*Tr term comes from -T*Bphi/R in the action.
    Js = state.J/state.epsilon;
    flux_derivative_kernel = ...
        (1-sigma).*state.goo./Js.*state.psir;
    flux_value_kernel = ...
        Js.*pressure.Pr./state.psir ...
        + Js.*state.T.*state.Tr./ ...
          (state.R.^2.*(1-sigma).*state.psir);
    row = block_rows(L, 1);
    residual(row) = L.B1{1}'*(w.*mean(flux_derivative_kernel, 2)) ...
                  - L.B0{1}'*(w.*mean(flux_value_kernel, 2));
    edge_flux = (1-edge.sigma).*state.edge.goo./ ...
        (state.edge.J/state.epsilon).*state.edge.psir;
    residual(row) = residual(row) ...
        - L.B0_edge{1}'*mean(edge_flux, 2);
    residual(row) = residual(row)/state.epsilon^2;

    % Compact first variation of Phi=J*L at fixed independent B.
    action_density = 0.5*B.^2-pressure.Pi-state.T.*Bphi./state.R;
    Tshift = state.T-state.a0;
    coeff_J = 0.5*state.Bp2-pressure.Pi ...
        -state.T.*state.Bp2./(state.R.*Bphi) ...
        +0.5*(Bphi-state.T./state.R).^2 ...
        -0.5*Tshift.*(state.T+state.a0)./state.R.^2;
    toroidal_R_kernel = state.T.*(state.R.*Bphi-state.T) ...
        +Tshift.*(state.T+state.a0);
    coeff_R = state.J.*(-pressure.PR+toroidal_R_kernel./state.R.^3);
    coeff_goo = state.T.*(state.epsilon*state.psir).^2./ ...
                (2*state.R.*Bphi.*state.J);

    % The constant 1/R toroidal vacuum field has identically zero bulk
    % force.  Remove its action variation analytically before quadrature;
    % otherwise O(epsilon) volume/edge terms cancel to O(epsilon^4), which
    % loses accuracy in aspect-ratio scans near epsilon=1e-3.

    % Delta variation: delta R=-v, delta Z=0.  The common epsilon^2
    % variation factor is divided from this equation as a row scaling.
    delta_value = -coeff_J.*state.JoverR-coeff_R;
    delta_derivative = -coeff_J.*state.R.*state.Zw;
    row = block_rows(L, 2);
    residual(row) = L.B0{2}'*(w.*mean(delta_value, 2)) ...
                  + L.B1{2}'*(w.*mean(delta_derivative, 2));
    edge_Tshift = state.edge.T-state.a0;
    edge_coeff_J = 0.5*state.edge.Bp2-edge.pressure.Pi ...
        -state.edge.T.*state.edge.Bp2./(state.edge.R.*edge.Bphi) ...
        +0.5*(edge.Bphi-state.edge.T./state.edge.R).^2 ...
        -0.5*edge_Tshift.*(state.edge.T+state.a0)./state.edge.R.^2;
    edge_delta_derivative = ...
        -edge_coeff_J.*state.edge.R.*state.edge.Zw;
    residual(row) = residual(row) ...
        - L.B0_edge{2}'*mean(edge_delta_derivative, 2);
    residual(row) = residual(row)/state.epsilon^4;

    % S_m variations assembled from delta J and delta g_omegaomega.
    % the Euler-Lagrange equations or taking radial second derivatives.
    for is = 1:L.P.Ns
        n = is;
        cn = cos(n*omega);
        sn = sin(n*omega);
        deltaJ_value = state.JoverR.*cn + state.R.*( ...
            -n*state.Rr.*cn+n*state.Zr.*sn);
        deltaJ_derivative = state.R.*(state.Zw.*cn+state.Rw.*sn);
        delta_goo_value = -2*n*(state.Rw.*sn+state.Zw.*cn);

        shape_value = coeff_J.*deltaJ_value+coeff_R.*cn ...
                    + coeff_goo.*delta_goo_value;
        shape_derivative = coeff_J.*deltaJ_derivative;
        profile = 3+is;
        row = block_rows(L, profile);
        residual(row) = L.B0{profile}'*(w.*mean(shape_value, 2)) ...
                      + L.B1{profile}'*(w.*mean(shape_derivative, 2));
        residual(row) = residual(row)/state.epsilon^4;
    end

    % Up-down asymmetric A_m variations, n=m-1: delta R=sin(n*omega)*v
    % and delta Z=+cos(n*omega)*v after dividing out epsilon^2.
    for ia = 1:L.P.Na
        n = L.P.A_modes(ia)-1;
        cn = cos(n*omega);
        sn = sin(n*omega);
        deltaJ_value = state.JoverR.*sn + state.R.*( ...
            -n*state.Rr.*sn-n*state.Zr.*cn);
        deltaJ_derivative = state.R.*(state.Zw.*sn-state.Rw.*cn);
        delta_goo_value = 2*n*(state.Rw.*cn-state.Zw.*sn);

        shape_value = coeff_J.*deltaJ_value+coeff_R.*sn ...
                    +coeff_goo.*delta_goo_value;
        shape_derivative = coeff_J.*deltaJ_derivative;
        profile = 3+L.P.Ns+ia;
        row = block_rows(L,profile);
        residual(row) = L.B0{profile}'*(w.*mean(shape_value,2)) ...
                      +L.B1{profile}'*(w.*mean(shape_derivative,2));
        residual(row) = residual(row)/state.epsilon^4;
    end

    % Optional vertical-center variation: delta R=0, delta Z=v(r).
    % It supplies the missing mean-Z boundary degree of freedom without
    % entering (or changing) the symmetric formulation when disabled.
    if L.P.vertical_shift
        deltaJ_derivative = -state.R.*state.Rw;
        shape_derivative = coeff_J.*deltaJ_derivative;
        profile = 4+L.P.Ns+L.P.Na;
        row = block_rows(L,profile);
        residual(row) = L.B1{profile}'*( ...
            w.*mean(shape_derivative,2));
        residual(row) = residual(row)/state.epsilon^4;
    end

    % Exact anisotropic SFL periodicity constraint.  Division by epsilon^2
    % follows the normalized paper convention and improves row scaling.
    gauge = mean(state.J./(state.R.^2.*(1-sigma)), 2) ...
            /state.epsilon^2-state.r/state.a0;
    row = block_rows(L, 3);
    residual(row) = L.B0{3}'*(w.*gauge)/state.epsilon;

    fields = fixed;
    fields.edge = edge;
    fields.action_density = action_density;
    fields.gauge = gauge;
    fields.flux_derivative_kernel = flux_derivative_kernel;
    fields.flux_value_kernel = flux_value_kernel;
end

function rows = block_rows(L, profile)
    first = L.profile_starts(profile);
    rows = first:first+L.profile_lengths(profile)-1;
end
