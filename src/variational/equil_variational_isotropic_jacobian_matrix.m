function jacobian = equil_variational_isotropic_jacobian_matrix(L,LX,x,cache)
%EQUIL_VARIATIONAL_ISOTROPIC_JACOBIAN_MATRIX Matrix-form isotropic tangent.
%   Each physical variation is represented as
%
%       dot(f) = f_value .* B0 + f_derivative .* B1.
%
%   The angular average is therefore taken before multiplication by the
%   radial trial matrices.  This avoids the Nq-by-Nomega-by-Ncoeff arrays
%   used by the direct vectorized tangent while retaining block assembly.

    if nargin < 4 || isempty(cache) ...
            || ~isfield(cache,'isotropic_reduced')
        [~,cache] = equil_variational_isotropic_residual(L,LX,x,[]);
    end

    state = cache.state;
    epsilon = state.epsilon;
    nprofiles = 3+L.P.Ns+L.P.Na+double(L.P.vertical_shift);
    jacobian = zeros(L.total_dofs,L.total_dofs);

    pressure = cache.fields.pressure;
    coeff_J = cache.fields.coeff_J;
    coeff_goo = cache.fields.coeff_goo;
    edge_coeff_J = isotropic_coeff_J(state.edge, ...
        cache.fields.edge.pressure);

    flux_source = pressure.Pr+state.T.*state.Tr./state.R.^2;
    flux_prefactor = (state.J/epsilon)./state.psir;

    for trial_profile = 1:nprofiles
        columns = block_rows(L,trial_profile);
        direction = profile_direction(L,trial_profile,epsilon,state.omega);
        edge_direction = profile_direction( ...
            L,trial_profile,epsilon,state.edge.omega);
        tangent = kinematic_coefficients(state,direction);
        edge_tangent = kinematic_coefficients(state.edge,edge_direction);

        % Flux-equation row block.
        dflux_derivative = linear_flux_derivative(state,tangent);
        dflux_prefactor.v = tangent.Jv./(epsilon*state.psir) ...
            -state.J.*tangent.psirv./(epsilon*state.psir.^2);
        dflux_prefactor.d = tangent.Jd./(epsilon*state.psir);
        dflux_source.v = state.Tr.*direction.Tv./state.R.^2 ...
            -2*state.T.*state.Tr.*direction.Rv./state.R.^3;
        dflux_source.d = state.T.*direction.Tv./state.R.^2;
        dflux_value.v = flux_source.*dflux_prefactor.v ...
            +flux_prefactor.*dflux_source.v;
        dflux_value.d = flux_source.*dflux_prefactor.d ...
            +flux_prefactor.*dflux_source.d;
        edge_dflux = linear_flux_derivative(state.edge,edge_tangent);
        rows = block_rows(L,1);
        jacobian(rows,columns) = ( ...
            project_linear(L,1,dflux_derivative,trial_profile,true) ...
            -project_linear(L,1,dflux_value,trial_profile,false) ...
            -project_edge_linear(L,1,edge_dflux,trial_profile))/epsilon^2;

        % Gauge-equation row block.
        dgauge.v = tangent.Jv./state.R.^2 ...
            -2*state.J.*direction.Rv./state.R.^3;
        dgauge.d = tangent.Jd./state.R.^2;
        rows = block_rows(L,3);
        jacobian(rows,columns) = ...
            project_linear(L,3,dgauge,trial_profile,false)/epsilon^3;

        % Vacuum-subtracted reduced-action coefficient derivatives.
        [dcoeff_J,dcoeff_R,dcoeff_goo] = action_coefficients( ...
            state,direction,tangent);
        [edge_dcoeff_J,~,~] = action_coefficients( ...
            state.edge,edge_direction,edge_tangent);

        % Delta test block.
        test_J_value = -state.JoverR;
        test_J_derivative = -state.R.*state.Zw;
        dtest_J_value.v = -tangent.Hv;
        dtest_J_value.d = -tangent.Hd;
        dtest_J_derivative.v = -( ...
            state.Zw.*direction.Rv+state.R.*direction.Zwv);
        dtest_J_derivative.d = zeros(size(state.R));
        dvalue = product_rule(test_J_value,dcoeff_J, ...
            coeff_J,dtest_J_value);
        dvalue.v = dvalue.v-dcoeff_R.v;
        dvalue.d = dvalue.d-dcoeff_R.d;
        dderivative = product_rule(test_J_derivative,dcoeff_J, ...
            coeff_J,dtest_J_derivative);

        edge_test_J_derivative = -state.edge.R.*state.edge.Zw;
        edge_dtest_J_derivative.v = -( ...
            state.edge.Zw.*edge_direction.Rv ...
            +state.edge.R.*edge_direction.Zwv);
        edge_dtest_J_derivative.d = zeros(size(state.edge.R));
        edge_dderivative = product_rule(edge_test_J_derivative, ...
            edge_dcoeff_J,edge_coeff_J,edge_dtest_J_derivative);
        rows = block_rows(L,2);
        jacobian(rows,columns) = ( ...
            project_linear(L,2,dvalue,trial_profile,false) ...
            +project_linear(L,2,dderivative,trial_profile,true) ...
            -project_edge_linear(L,2,edge_dderivative,trial_profile)) ...
            /epsilon^4;

        % Symmetric shape test blocks.
        for is = 1:L.P.Ns
            n = is;
            cn = cos(n*state.omega);
            sn = sin(n*state.omega);
            radial_part = -n*state.Rr.*cn+n*state.Zr.*sn;
            test_J_value = state.JoverR.*cn+state.R.*radial_part;
            test_J_derivative = state.R.*(state.Zw.*cn+state.Rw.*sn);
            test_goo_value = -2*n*(state.Rw.*sn+state.Zw.*cn);

            dtest_J_value.v = tangent.Hv.*cn ...
                +direction.Rv.*radial_part;
            dtest_J_value.d = tangent.Hd.*cn+state.R.*( ...
                -n*direction.Rv.*cn+n*direction.Zv.*sn);
            dtest_J_derivative.v = direction.Rv.*( ...
                state.Zw.*cn+state.Rw.*sn)+state.R.*( ...
                direction.Zwv.*cn+direction.Rwv.*sn);
            dtest_J_derivative.d = zeros(size(state.R));
            dtest_goo_value.v = -2*n*( ...
                direction.Rwv.*sn+direction.Zwv.*cn);
            dtest_goo_value.d = zeros(size(state.R));

            dvalue = product_rule(test_J_value,dcoeff_J, ...
                coeff_J,dtest_J_value);
            dvalue.v = dvalue.v+cn.*dcoeff_R.v ...
                +test_goo_value.*dcoeff_goo.v ...
                +coeff_goo.*dtest_goo_value.v;
            dvalue.d = dvalue.d+cn.*dcoeff_R.d ...
                +test_goo_value.*dcoeff_goo.d;
            dderivative = product_rule(test_J_derivative,dcoeff_J, ...
                coeff_J,dtest_J_derivative);
            profile = 3+is;
            rows = block_rows(L,profile);
            jacobian(rows,columns) = ( ...
                project_linear(L,profile,dvalue,trial_profile,false) ...
                +project_linear(L,profile,dderivative,trial_profile,true)) ...
                /epsilon^4;
        end

        % Up-down asymmetric shape test blocks.
        for ia = 1:L.P.Na
            n = L.P.A_modes(ia)-1;
            cn = cos(n*state.omega);
            sn = sin(n*state.omega);
            radial_part = -n*state.Rr.*sn-n*state.Zr.*cn;
            test_J_value = state.JoverR.*sn+state.R.*radial_part;
            test_J_derivative = state.R.*(state.Zw.*sn-state.Rw.*cn);
            test_goo_value = 2*n*(state.Rw.*cn-state.Zw.*sn);

            dtest_J_value.v = tangent.Hv.*sn ...
                +direction.Rv.*radial_part;
            dtest_J_value.d = tangent.Hd.*sn+state.R.*( ...
                -n*direction.Rv.*sn-n*direction.Zv.*cn);
            dtest_J_derivative.v = direction.Rv.*( ...
                state.Zw.*sn-state.Rw.*cn)+state.R.*( ...
                direction.Zwv.*sn-direction.Rwv.*cn);
            dtest_J_derivative.d = zeros(size(state.R));
            dtest_goo_value.v = 2*n*( ...
                direction.Rwv.*cn-direction.Zwv.*sn);
            dtest_goo_value.d = zeros(size(state.R));

            dvalue = product_rule(test_J_value,dcoeff_J, ...
                coeff_J,dtest_J_value);
            dvalue.v = dvalue.v+sn.*dcoeff_R.v ...
                +test_goo_value.*dcoeff_goo.v ...
                +coeff_goo.*dtest_goo_value.v;
            dvalue.d = dvalue.d+sn.*dcoeff_R.d ...
                +test_goo_value.*dcoeff_goo.d;
            dderivative = product_rule(test_J_derivative,dcoeff_J, ...
                coeff_J,dtest_J_derivative);
            profile = 3+L.P.Ns+ia;
            rows = block_rows(L,profile);
            jacobian(rows,columns) = ( ...
                project_linear(L,profile,dvalue,trial_profile,false) ...
                +project_linear(L,profile,dderivative,trial_profile,true)) ...
                /epsilon^4;
        end

        if L.P.vertical_shift
            profile = 4+L.P.Ns+L.P.Na;
            test_J_derivative = -state.R.*state.Rw;
            dtest_J_derivative.v = -( ...
                direction.Rv.*state.Rw+state.R.*direction.Rwv);
            dtest_J_derivative.d = zeros(size(state.R));
            dderivative = product_rule(test_J_derivative,dcoeff_J, ...
                coeff_J,dtest_J_derivative);
            rows = block_rows(L,profile);
            jacobian(rows,columns) = project_linear( ...
                L,profile,dderivative,trial_profile,true)/epsilon^4;
        end
    end

    jacobian = sparse(jacobian);
end

function direction = profile_direction(L,profile,epsilon,omega)
    zero = zeros(size(omega));
    direction.Tv = zero;
    direction.Rv = zero;
    direction.Zv = zero;
    direction.Rwv = zero;
    direction.Zwv = zero;

    if profile == 1
        direction.Tv = epsilon^2+zero;
    elseif profile == 2
        direction.Rv = -epsilon^2+zero;
    elseif profile == 3
        direction.Rv = epsilon^3*cos(omega);
        direction.Zv = epsilon^3*sin(omega);
        direction.Rwv = -epsilon^3*sin(omega);
        direction.Zwv = epsilon^3*cos(omega);
    elseif profile <= 3+L.P.Ns
        n = profile-3;
        direction.Rv = epsilon^2*cos(n*omega);
        direction.Zv = -epsilon^2*sin(n*omega);
        direction.Rwv = -epsilon^2*n*sin(n*omega);
        direction.Zwv = -epsilon^2*n*cos(n*omega);
    elseif profile <= 3+L.P.Ns+L.P.Na
        ia = profile-(3+L.P.Ns);
        n = L.P.A_modes(ia)-1;
        direction.Rv = epsilon^2*sin(n*omega);
        direction.Zv = epsilon^2*cos(n*omega);
        direction.Rwv = epsilon^2*n*cos(n*omega);
        direction.Zwv = -epsilon^2*n*sin(n*omega);
    else
        direction.Zv = epsilon^2+zero;
    end
end

function tangent = kinematic_coefficients(state,direction)
    tangent.Hv = state.Rr.*direction.Zwv ...
        -state.Zr.*direction.Rwv;
    tangent.Hd = state.Zw.*direction.Rv ...
        -state.Rw.*direction.Zv;
    tangent.Jv = state.JoverR.*direction.Rv+state.R.*tangent.Hv;
    tangent.Jd = state.R.*tangent.Hd;
    tangent.goov = 2*(state.Rw.*direction.Rwv ...
        +state.Zw.*direction.Zwv);
    tangent.good = zeros(size(state.R));
    tangent.psirv = state.epsilon*state.r./(state.q*state.a0) ...
        .*direction.Tv+zeros(size(state.R));
    tangent.psird = zeros(size(state.R));

    common = state.psir.^2.*state.goo./state.J.^3;
    tangent.Bp2v = 2*state.epsilon^2*( ...
        state.psir.*state.goo./state.J.^2.*tangent.psirv ...
        +0.5*state.psir.^2./state.J.^2.*tangent.goov ...
        -common.*tangent.Jv);
    tangent.Bp2d = -2*state.epsilon^2*common.*tangent.Jd;
end

function linear = linear_flux_derivative(state,tangent)
    linear.v = state.epsilon*( ...
        state.goo./state.J.*tangent.psirv ...
        +state.psir./state.J.*tangent.goov ...
        -state.psir.*state.goo./state.J.^2.*tangent.Jv);
    linear.d = -state.epsilon*state.psir.*state.goo./state.J.^2 ...
        .*tangent.Jd;
end

function [dcoeff_J,dcoeff_R,dcoeff_goo] = action_coefficients( ...
        state,direction,tangent)
    shifted_T2 = state.T.^2-state.a0^2;
    dcoeff_J.v = -0.5*tangent.Bp2v ...
        -state.T./state.R.^2.*direction.Tv ...
        +shifted_T2./state.R.^3.*direction.Rv;
    dcoeff_J.d = -0.5*tangent.Bp2d;
    dcoeff_R.v = shifted_T2./state.R.^3.*tangent.Jv ...
        +2*state.J.*state.T./state.R.^3.*direction.Tv ...
        -3*state.J.*shifted_T2./state.R.^4.*direction.Rv;
    dcoeff_R.d = shifted_T2./state.R.^3.*tangent.Jd;
    dcoeff_goo.v = state.epsilon^2*( ...
        state.psir./state.J.*tangent.psirv ...
        -0.5*state.psir.^2./state.J.^2.*tangent.Jv);
    dcoeff_goo.d = -0.5*state.epsilon^2*state.psir.^2./state.J.^2 ...
        .*tangent.Jd;
end

function result = product_rule(a,db,b,da)
    result.v = a.*db.v+b.*da.v;
    result.d = a.*db.d+b.*da.d;
end

function block = project_linear( ...
        L,test_profile,linear,trial_profile,use_test_derivative)
    average_value = mean(linear.v,2);
    average_derivative = mean(linear.d,2);
    weighted_trials = spdiags(L.w_r.*average_value,0, ...
            length(L.w_r),length(L.w_r))*L.B0{trial_profile} ...
        +spdiags(L.w_r.*average_derivative,0, ...
            length(L.w_r),length(L.w_r))*L.B1{trial_profile};
    if use_test_derivative
        block = L.B1{test_profile}'*weighted_trials;
    else
        block = L.B0{test_profile}'*weighted_trials;
    end
end

function block = project_edge_linear(L,test_profile,linear,trial_profile)
    average_value = mean(linear.v,2);
    average_derivative = mean(linear.d,2);
    trial = average_value*L.B0_edge{trial_profile} ...
        +average_derivative*L.B1_edge{trial_profile};
    block = L.B0_edge{test_profile}'*trial;
end

function coeff_J = isotropic_coeff_J(state,pressure)
    shifted_T2 = state.T.^2-state.a0^2;
    coeff_J = -0.5*state.Bp2-pressure.Pi ...
        -0.5*shifted_T2./state.R.^2;
end

function rows = block_rows(L,profile)
    first = L.profile_starts(profile);
    rows = first:first+L.profile_lengths(profile)-1;
end
