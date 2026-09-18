function basis = assemble_variational_bspline_legacy( ...
        m,nq,p,Ns,Na,vertical_shift,second_derivatives)
%ASSEMBLE_VARIATIONAL_BSPLINE_LEGACY Legacy EQUIL radial spaces.
%   This is the weak-form view of assemble_FE_matrices_bspline_neumann.
%   It uses the same open uniform knots, Gauss points, and endpoint
%   coefficient eliminations, without imposing polar Taylor parity.

    if nargin < 6 || isempty(vertical_shift)
        vertical_shift = false;
    end
    if nargin < 7 || isempty(second_derivatives)
        second_derivatives = false;
    end
    if ~isscalar(m) || m < 1 || m ~= round(m) ...
            || ~isscalar(nq) || nq < 1 || nq ~= round(nq) ...
            || ~isscalar(p) || p < 1 || p ~= round(p)
        error('Use positive integer m, nq, and spline_p.');
    end

    knots = [zeros(1,p+1),(1:m-1)/m,ones(1,p+1)];
    nb = m+p;
    [gauss_points,gauss_weights] = lgwt(nq,0,1);
    gauss_points = flip(gauss_points(:));
    gauss_weights = gauss_weights(:);
    r = zeros(m*nq,1);
    w = zeros(m*nq,1);
    N0 = zeros(nb,m*nq);
    N1 = zeros(nb,m*nq);
    if second_derivatives
        N2 = zeros(nb,m*nq);
    end
    for element = 1:m
        rows = (element-1)*nq+(1:nq);
        left = knots(element+p);
        right = knots(element+p+1);
        element_length = right-left;
        r(rows) = left+element_length*gauss_points;
        w(rows) = element_length*gauss_weights;
        if second_derivatives
            [element_N0,element_N1,element_N2] = ...
                bspline_eval_all(knots,p,r(rows).');
            N2(:,rows) = element_N2;
        else
            [element_N0,element_N1] = ...
                bspline_eval_all(knots,p,r(rows).');
        end
        N0(:,rows) = element_N0;
        N1(:,rows) = element_N1;
    end

    if second_derivatives
        [E0,E1,E2] = bspline_eval_all(knots,p,[0,1]);
    else
        [E0,E1] = bspline_eval_all(knots,p,[0,1]);
    end

    % The first clamped B-spline is removed for every legacy profile.
    % Delta additionally removes the next coefficient.  Fixed-boundary
    % shapes remove the final unknown and use it as a Dirichlet lift.
    full_columns = 2:nb;
    delta_columns = 3:nb;
    shape_columns = 2:nb-1;
    lift_column = nb;

    nshapes = Ns+Na+double(vertical_shift);
    nprofiles = 3+nshapes;
    profile_columns = cell(nprofiles,1);
    profile_columns{1} = full_columns;
    profile_columns{2} = delta_columns;
    profile_columns{3} = full_columns;
    for profile = 4:nprofiles
        profile_columns{profile} = shape_columns;
    end

    basis.B0 = cell(nprofiles,1);
    basis.B1 = cell(nprofiles,1);
    basis.B0_axis = cell(nprofiles,1);
    basis.B1_axis = cell(nprofiles,1);
    basis.B0_edge = cell(nprofiles,1);
    basis.B1_edge = cell(nprofiles,1);
    if second_derivatives
        basis.B2 = cell(nprofiles,1);
        basis.B2_axis = cell(nprofiles,1);
        basis.B2_edge = cell(nprofiles,1);
    end
    basis.profile_lengths = zeros(1,nprofiles);
    for profile = 1:nprofiles
        columns = profile_columns{profile};
        basis.B0{profile} = sparse(N0(columns,:).');
        basis.B1{profile} = sparse(N1(columns,:).');
        basis.B0_axis{profile} = sparse(E0(columns,1).');
        basis.B1_axis{profile} = sparse(E1(columns,1).');
        basis.B0_edge{profile} = sparse(E0(columns,2).');
        basis.B1_edge{profile} = sparse(E1(columns,2).');
        if second_derivatives
            basis.B2{profile} = sparse(N2(columns,:).');
            basis.B2_axis{profile} = sparse(E2(columns,1).');
            basis.B2_edge{profile} = sparse(E2(columns,2).');
        end
        basis.profile_lengths(profile) = numel(columns);
    end

    lift0 = sparse(N0(lift_column,:).');
    lift1 = sparse(N1(lift_column,:).');
    lift0_axis = E0(lift_column,1);
    lift1_axis = E1(lift_column,1);
    lift0_edge = E0(lift_column,2);
    lift1_edge = E1(lift_column,2);
    basis.lift0 = repmat({lift0},nshapes,1);
    basis.lift1 = repmat({lift1},nshapes,1);
    basis.lift0_axis = lift0_axis*ones(nshapes,1);
    basis.lift1_axis = lift1_axis*ones(nshapes,1);
    basis.lift0_edge = lift0_edge*ones(nshapes,1);
    basis.lift1_edge = lift1_edge*ones(nshapes,1);
    if second_derivatives
        lift2 = sparse(N2(lift_column,:).');
        lift2_axis = E2(lift_column,1);
        lift2_edge = E2(lift_column,2);
        basis.lift2 = repmat({lift2},nshapes,1);
        basis.lift2_axis = lift2_axis*ones(nshapes,1);
        basis.lift2_edge = lift2_edge*ones(nshapes,1);
    end

    basis.r = r;
    basis.w = w;
    basis.profile_starts = ...
        cumsum([1,basis.profile_lengths(1:end-1)]);
end
