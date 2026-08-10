function [L, LX, LY] = equilVariationalSol(varargin)
%EQUILVARIATIONALSOL Variational EQUIL proof-of-principle solver.
%   The global unknowns are [t2, Delta, P, S_2, ..., S_(Ns+1),
%   V_2, ..., V_(Nh+1)]. Magnetic field strength is eliminated locally.

    setup_variational_paths();
    P = variational_parameters(varargin{:});
    L = variational_discretization(P);

    if nargout > 1
        LX = equilVariationalX(L);
        if nargout > 2
            LY = equilVariationalY(L, LX);
        end
    end
end

function P = variational_parameters(varargin)
    P.m = 8;
    P.nq = 6;
    P.spline_p = 4;
    P.om_pts = 96;
    P.Ns = 3;
    P.Nh = 0;

    P.equation_of_state = @isotropic;
    P.beta = 0.3;
    P.q0 = 1;
    P.s0 = 4;
    P.A0 = 1;
    P.Bc0 = 1.3;
    P.mach20 = 0;
    P.Theta0 = 0;
    P.gamma0 = 1;
    P.do_ana = false;
    P.do_shift_NLO = false;

    P.nk = 30;
    P.NLtol = 1e-10;
    P.NLstepTol = 1e-11;
    P.damping = 1;
    P.beta_continuation = true;
    P.jacobian_step = eps^(1/3);
    P.local_B_tol = 5e-13;
    P.local_B_maxit = 20;
    P.min_one_minus_sigma = 1e-8;
    P.min_GB = 1e-10;
    P.debug = 0;

    if mod(numel(varargin), 2) ~= 0
        error('Options must be supplied as name/value pairs.');
    end
    for k = 1:2:numel(varargin)
        name = varargin{k};
        if ~isfield(P,name) && ~strcmp(name,'Mach20')
            error('Unknown variational option %s.',name);
        end
        P.(name) = varargin{k+1};
    end
    if isfield(P, 'Mach20')
        P.mach20 = P.Mach20;
    end
end

function L = variational_discretization(P)
    L.r_nodes = linspace(0, 1, P.m + 1);
    L.omega = (0:P.om_pts-1) * (2*pi/P.om_pts);
    basis = assemble_variational_bspline_axis_regular( ...
        P.m,P.nq,P.spline_p,P.Ns,P.Nh);
    L.r_q = basis.r;
    L.w_r = basis.w;
    L.B0 = basis.B0;
    L.B1 = basis.B1;
    L.B0_axis = basis.B0_axis;
    L.B1_axis = basis.B1_axis;
    L.B0_edge = basis.B0_edge;
    L.B1_edge = basis.B1_edge;
    L.profile_lengths = basis.profile_lengths;
    L.profile_starts = basis.profile_starts;
    L.Sbc0 = basis.lift0(1:P.Ns);
    L.Sbc1 = basis.lift1(1:P.Ns);
    L.Sbc0_axis = basis.lift0_axis(1:P.Ns);
    L.Sbc1_axis = basis.lift1_axis(1:P.Ns);
    L.Sbc0_edge = basis.lift0_edge(1:P.Ns);
    L.Sbc1_edge = basis.lift1_edge(1:P.Ns);
    shape = P.Ns+(1:P.Nh);
    L.Vbc0 = basis.lift0(shape);
    L.Vbc1 = basis.lift1(shape);
    L.Vbc0_axis = basis.lift0_axis(shape);
    L.Vbc1_axis = basis.lift1_axis(shape);
    L.Vbc0_edge = basis.lift0_edge(shape);
    L.Vbc1_edge = basis.lift1_edge(shape);
    L.Nq = numel(L.r_q);
    L.dof_count = L.profile_lengths(1);
    L.total_dofs = sum(L.profile_lengths);
    L.P = P;
end

function setup_variational_paths()
    here = fileparts(mfilename('fullpath'));
    wanted = {here, fullfile(here, 'utils'), ...
              fullfile(here, 'equations_of_state'), ...
              fullfile(here, 'variational')};
    current = strsplit(path, pathsep);
    for k = 1:numel(wanted)
        if exist(wanted{k}, 'dir') && ~any(strcmp(current, wanted{k}))
            addpath(wanted{k});
        end
    end
end
