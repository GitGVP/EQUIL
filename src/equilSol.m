function [L, LX, LY] = equilSol(varargin)
% Top level function, inspired by LIU from MEQ.
    
    setup_EQUIL_paths();
    % defaults and parameters
    P = equilP(varargin{:});
    % matrices assembly, other utilities
    L = equilL(P);
    % code inputs (BCs, profiles, aspect ratio)
    if nargout > 1
        LX = equilX(L);
        if nargout > 2
            LY = equilY(L,LX);
        end
    end
end

function P = equilP(varargin)
 % default parameters
 P.m = 15; 
 P.nq = 10; % adapt to p for B-spline
 P.spline_p = P.nq; % might be a bit risky to set like this
 P.axis_regular_basis = false; % profile-specific radial powers/parity
 P.axis_regular_max_power = []; % [] -> m+spline_p
 P.om_pts = 300;               % quadrature points in omega
 P.equation_of_state = @isotropic;
 P.residuals_fun = @residuals_noRepl;
 P.jacobian_fun = @jacobian_noRepl;
 P.Ns = 3; % (R,Z) cosine harmonics
 P.Nb = 1; % B cosine harmonics 
 P.nk = 30; % Newton iterations
 % input default parameters
 P.debug = 0;
 P.damping = 1;
 P.NLtol = 5e-14;
 
 P.do_SFL = false;
 P.do_ana = true;
 P.do_shift_NLO = false;
 P.q0 = 1;
 P.s0 = 4;
 
 % kinetic profiles' defaults
 P.beta = 0.3;
 P.A0 = 1;
 P.Bc0 = 1+0.3; 
 P.mach20 = 0;
 P.gamma0 = 1;
 P.Theta0 = 0;
 
 % overwrites
 for k = 1:2:length(varargin)
  P.(varargin{k}) = varargin{k+1};
 end
end

function L = equilL(P)
  L.r_nodes = linspace(0, 1, 2*P.m+1);   % node positions
  L.omega = (0:P.om_pts-1) * 2*pi / P.om_pts;        % nodes: 0, 2pi/N, ..., 2pi*(N-1)/N
  
  
  if P.axis_regular_basis
      current_residuals = ismember(func2str(P.residuals_fun), ...
          {'residuals_noRepl', 'residuals_iso_static'});
      current_jacobians = ismember(func2str(P.jacobian_fun), ...
          {'jacobian_noRepl', 'jacobian_iso_static'});
      if ~current_residuals || ~current_jacobians
          error(['axis_regular_basis is implemented only for the current ', ...
                 'noRepl and iso_static residual/Jacobian paths.']);
      end
      [L.r_q, L.P0, L.P1, L.P2, L.P3, ...
       L.M_profiles, L.M_extended, L.P_templates, L.P_extended, ...
       L.A_global, L.profile_lengths, L.profile_starts, L.P0_end] = ...
          assemble_FE_matrices_polynomial_axis_regular( ...
              L.r_nodes, P.Nb, P.Ns, P.nq, P.spline_p, ...
              P.axis_regular_max_power);
  else
      [L.r_q, L.P0, L.P1, L.P2, L.P3, ...
       L.M_profiles, L.M_extended, L.P_templates, L.P_extended, ...
       L.A_global, L.profile_lengths, L.profile_starts, L.P0_end] = ...
          assemble_FE_matrices_bspline_neumann( ...
              L.r_nodes, P.Nb, P.Ns, P.nq, P.spline_p);
  end
  % numeric check (coerce if necessary)
  if iscell(L.profile_lengths), L.profile_lengths = cell2mat(L.profile_lengths); end

  L.dof_count = L.profile_lengths(1);   % first profile block size
  L.Nq        = numel(L.r_q);          % number of quadrature points per profile
  
  [lu_L, lu_U, lu_P, lu_Q] = lu(L.A_global);
  L.lu = struct('L', lu_L, 'U', lu_U, 'P', lu_P, 'Q', lu_Q);
  L.P = P;
end


function setup_EQUIL_paths()
% Add project subfolders if they exist and are not already on the path.
% Silent on subsequent calls unless new folders are added.

    root = locate_project_root();
    if isempty(root)
        root = pwd;
    end

    want = { fullfile(root,'src'), ...
             fullfile(root,'src','utils'), ...
             fullfile(root,'src','equations_of_state'), ...
             fullfile(root,'src','residuals_jacobians'), ...
             fullfile(root,'examples') };

    curpath = strsplit(path, pathsep);
    added = false;

    for k = 1:numel(want)
        p = want{k};
        if ~exist(p,'dir')
            continue
        end
        % check if already on path (exact directory match)
        if ~any(strcmp(curpath, p))
            addpath(p);
            added = true;
        end
    end

    if added
        fprintf('Found project root at %s\nAdded paths for EQUIL\n', root);
    end
end


function root = locate_project_root()
    % Try to locate a reasonable project root relative to this file.
    % Heuristics: look upwards for 'src' dir, README.md or .git folder.
    mf = mfilename('fullpath');
    if isempty(mf)
        root = '';
        return
    end
    cur = fileparts(mf);
    maxdepth = 6;
    root = '';
    for i = 1:maxdepth
        if exist(fullfile(cur,'src'),'dir') || exist(fullfile(cur,'README.md'),'file') || exist(fullfile(cur,'.git'),'dir')
            root = cur;
            return
        end
        parent = fileparts(cur);
        if isempty(parent) || strcmp(parent,cur)
            break
        end
        cur = parent;
    end
end
