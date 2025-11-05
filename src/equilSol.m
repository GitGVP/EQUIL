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
 P.hot_restart = false;
 P.m = 15; 
 P.nq = 10; % adapt to p for B-spline
 P.spline_p = P.nq; % might be a bit risky to set like this
 P.om_pts = 300;               % quadrature points in omega
 P.equation_of_state = @anisotropic_minB;
 P.residuals_fun = @residuals_noEL;
 P.jacobian_fun = @jacobian_anyDF;
 P.Ns = 3; % (R,Z) cosine harmonics
 P.Nb = 1; % B cosine harmonics 
 P.nk = 30; % Newton iterations
 % input default parameters
 P.debug = 0;

 P.beta = 0;
 P.q0 = 1;
 P.s0 = 0.5;
 P.A0 = 1;
 
 % overwrites
 for k = 1:2:length(varargin)
  P.(varargin{k}) = varargin{k+1};
 end
end

function L = equilL(P)
  L.r_nodes = linspace(0, 1, 2*P.m+1);   % node positions
  L.omega = (0:P.om_pts-1) * 2*pi / P.om_pts;        % nodes: 0, 2pi/N, ..., 2pi*(N-1)/N
  
  
  [L.r_q, L.P0, L.P1, L.P2, L.M, ...
  L.P0_S, L.P1_S, L.P2_S, L.M_S, ...
  L.M_profiles, L.M_extended, L.P_templates, L.P_extended, L.M_template, ...
  L.A_global, L.profile_lengths, L.profile_starts, L.P0_end] = ...
  assemble_FE_matrices_bspline(L.r_nodes, P.Nb, P.Ns, P.nq, P.spline_p);
  % numeric check (coerce if necessary)
  if iscell(L.profile_lengths), L.profile_lengths = cell2mat(L.profile_lengths); end

  L.dof_count = L.profile_lengths(1);   % first profile block size (2*m normally)
  L.Nq        = numel(L.r_q);          % number of quadrature points per profile
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