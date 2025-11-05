function [S_fine, Sp_fine] = solve_S_eq(r_fine, Sbc, qfun, qpfun, mM, RHS)
% Solve the homogeneous ODE for m=2:
% d/dr (r^3 / q(r)^2 * dS/dr) - r * (m^2-1)/q(r)^2 * S = 0
% with BCs: S(0)=0, S(1)=S2bc
%
% INPUTS:
%   r_fine : vector of points where solution is evaluated (0 excluded)
%   S2bc   : boundary value at r=1
%   qfun   : function handle for q(r)
%
% OUTPUTS:
%   S2_fine  : solution S(r) on r_fine
%   S2p_fine : derivative dS/dr on r_fine
if nargin < 6, RHS = @(r) 0; end

% Define the ODE as a first-order system

odefun = @(r, Y) [Y(2); ...
     -(3 / r - 2 * qpfun(r) / qfun(r)) * Y(2)+(mM^2-1)/r^2 * Y(1)+RHS(r)];

% Boundary conditions
bcfun = @(Ya, Yb) [Ya(1); Yb(1) - Sbc];

% Initial guess: linear profile
solinit = bvpinit(linspace(r_fine(1), r_fine(end), 10), @(r) [Sbc*r; Sbc]);

% Solve with bvp4c
options = bvpset('RelTol',1e-8,'AbsTol',1e-10);
sol = bvp4c(odefun, bcfun, solinit, options);

% Evaluate solution on r_fine
S_fine = deval(sol, r_fine, 1);
Sp_fine = deval(sol, r_fine, 2);
end