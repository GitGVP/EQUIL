function [Xi_fine, Xip_fine] = solve_Xi_eq(r_fine, Xibc_left, Xibc_right, qfun, qpfun, V, neumann_out)
% Solve the ODE:
%   d/dr [r^3*(2mu-1)^2 * dXi/dr] - 3*r*(2mu-1)^2 * Xi = 0
% where mu(r) = 1/q(r)
%
% INPUTS:
%   r_fine      : vector of r points for evaluation
%   Xibc_left   : BC value at left endpoint (r_fine(1))
%   Xibc_right  : BC value at right endpoint (r_fine(end))
%   qfun        : function handle q(r)
%   qpfun       : function handle dq/dr
%
% OUTPUTS:
%   Xi_fine     : solution xi(r) on r_fine
%   Xip_fine    : derivative dxi/dr on r_fine
if  nargin < 7
    neumann_out = false;
end
if  nargin < 6
    V = @(r) 0; 
    neumann_out = false;
end
% mu and its derivative
mufun  = @(r)  1 ./ qfun(r);
mupfun = @(r) -qpfun(r) ./ qfun(r).^2;
ml=2;
d2coeff = @(r) (-1 + ml./qfun(r)).^2.*r.^2;
d1coeff = @(r) -(((ml - qfun(r)).*r.*(-3.*ml.*qfun(r) + ...
3.*qfun(r).^2 + 2.*ml.*qpfun(r).*r))./qfun(r).^3);
d0coeff = @(r) (1 - ml.^2).*(-1 + ml./qfun(r)).^2;
% Coefficient w(r) = r^3*(2mu-1)^2 and its log-derivative w'/w
% w'/w = 3/r + 2*2*mu'/(2*mu-1)
log_wp = @(r) 3/r + 4*mupfun(r) / (2*mufun(r) - 1);

% ODE as first-order system: Y = [Xi; Xi']
% Xi'' = -(w'/w)*Xi' + 3/r^2 * Xi   [since 3r(2mu-1)^2/w = 3/r^2]
%odefun = @(r, Y) [ Y(2); -log_wp(r) * Y(2) + (3 / r^2) * Y(1) + V(r)];
odefun = @(r, Y) [ Y(2); -(d1coeff(r) * Y(2) + d0coeff(r) * Y(1) + V(r)) ./ d2coeff(r)];
% Boundary conditions: Xi(r_left) = Xibc_left, Xi(r_right) = Xibc_right
if neumann_out
    bcfun = @(Ya, Yb) [ Ya(1) - Xibc_left; Yb(2) - Xibc_right ];
else
    bcfun = @(Ya, Yb) [ Ya(1) - Xibc_left; Yb(1) - Xibc_right ];
end
% Initial guess: linear interpolation between BCs
solinit = bvpinit(linspace(r_fine(1), r_fine(end), 10), @(r) [Xibc_right * r; Xibc_right]);

% Solve
options = bvpset('RelTol', 1e-8, 'AbsTol', 1e-10);
sol = bvp4c(odefun, bcfun, solinit, options);

% Evaluate on fine grid
Xi_fine  = deval(sol, r_fine, 1);
Xip_fine = deval(sol, r_fine, 2);

end