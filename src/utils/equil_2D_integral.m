function I = equil_2D_integral(integrand, r, theta, Jac, include_axis)
% INT2D  Area integral on a polar-like grid using midpoint rule.
%   I = int2D(integrand, r, theta, Jac, include_axis)
%
%   integrand : Nr x Ntheta
%   r         : Nr x 1
%   theta     : 1  x Ntheta
%   Jac       : Nr x Ntheta   (e.g. J or area jacobian, J/R)
%   include_axis : if false, skip the first radial point (r=0 singularity)

if ~include_axis
    integrand = integrand(2:end, :);
    Jac       = Jac(2:end, :);
    r         = r(2:end);
end

dr     = diff(r);
dtheta = diff(theta);

dA  = Jac(1:end-1, 1:end-1) .* dr(:) .* dtheta(:).';
I   = sum(integrand(1:end-1, 1:end-1)  .* dA, 'all');

end