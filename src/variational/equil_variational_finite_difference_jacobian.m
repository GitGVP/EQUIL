function jacobian = equil_variational_finite_difference_jacobian( ...
        L,LX,x,cache) %#ok<INUSD>
%EQUIL_VARIATIONAL_FINITE_DIFFERENCE_JACOBIAN Differentiate a residual.
%   This is the generic centered-difference reference used to check an
%   analytical Jacobian.  The selected residual function is responsible
%   for any internal-variable elimination.

    ndof = numel(x);
    jacobian = zeros(ndof,ndof);
    for column = 1:ndof
        step = L.P.jacobian_step*(1+abs(x(column)));
        xp = x;
        xm = x;
        xp(column) = xp(column)+step;
        xm(column) = xm(column)-step;
        Rp = L.P.residuals_fun(L,LX,xp,[]);
        Rm = L.P.residuals_fun(L,LX,xm,[]);
        jacobian(:,column) = (Rp-Rm)/(2*step);
    end
end
