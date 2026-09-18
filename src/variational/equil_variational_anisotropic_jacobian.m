function jacobian = equil_variational_anisotropic_jacobian(L,LX,x,cache)
%EQUIL_VARIATIONAL_ANISOTROPIC_JACOBIAN Compatibility wrapper.

    if nargin < 4
        cache = [];
    end
    jacobian = equil_variational_general_jacobian(L,LX,x,cache);
end
