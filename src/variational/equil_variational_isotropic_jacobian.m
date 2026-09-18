function jacobian = equil_variational_isotropic_jacobian(L,LX,x,cache)
%EQUIL_VARIATIONAL_ISOTROPIC_JACOBIAN Exact reduced isotropic tangent.
%   The native isotropic formulation uses a matrix-form block assembly
%   which averages angular coefficients before radial projection.  This
%   avoids allocating a quadrature-by-angle-by-trial-coefficient array.

    if nargin < 4
        cache = [];
    end
    if ~is_native_isotropic(L,LX)
        jacobian = equil_variational_general_jacobian(L,LX,x,cache);
        return
    end
    jacobian = equil_variational_isotropic_jacobian_matrix(L,LX,x,cache);
end

function tf = is_native_isotropic(L,LX)
    tf = strcmp(func2str(L.P.equation_of_state),'isotropic') ...
        && (~isfield(LX,'pressure_map') || isempty(LX.pressure_map));
end
