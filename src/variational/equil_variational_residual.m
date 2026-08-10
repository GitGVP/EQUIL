function [residual, cache] = equil_variational_residual(L, LX, x, Bguess)
%EQUIL_VARIATIONAL_RESIDUAL Eliminate local B and assemble global residual.

    if nargin < 4
        Bguess = [];
    end
    state = equil_variational_state(L, LX, x);
    local = equil_variational_local_B(L, LX, state, Bguess);
    [residual, fields] = ...
        equil_variational_assemble(L, LX, state, local.B);

    cache.state = state;
    cache.local = local;
    cache.fields = fields;
    cache.residual = residual;
end
