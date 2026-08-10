function jacobian = equil_variational_jacobian(L, LX, x, cache)
%EQUIL_VARIATIONAL_JACOBIAN Consistent locally condensed Jacobian.
%   At fixed B, centered directional differences give G_x.  The internal
%   field response is then formed explicitly as B_x=-G_x/G_B.  Centered
%   differences of the residual along (x,B(x)) produce the reduced
%   Jacobian R_x-R_B*G_B^{-1}*G_x without global B degrees of freedom.

    if nargin < 4 || isempty(cache)
        [~, cache] = equil_variational_residual(L, LX, x, []);
    end
    B = cache.local.B;
    GB = cache.local.GB;
    if any(GB(:) <= L.P.min_GB)
        error('Cannot condense the Jacobian at or beyond the mirror boundary.');
    end

    ndof = numel(x);
    jacobian = zeros(ndof, ndof);
    for column = 1:ndof
        base_step = L.P.jacobian_step*(1+abs(x(column)));
        differentiated = false;
        for attempt = 1:12
            step = base_step/2^(attempt-1);
            xp = x;
            xm = x;
            xp(column) = xp(column)+step;
            xm(column) = xm(column)-step;
            try
                state_p = equil_variational_state(L,LX,xp);
                state_m = equil_variational_state(L,LX,xm);
                Gp = equil_variational_B_constraint(L,LX,state_p,B);
                Gm = equil_variational_B_constraint(L,LX,state_m,B);
                Gx = (Gp.G-Gm.G)/(2*step);
                Bx = -Gx./GB;
                Bplus = B+step*Bx;
                Bminus = B-step*Bx;
                if any(Bplus(:).^2 <= state_p.Bp2(:)) ...
                        || any(Bminus(:).^2 <= state_m.Bp2(:))
                    continue
                end
                Rp = equil_variational_assemble(L,LX,state_p,Bplus);
                Rm = equil_variational_assemble(L,LX,state_m,Bminus);
                jacobian(:,column) = (Rp-Rm)/(2*step);
                differentiated = true;
                break
            catch
                % Reduce the directional step near geometric/B constraints.
            end
        end
        if ~differentiated
            error('Could not differentiate condensed column %d.',column);
        end
    end
end
