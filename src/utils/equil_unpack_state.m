function [t2, t2p, delta, deltap, deltapp, P, Pp, Ppp, Bs, Bsp, S, Sp, Spp] = ...
          equil_unpack_state(x, dof_count, Nb, Sbc, P0, P1, P2)

Ns = numel(Sbc);
idx = 1;

t2_c = x(idx:idx+dof_count-1);                idx = idx + dof_count;
delta_c = x(idx:idx+dof_count-2);             idx = idx + (dof_count-1);
delta_cf = [0; delta_c];
P_c = x(idx:idx+dof_count-1);                 idx = idx + dof_count;
Bs_c = reshape(x(idx:idx+Nb*dof_count-1), [dof_count, 1, Nb]);  idx = idx + Nb*dof_count;
S_unkn = reshape(x(idx:end), [dof_count-1, 1, Ns]);

S_c = cat(1, S_unkn, reshape(Sbc, [1, 1, Ns]));
Bs_flat = reshape(Bs_c, size(Bs_c,1), []);
S_flat  = reshape(S_c,  size(S_c,1),  []);

t2      = P0 * t2_c;
t2p     = P1 * t2_c;

delta   = P0 * delta_cf;
deltap  = P1 * delta_cf;
deltapp = P2 * delta_cf;

P       = P0 * P_c;
Pp      = P1 * P_c;
Ppp     = P2 * P_c;

Bs      = reshape(P0 * Bs_flat, size(P0,1), 1, []);
Bsp     = reshape(P1 * Bs_flat, size(P1,1), 1, []);

S       = reshape(P0 * S_flat, size(P0,1), 1, []);
Sp      = reshape(P1 * S_flat, size(P1,1), 1, []);
Spp     = reshape(P2 * S_flat, size(P2,1), 1, []);

end

% -------------------------------------------------------------------------
% Old version (no Neumann at axis):
%
% function [t2_c, delta_c, delta_cf, P_c, Bs_c, Bs_flat, S_unkn, S_c, S_flat, ...
%           t2, t2p, delta, deltap, deltapp, P, Pp, Ppp, Bs, Bsp, S, Sp, Spp] = ...
%           equil_unpack_state(x, dof_count, Nb, Sbc, P0, P1, P2)
%
% Ns = numel(Sbc);
% idx = 1;
%
% t2_c = x(idx:idx+dof_count-1);                idx = idx + dof_count;
% delta_c = x(idx:idx+dof_count-1);             idx = idx + dof_count;
% delta_cf = delta_c;
% P_c = x(idx:idx+dof_count-1);                 idx = idx + dof_count;
% Bs_c = reshape(x(idx:idx+Nb*dof_count-1), [dof_count, 1, Nb]);  idx = idx + Nb*dof_count;
% S_unkn = reshape(x(idx:end), [dof_count-1, 1, Ns]);
%
% S_c = cat(1, S_unkn, reshape(Sbc, [1, 1, Ns]));
% Bs_flat = reshape(Bs_c, size(Bs_c,1), []);
% S_flat  = reshape(S_c,  size(S_c,1),  []);
%
% t2      = P0 * t2_c;
% t2p     = P1 * t2_c;
% delta   = P0 * delta_cf;
% deltap  = P1 * delta_cf;
% deltapp = P2 * delta_cf;
% P       = P0 * P_c;
% Pp      = P1 * P_c;
% Ppp     = P2 * P_c;
% Bs      = reshape(P0 * Bs_flat, size(P0,1), 1, []);
% Bsp     = reshape(P1 * Bs_flat, size(P1,1), 1, []);
% S       = reshape(P0 * S_flat, size(P0,1), 1, []);
% Sp      = reshape(P1 * S_flat, size(P1,1), 1, []);
% Spp     = reshape(P2 * S_flat, size(P2,1), 1, []);
%
% end