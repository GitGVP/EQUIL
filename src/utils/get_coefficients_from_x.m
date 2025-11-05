function [t2_c, delta_c, P_c, S_c, Bs_c] = get_coefficients_from_x(x, dof_count, Nb, Sbc)
    Ns = numel(Sbc);
    t2_c    = x(1:dof_count);
	delta_c = x(dof_count+1:2*dof_count);
	P_c = x(2*dof_count+1:3*dof_count);
	Bs_c    = reshape(x(3*dof_count+1:(3+Nb)*dof_count), [dof_count, 1, Nb]);
	S_unkn     = reshape(x((3+Nb)*dof_count+1:end), [dof_count-1, 1, Ns]);
	S_c     = cat(1, S_unkn, reshape(Sbc, [1, 1, Ns]));
end