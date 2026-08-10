function constraint = equil_variational_B_constraint(L, LX, state, B)
%EQUIL_VARIATIONAL_B_CONSTRAINT G(B)=0 and its local B derivative.

    if any(B(:).^2 <= state.Bp2(:))
        error('A fixed-B evaluation has B^2 <= Bp2.');
    end
    pressure = equil_variational_pressure(L, LX, state.r, state.R, B);
    Bphi = sqrt(B.^2-state.Bp2);
    sigma = pressure.PB./B;
    sigma_B = pressure.PBB./B-pressure.PB./B.^2;

    constraint.G = 1-sigma-state.T./(state.R.*Bphi);
    constraint.GB = -sigma_B+state.T.*B./(state.R.*Bphi.^3);
    constraint.Bphi = Bphi;
    constraint.sigma = sigma;
    constraint.sigma_B = sigma_B;
    constraint.pressure = pressure;
end
