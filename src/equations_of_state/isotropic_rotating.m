function [dbetapardB0, dbetapardr, dbetapardB, dbetapardR, d2betapardB2, ...
    d2betapardrdB, d2betapardRdB, d2betapardrdR, d2betapardR2, d3betapardrdB2, ...
    d3betapardB3, d3betapardBdR2, d3betapardB2dR, d3betapardrdRdB, betapar, ...
    betaperp] = isotropic_rotating(kinetic_profiles,r,RR,BB)
    beta = kinetic_profiles.beta(r);
    betap = kinetic_profiles.betap(r);
    mach2 = kinetic_profiles.mach2(r);
    mach2p = kinetic_profiles.mach2p(r);
    dbetapardr = exp(mach2.*(-1 + RR.^2)).*(betap + beta.*mach2p.*(-1 + RR.^2));
    dbetapardB = zeros(size(RR));
    dbetapardR = 2.*beta.*exp(mach2.*(-1 + RR.^2)).*mach2.*RR;
    d2betapardB2 = zeros(size(RR));
    d2betapardrdB = zeros(size(RR));
    d2betapardRdB = zeros(size(RR));
    dbetapardB0 = dbetapardB(1,1);
    if nargout > 7 % Jacobian computation
        d2betapardrdR = 2.*exp(mach2.*(-1 + RR.^2)).*RR.*(beta.*mach2p + mach2.*(betap + beta.*mach2p.*(-1 + RR.^2)));
        d2betapardR2 = 2.*beta.*exp(mach2.*(-1 + RR.^2)).*mach2.*(1 + 2.*mach2.*RR.^2);
        d3betapardrdB2 = zeros(size(RR));
        d3betapardB3 = zeros(size(RR));
        d3betapardBdR2 = zeros(size(RR));
        d3betapardB2dR = zeros(size(RR));
        d3betapardrdRdB = zeros(size(RR));
        if nargout > 13 % postprocessing
            betapar   = beta.*exp(mach2.*(-1 + RR.^2));
            betaperp  = beta.*exp(mach2.*(-1 + RR.^2));
        end
    end
end
