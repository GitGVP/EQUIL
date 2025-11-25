function [dbetapardB0, dbetapardr, dbetapardB, dbetapardR, d2betapardB2, ...
    d2betapardrdB, d2betapardRdB, d2betapardrdR, d2betapardR2, d3betapardrdB2, ...
    d3betapardB3, d3betapardBdR2, d3betapardB2dR, d3betapardrdRdB, ...
    d3betapardrdR2, betapar, betaperp] = parallel_rotating(kinetic_profiles,r,RR,BB)
    beta = kinetic_profiles.beta(r);
    betap = kinetic_profiles.betap(r);
    mach2 = kinetic_profiles.mach2(r);
    mach2p = kinetic_profiles.mach2p(r);
    dbetapardr = BB.*exp(mach2.*(-1 + RR.^2)).*(betap + beta.*mach2p.*(-1 + RR.^2));
    dbetapardB = beta.*exp(mach2.*(-1 + RR.^2));
    dbetapardR = 2.*BB.*beta.*exp(mach2.*(-1 + RR.^2)).*mach2.*RR;
    d2betapardB2 = zeros(size(RR));
    d2betapardrdB = exp(mach2.*(-1 + RR.^2)).*(betap + beta.*mach2p.*(-1 + RR.^2));
    d2betapardRdB = 2.*beta.*exp(mach2.*(-1 + RR.^2)).*mach2.*RR;
    dbetapardB0 = dbetapardB(1,1);
    if nargout > 7 % Jacobian computation
        d2betapardrdR = 2.*BB.*exp(mach2.*(-1 + RR.^2)).*RR.*(beta.*mach2p + mach2.*(betap + beta.*mach2p.*(-1 + RR.^2)));
        d2betapardR2 = 2.*BB.*beta.*exp(mach2.*(-1 + RR.^2)).*mach2.*(1 + 2.*mach2.*RR.^2);
        d3betapardrdB2 = zeros(size(RR));
        d3betapardB3 = zeros(size(RR));
        d3betapardBdR2 = 2.*beta.*exp(mach2.*(-1 + RR.^2)).*mach2.*(1 + 2.*mach2.*RR.^2);
        d3betapardB2dR = zeros(size(RR));
        d3betapardrdRdB = 2.*exp(mach2.*(-1 + RR.^2)).*RR.*(beta.*mach2p + mach2.*(betap + beta.*mach2p.*(-1 + RR.^2)));
        d3betapardrdR2 = 2.*BB.*exp(mach2.*(-1 + RR.^2)).*(beta.*mach2p + 2.*mach2.^2.*RR.^2.*(betap + beta.*mach2p.*(-1 + RR.^2)) + mach2.*(betap + beta.*mach2p.*(-1 + 5.*RR.^2)));
        if nargout > 14 % postprocessing
            betapar   = BB.*beta.*exp(mach2.*(-1 + RR.^2));
            betaperp  = 0;
        end
    end
end
