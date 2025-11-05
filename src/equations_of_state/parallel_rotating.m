function [dbetapardB0, dbetapardr, dbetapardB, dbetapardR, d2betapardB2, ...
    d2betapardrdB, d2betapardRdB, d2betapardrdR, d2betapardR2, d3betapardrdB2, ...
    d3betapardB3, d3betapardBdR2, d3betapardB2dR, d3betapardrdRdB] = parallel_rotating(kinetic_profiles,r,RR,BB)
    beta = kinetic_profiles.beta(r);
    betap = kinetic_profiles.betap(r);
    mach2 = kinetic_profiles.mach2(r);
    mach2p = kinetic_profiles.mach2p(r);
    dbetapardr = BB.*betap.*exp(mach2.*(-1 + RR.^2)) + BB.*beta.*exp(mach2.*(-1 + RR.^2)).*mach2p.*(-1 + RR.^2);
    dbetapardB = beta.*exp(mach2.*(-1 + RR.^2));
    dbetapardR = 2.*BB.*beta.*exp(mach2.*(-1 + RR.^2)).*mach2.*RR;
    d2betapardB2 = zeros(size(RR));
    d2betapardrdB = betap.*exp(mach2.*(-1 + RR.^2)) + beta.*exp(mach2.*(-1 + RR.^2)).*mach2p.*(-1 + RR.^2);
    d2betapardRdB = 2.*beta.*exp(mach2.*(-1 + RR.^2)).*mach2.*RR;
    dbetapardB0 = dbetapardB(1,1);
    if nargout > 7
        d2betapardrdR = 2.*BB.*betap.*exp(mach2.*(-1 + RR.^2)).*mach2.*RR + 2.*BB.*beta.*exp(mach2.*(-1 + RR.^2)).*mach2p.*RR + 2.*BB.*beta.*exp(mach2.*(-1 + RR.^2)).*mach2.*mach2p.*RR.*(-1 + RR.^2);
        d2betapardR2 = 2.*BB.*beta.*exp(mach2.*(-1 + RR.^2)).*mach2 + 4.*BB.*beta.*exp(mach2.*(-1 + RR.^2)).*mach2.^2.*RR.^2;
        d3betapardrdB2 = zeros(size(RR));
        d3betapardB3 = zeros(size(RR));
        d3betapardBdR2 = 2.*beta.*exp(mach2.*(-1 + RR.^2)).*mach2 + 4.*beta.*exp(mach2.*(-1 + RR.^2)).*mach2.^2.*RR.^2;
        d3betapardB2dR = zeros(size(RR));
        d3betapardrdRdB = 2.*betap.*exp(mach2.*(-1 + RR.^2)).*mach2.*RR + 2.*beta.*exp(mach2.*(-1 + RR.^2)).*mach2p.*RR + 2.*beta.*exp(mach2.*(-1 + RR.^2)).*mach2.*mach2p.*RR.*(-1 + RR.^2);
    end
end
