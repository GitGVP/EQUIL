function [dbetapardB0, dbetapardr, dbetapardB, dbetapardR, d2betapardB2, ...
    d2betapardrdB, d2betapardRdB, d2betapardrdR, d2betapardR2, d3betapardrdB2, ...
    d3betapardB3, d3betapardBdR2, d3betapardB2dR, d3betapardrdRdB, ...
    d3betapardrdR2, betapar, betaperp] = de_Blank(kinetic_profiles,r,RR,BB)
    beta = kinetic_profiles.beta(r);
    betap = kinetic_profiles.betap(r);
    mach2 = kinetic_profiles.mach2(r);
    mach2p = kinetic_profiles.mach2p(r);
    Theta = kinetic_profiles.Theta(r);
    Thetap = kinetic_profiles.Thetap(r);
    dbetapardr = -((BB.*exp(mach2.*(-1 + RR.^2)).*(betap.*(BB - Theta).*(-1 + Theta) + beta.*(mach2p.*(-1 + RR.^2).*(BB - Theta).*(-1 + Theta) + (-1 + BB).*Thetap)))./(BB - Theta).^2);
    dbetapardB = (beta.*exp(mach2.*(-1 + RR.^2)).*(-1 + Theta).*Theta)./(BB - Theta).^2;
    dbetapardR = (-2.*BB.*beta.*exp(mach2.*(-1 + RR.^2)).*mach2.*RR.*(-1 + Theta))./(BB - Theta);
    d2betapardB2 = (2.*beta.*exp(mach2.*(-1 + RR.^2)).*(-1 + Theta).*Theta)./(-BB + Theta).^3;
    d2betapardrdB = (exp(mach2.*(-1 + RR.^2)).*((1 + BB).*(betap + beta.*mach2p.*(-1 + RR.^2)).*Theta.^2 - (betap + beta.*mach2p.*(-1 + RR.^2)).*Theta.^3 - BB.*beta.*Thetap - Theta.*(BB.*betap + beta.*(BB.*mach2p.*(-1 + RR.^2) + (1 - 2.*BB).*Thetap))))./(BB - Theta).^3;
    d2betapardRdB = (2.*beta.*exp(mach2.*(-1 + RR.^2)).*mach2.*RR.*(-1 + Theta).*Theta)./(BB - Theta).^2;
    dbetapardB0 = dbetapardB(1,1);
    if nargout > 7 % Jacobian computation
        d2betapardrdR = (-2.*BB.*exp(mach2.*(-1 + RR.^2)).*RR.*(beta.*mach2p.*(BB - Theta).*(-1 + Theta) + mach2.*(betap.*(BB - Theta).*(-1 + Theta) + beta.*(mach2p.*(-1 + RR.^2).*(BB - Theta).*(-1 + Theta) + (-1 + BB).*Thetap))))./(BB - Theta).^2;
        d2betapardR2 = (-2.*BB.*beta.*exp(mach2.*(-1 + RR.^2)).*mach2.*(1 + 2.*mach2.*RR.^2).*(-1 + Theta))./(BB - Theta);
        d3betapardrdB2 = (-2.*exp(mach2.*(-1 + RR.^2)).*(-((betap + beta.*mach2p.*(-1 + RR.^2)).*Theta.^3) - BB.*beta.*Thetap + Theta.^2.*((1 + BB).*betap + beta.*((1 + BB).*mach2p.*(-1 + RR.^2) + Thetap)) - Theta.*(BB.*betap + beta.*(BB.*mach2p.*(-1 + RR.^2) - 2.*(-1 + BB).*Thetap))))./(BB - Theta).^4;
        d3betapardB3 = (6.*beta.*exp(mach2.*(-1 + RR.^2)).*(-1 + Theta).*Theta)./(BB - Theta).^4;
        d3betapardBdR2 = (2.*beta.*exp(mach2.*(-1 + RR.^2)).*mach2.*(1 + 2.*mach2.*RR.^2).*(-1 + Theta).*Theta)./(BB - Theta).^2;
        d3betapardB2dR = (4.*beta.*exp(mach2.*(-1 + RR.^2)).*mach2.*RR.*(-1 + Theta).*Theta)./(-BB + Theta).^3;
        d3betapardrdRdB = (2.*exp(mach2.*(-1 + RR.^2)).*RR.*(-((1 + BB).*(beta.*mach2p + mach2.*(betap + beta.*mach2p.*(-1 + RR.^2))).*Theta.^2) + (beta.*mach2p + mach2.*(betap + beta.*mach2p.*(-1 + RR.^2))).*Theta.^3 + BB.*beta.*mach2.*Thetap + Theta.*(BB.*beta.*mach2p + mach2.*(BB.*betap + beta.*(BB.*mach2p.*(-1 + RR.^2) + (1 - 2.*BB).*Thetap)))))./(-BB + Theta).^3;
        d3betapardrdR2 = (-2.*BB.*exp(mach2.*(-1 + RR.^2)).*(beta.*mach2p.*(BB - Theta).*(-1 + Theta) + 2.*mach2.^2.*RR.^2.*(betap.*(BB - Theta).*(-1 + Theta) + beta.*(mach2p.*(-1 + RR.^2).*(BB - Theta).*(-1 + Theta) + (-1 + BB).*Thetap)) + mach2.*(betap.*(BB - Theta).*(-1 + Theta) + beta.*(mach2p.*(-1 + 5.*RR.^2).*(BB - Theta).*(-1 + Theta) + (-1 + BB).*Thetap))))./(BB - Theta).^2;
        if nargout > 14 % postprocessing
            betapar   = -((BB.*beta.*exp(mach2.*(-1 + RR.^2)).*(-1 + Theta))./(BB - Theta));
            betaperp  = -((BB.^2.*beta.*exp(mach2.*(-1 + RR.^2)).*(-1 + Theta))./(BB - Theta).^2);
        end
    end
end
