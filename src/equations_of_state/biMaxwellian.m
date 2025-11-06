function [dbetapardB0, dbetapardr, dbetapardB, dbetapardR, d2betapardB2, ...
    d2betapardrdB, d2betapardRdB, d2betapardrdR, d2betapardR2, d3betapardrdB2, ...
    d3betapardB3, d3betapardBdR2, d3betapardB2dR, d3betapardrdRdB, betapar, ...
    betaperp] = biMaxwellian(kinetic_profiles,r,RR,BB)
    beta = kinetic_profiles.beta(r) * ones(1,size(RR,2));
    betap = kinetic_profiles.betap(r) * ones(1,size(RR,2));
    Ah = kinetic_profiles.Ah(r) * ones(1,size(RR,2));
    Ahp = kinetic_profiles.Ahp(r) * ones(1,size(RR,2));
    Bc = kinetic_profiles.Bc(r) * ones(1,size(RR,2));
    Bcp = kinetic_profiles.Bcp(r) * ones(1,size(RR,2));
    idx = BB < Bc;
    dbetapardr = (BB.*((Ahp.*(-BB + Bc) + (-1 + Ah).*Bcp).*beta + (Ah.*(BB - Bc) + Bc).*betap))./(Ah.*(BB - Bc) + Bc).^2;
    dbetapardr(idx) = (BB(idx).*((Ahp(idx).*(-BB(idx) + Bc(idx)) + (-1 + Ah(idx)).*Bcp(idx)).*beta(idx) + (Ah(idx).*(BB(idx) - Bc(idx)) + Bc(idx)).*betap(idx)))./(Ah(idx).*(BB(idx) - Bc(idx)) + Bc(idx)).^2;
    dbetapardB = -(((-1 + Ah).*Bc.*beta)./(Ah.*(BB - Bc) + Bc).^2);
    dbetapardB(idx) = -(((-1 + Ah(idx)).*Bc(idx).*beta(idx))./(Ah(idx).*(BB(idx) - Bc(idx)) + Bc(idx)).^2);
    dbetapardR = zeros(size(RR));
    d2betapardB2 = (2.*(-1 + Ah).*Ah.*Bc.*beta)./(Ah.*(BB - Bc) + Bc).^3;
    d2betapardB2(idx) = (2.*(-1 + Ah(idx)).*Ah(idx).*Bc(idx).*beta(idx))./(Ah(idx).*(BB(idx) - Bc(idx)) + Bc(idx)).^3;
    d2betapardrdB = (-((-1 + Ah).*Ah.*BB.*Bcp.*beta) - (-1 + Ah).*Bc.^2.*(Ahp.*beta - (-1 + Ah).*betap) + Bc.*(((-2 + Ah).*Ahp.*BB - (-1 + Ah).^2.*Bcp).*beta - (-1 + Ah).*Ah.*BB.*betap))./(Ah.*(BB - Bc) + Bc).^3;
    d2betapardrdB(idx) = (-((-1 + Ah(idx)).*Ah(idx).*BB(idx).*Bcp(idx).*beta(idx)) - (-1 + Ah(idx)).*Bc(idx).^2.*(Ahp(idx).*beta(idx) - (-1 + Ah(idx)).*betap(idx)) + Bc(idx).*(((-2 + Ah(idx)).*Ahp(idx).*BB(idx) - (-1 + Ah(idx)).^2.*Bcp(idx)).*beta(idx) - (-1 + Ah(idx)).*Ah(idx).*BB(idx).*betap(idx)))./(Ah(idx).*(BB(idx) - Bc(idx)) + Bc(idx)).^3;
    d2betapardRdB = zeros(size(RR));
    dbetapardB0 = dbetapardB(1,1);
    if nargout > 7 % Jacobian computation
        d2betapardrdR = zeros(size(RR));
        d2betapardR2 = zeros(size(RR));
        d3betapardrdB2 = (2.*(-(Ahp.*Bc.^2.*beta) + Ah.*Bc.*(2.*(Ahp.*BB + Bcp).*beta - Bc.*betap) + Ah.^3.*((BB + 2.*Bc).*Bcp.*beta + (BB - Bc).*Bc.*betap) - Ah.^2.*(BB.*Bcp.*beta - Bc.^2.*(Ahp.*beta + 2.*betap) + Bc.*(Ahp.*BB.*beta + 4.*Bcp.*beta + BB.*betap))))./(Ah.*(BB - Bc) + Bc).^4;
        d3betapardrdB2(idx) = (2.*(-(Ahp(idx).*Bc(idx).^2.*beta(idx)) + Ah(idx).*Bc(idx).*(2.*(Ahp(idx).*BB(idx) + Bcp(idx)).*beta(idx) - Bc(idx).*betap(idx)) + Ah(idx).^3.*((BB(idx) + 2.*Bc(idx)).*Bcp(idx).*beta(idx) + (BB(idx) - Bc(idx)).*Bc(idx).*betap(idx)) - Ah(idx).^2.*(BB(idx).*Bcp(idx).*beta(idx) - Bc(idx).^2.*(Ahp(idx).*beta(idx) + 2.*betap(idx)) + Bc(idx).*(Ahp(idx).*BB(idx).*beta(idx) + 4.*Bcp(idx).*beta(idx) + BB(idx).*betap(idx)))))./(Ah(idx).*(BB(idx) - Bc(idx)) + Bc(idx)).^4;
        d3betapardB3 = (-6.*(-1 + Ah).*Ah.^2.*Bc.*beta)./(Ah.*(BB - Bc) + Bc).^4;
        d3betapardB3(idx) = (-6.*(-1 + Ah(idx)).*Ah(idx).^2.*Bc(idx).*beta(idx))./(Ah(idx).*(BB(idx) - Bc(idx)) + Bc(idx)).^4;
        d3betapardBdR2 = zeros(size(RR));
        d3betapardB2dR = zeros(size(RR));
        d3betapardrdRdB = zeros(size(RR));
        if nargout > 13 % postprocessing
            betapar = (BB.*beta)./(Ah.*(BB - Bc) + Bc);
            betapar(idx) = (2.*Ah(idx).^2.*sqrt(Ah(idx).*(1 - BB(idx)./Bc(idx))).*(BB(idx) - Bc(idx)).^2 + Ah(idx).*(BB(idx) - Bc(idx)).*Bc(idx) - Bc(idx).^2)./(Bc(idx).*(Ah(idx).^2.*(BB(idx) - Bc(idx)).^2 - Bc(idx).^2));
            betaperp = betapar - BB.*(-(((-1 + Ah).*Bc.*beta)./(Ah.*(BB - Bc) + Bc).^2));
            betaperp(idx) = betapar(idx) - BB(idx).*(-(((-1 + Ah(idx)).*Bc(idx).*beta(idx))./(Ah(idx).*(BB(idx) - Bc(idx)) + Bc(idx)).^2));
        end
    end
end
