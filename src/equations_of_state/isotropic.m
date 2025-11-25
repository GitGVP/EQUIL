function [dbetapardB0, dbetapardr, dbetapardB, dbetapardR, d2betapardB2, ...
    d2betapardrdB, d2betapardRdB, d2betapardrdR, d2betapardR2, d3betapardrdB2, ...
    d3betapardB3, d3betapardBdR2, d3betapardB2dR, d3betapardrdRdB, ...
    d3betapardrdR2, betapar, betaperp] = isotropic(kinetic_profiles,r,RR,BB)
    beta = kinetic_profiles.beta(r);
    betap = kinetic_profiles.betap(r);
    dbetapardr = betap;
    dbetapardB = zeros(size(RR));
    dbetapardR = zeros(size(RR));
    d2betapardB2 = zeros(size(RR));
    d2betapardrdB = zeros(size(RR));
    d2betapardRdB = zeros(size(RR));
    dbetapardB0 = dbetapardB(1,1);
    if nargout > 7 % Jacobian computation
        d2betapardrdR = zeros(size(RR));
        d2betapardR2 = zeros(size(RR));
        d3betapardrdB2 = zeros(size(RR));
        d3betapardB3 = zeros(size(RR));
        d3betapardBdR2 = zeros(size(RR));
        d3betapardB2dR = zeros(size(RR));
        d3betapardrdRdB = zeros(size(RR));
        d3betapardrdR2 = zeros(size(RR));
        if nargout > 14 % postprocessing
            betapar   = beta;
            betaperp  = beta;
        end
    end
end
