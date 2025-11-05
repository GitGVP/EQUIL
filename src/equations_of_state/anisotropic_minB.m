function [dbetapardr, dbetapardB, d2betapardB2, d2betapardrdB] = anisotropic_minB(beta,betap,Ah,Ahp,BB,dBBdr)
    Nr = size(BB,1);
    rowIdx = (1:Nr).';
    % Test to see if we can reproduce the MEQ observations
    
    [minB, idx] = min(BB, [], 2); % This acts as Bc in the modified bi-Maxwellian
    % no dminBdr as Bc is cst
    dminBdr = dBBdr(sub2ind(size(dBBdr), rowIdx, idx));

    % test to see if this is what matters
    
    minB = 1 - 0.3 ; %! Have to be careful to choose that this be less than B 
    % everywhere to avoid the discontinuity
    dminBdr = 0;
    d2betapardrdB = (2.*Ah.*BB.*beta.*(dminBdr - Ah.*dminBdr + Ahp.*(BB - ...
    minB)))./(Ah.*(BB - minB) + minB).^3 - (Ahp.*BB.*beta)./(Ah.*(BB - ...
    minB) + minB).^2 - (Ah.*BB.*betap)./(Ah.*(BB - minB) + minB).^2 - ...
    (beta.*(dminBdr - Ah.*dminBdr + Ahp.*(BB - minB)))./(Ah.*(BB - minB) ...
    + minB).^2 + betap./(Ah.*(BB - minB) + minB);

    dbetapardB = -((Ah.*BB.*beta)./(Ah.*(BB - minB) + minB).^2) + beta./(Ah.*(BB - minB) + minB);
    d2betapardB2 = (2.*Ah.^2.*BB.*beta)./(Ah.*(BB - minB) + minB).^3 - (2.*Ah.*beta)./(Ah.*(BB - minB) + minB).^2;
    dbetapardr = -((BB.*beta.*(dminBdr - Ah.*dminBdr + Ahp.*(BB - minB)))./(Ah.*(BB - minB) + minB).^2) + (BB.*betap)./(Ah.*(BB - minB) + minB);
end