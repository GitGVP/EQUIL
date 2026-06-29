function to_venus_infani(filename, L, LX, qpp, m, n, deltapress, deltapressp, dbetapardr, d2betapardr2, d2betapardB2, ...
    d2betapardrdB, d3betapardB2dr, d3betapardr2dB, d2betapardrdR)
%
 
    if exist(filename, 'file')
        delete(filename)
    end
    
    Nr = numel(L.r_q); Ntheta = 300-1;
    R0 = 1; B0 = 1;
    P0 = LX.eps_val^2*B0^2/4./pi/1.0E-07;
    Lref = 1; Bref = 1;
    
    
    
    %% GEOMETRY
    
    h5create(filename, '/geometry/R', [Nr, Ntheta]);
    h5write(filename, '/geometry/R', ones(Nr,Ntheta));
%     
%     h5create(filename, '/geometry/Z',[Nr, Ntheta]);
%     h5write(filename, '/geometry/Z', R0*LY.ZZ_sfl(:,1:end-1));
    
    %% NORMALISATION
    h5create(filename, '/normalisation/R0', 1);
    h5write(filename, '/normalisation/R0', R0);
    
    h5create(filename, '/normalisation/B0', 1);
    h5write(filename, '/normalisation/B0', B0);
    
    h5create(filename, '/normalisation/Bref', 1);
    h5write(filename, '/normalisation/Bref', Bref);

    h5create(filename, '/normalisation/Lref', 1);
    h5write(filename, '/normalisation/Lref', Lref);
    
    h5create(filename, '/normalisation/P0', 1);
    h5write(filename, '/normalisation/P0', P0);
    
    h5create(filename, '/normalisation/M02', 1);
    h5write(filename, '/normalisation/M02', 0);
    
    %% PROFILES

%     r=LY.r_fine;
%     beta_poloidal = LX.qfun(L.r_q).^2 ./ L.r_q .^ 3 .* cumtrapz(L.r_q.^2 .* (-2 * dbetapardr - d2betapardrdR + d2betapardrdB));
%     [~, I] = min(sqrt((LX.qfun(r)-1).^2));
%     beta_rs = beta_poloidal(I);

    deltap = LX.qfun(L.r_q).^2 ./ L.r_q .^ 3 .* cumtrapz(L.r_q, L.r_q.^3 ./ LX.qfun(L.r_q).^2 + L.r_q.^2 .* (-2 * dbetapardr - d2betapardrdR + d2betapardrdB));
    delta = cumtrapz(L.r_q, deltap);
    deltapp = 1 - 2 * LX.qfun(L.r_q).^2 .* LX.kinetic_profiles.betap(L.r_q) ./ L.r_q - deltap .* (3 ./L.r_q - 2 * LX.qpfun(L.r_q)./LX.qfun(L.r_q));
    
    
    h5create(filename, '/profiles/eps', 1);
    h5write(filename, '/profiles/eps', LX.eps_val);

    h5create(filename, '/profiles/mic', 1);
    h5write(filename, '/profiles/mic', m);

    h5create(filename, '/profiles/nic', 1);
    h5write(filename, '/profiles/nic', n);

%     h5create(filename, '/postprocessing/betap_rs', 1);
%     h5write(filename, '/postprocessing/betap_rs', beta_rs);
    
    h5create(filename, '/profiles/h', Nr);
    h5write(filename, '/profiles/h', L.r_q);
    
    h5create(filename, '/profiles/q', Nr);
    h5write(filename, '/profiles/q', LX.qfun(L.r_q));

    h5create(filename, '/profiles/qp', Nr);
    h5write(filename, '/profiles/qp', LX.qpfun(L.r_q));

    h5create(filename, '/profiles/qpp', Nr);
    h5write(filename, '/profiles/qpp', qpp(L.r_q));

    h5create(filename, '/profiles/delta', Nr);
    h5write(filename, '/profiles/delta', delta);

    h5create(filename, '/profiles/deltap', Nr);
    h5write(filename, '/profiles/deltap', deltap);
    
    h5create(filename, '/profiles/deltapp', Nr);
    h5write(filename, '/profiles/deltapp', deltapp);

    h5create(filename, '/profiles/dbetapardr', Nr);
    h5write(filename, '/profiles/dbetapardr', dbetapardr);

    h5create(filename, '/profiles/d2betapardr2', Nr);
    h5write(filename, '/profiles/d2betapardr2', d2betapardr2);

    h5create(filename, '/profiles/d2betapardB2', Nr);
    h5write(filename, '/profiles/d2betapardB2', d2betapardB2);

    h5create(filename, '/profiles/d2betapardrdB', Nr);
    h5write(filename, '/profiles/d2betapardrdB', d2betapardrdB);

    h5create(filename, '/profiles/d3betapardB2dr', Nr);
    h5write(filename, '/profiles/d3betapardB2dr', d3betapardB2dr);

    h5create(filename, '/profiles/d3betapardr2dB', Nr);
    h5write(filename, '/profiles/d3betapardr2dB', d3betapardr2dB);

    h5create(filename, '/profiles/deltapress', Nr);
    h5write(filename, '/profiles/deltapress', deltapress);

    h5create(filename, '/profiles/deltapressp', Nr);
    h5write(filename, '/profiles/deltapressp', deltapressp);


    fprintf('Done :)')



end