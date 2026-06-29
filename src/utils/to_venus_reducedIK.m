function to_venus_reducedIK(L, LX, LY, filename, qpp, m, n, ana)
%
 
    if exist(filename, 'file')
        delete(filename)
    end
    
    Nr = numel(LY.r_plt); Ntheta = numel(LY.omega_plt)-1;
    R0 = 1; B0 = 1;
%     R0 = 3.1923984520629225; B0 =  0.9617498446849172;
    P0 = LX.eps_val^2*B0^2/4./pi/1.0E-07;
    Lref = 1; Bref = 1;
    rho0=1;
    R0marg=3.0;

    % added constants 
    T0 = 1.0;
    
    
    %% GEOMETRY
    
    h5create(filename, '/geometry/R', [Nr, Ntheta]);
    h5write(filename, '/geometry/R', R0*LY.RR_sfl(:,1:end-1));
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

    betap = LX.kinetic_profiles.betap(LY.r_plt);
    betapp = LX.kinetic_profiles.betapp(LY.r_plt);

    r=LY.r_fine;
    beta_poloidal = LX.qfun(r).^2 ./ r .^ 4 .* cumtrapz(r, r.^2 .* (-2 * LX.kinetic_profiles.betap(r)));
    [~, I] = min(sqrt((LX.qfun(r)-1).^2));
    beta_rs = beta_poloidal(I);

    alpha = -2*LX.qfun(LY.r_plt).^2.*betap;
    alphap = -4*LX.qfun(LY.r_plt).*LX.qpfun(LY.r_plt).*betap - 2*LX.qfun(LY.r_plt).^2.*betapp;

    if (ana==true)
        delta = LY.delta_ana;
        deltap=LY.deltap_fine;
        deltapp = 1 - 2 * LX.qfun(LY.r_fine).^2 .* LX.kinetic_profiles.betap(LY.r_fine) ./ LY.r_fine - LY.deltap_fine .* (3 ./LY.r_fine - 2 * LX.qpfun(LY.r_fine)./LX.qfun(LY.r_fine));
        delta=interp1(L.r_q, delta, LY.r_plt, 'spline');
        deltap=interp1(LY.r_fine, deltap, LY.r_plt, 'spline');
        deltapp=interp1(LY.r_fine, deltapp, LY.r_plt, 'spline');
    else
        delta = cumtrapz(LY.r_plt, LY.deltap);
        deltap=LY.deltap;
        deltapp=LY.deltapp;
    end
    h5create(filename, '/profiles/eps', 1);
    h5write(filename, '/profiles/eps', LX.eps_val);

    h5create(filename, '/profiles/mic', 1);
    h5write(filename, '/profiles/mic', m);

    h5create(filename, '/profiles/nic', 1);
    h5write(filename, '/profiles/nic', n);

    h5create(filename, '/postprocessing/betap_rs', 1);
    h5write(filename, '/postprocessing/betap_rs', beta_rs);
    
    h5create(filename, '/profiles/h', Nr);
    h5write(filename, '/profiles/h', LY.r_plt);
    
    h5create(filename, '/profiles/q', Nr);
    h5write(filename, '/profiles/q', LX.qfun(LY.r_plt));

    h5create(filename, '/profiles/qp', Nr);
    h5write(filename, '/profiles/qp', LX.qpfun(LY.r_plt));

    h5create(filename, '/profiles/qpp', Nr);
    h5write(filename, '/profiles/qpp', qpp(LY.r_plt));

    h5create(filename, '/profiles/delta', Nr);
    h5write(filename, '/profiles/delta', delta);

    h5create(filename, '/profiles/deltap', Nr);
    h5write(filename, '/profiles/deltap', deltap);
    
    h5create(filename, '/profiles/deltapp', Nr);
    h5write(filename, '/profiles/deltapp', deltapp);

    h5create(filename, '/profiles/deltappp', Nr);
    h5write(filename, '/profiles/deltappp', LY.deltappp);

    h5create(filename, '/profiles/alpha', Nr);
    h5write(filename, '/profiles/alpha', alpha);

    h5create(filename, '/profiles/alphap', Nr);
    h5write(filename, '/profiles/alphap', alphap);

    h5create(filename, '/profiles/betap', Nr);
    h5write(filename, '/profiles/betap', betap);

    h5create(filename, '/profiles/betapp', Nr);
    h5write(filename, '/profiles/betapp', betapp);

    h5create(filename, '/profiles/P', Nr);
    h5write(filename, '/profiles/P', LY.P);

    h5create(filename, '/profiles/Pp', Nr);
    h5write(filename, '/profiles/Pp', LY.Pp);

    h5create(filename, '/profiles/Ppp', Nr);
    h5write(filename, '/profiles/Ppp', LY.Ppp);

    h5create(filename, '/profiles/Pppp', Nr);
    h5write(filename, '/profiles/Pppp', LY.Pppp);

    h5create(filename, '/profiles/t2p', Nr);
    h5write(filename, '/profiles/t2p', LY.t2p);

    h5create(filename, '/profiles/t2pp', Nr);
    h5write(filename, '/profiles/t2pp', LY.t2pp);

    h5create(filename, '/profiles/epsmarg', Nr);
    h5write(filename, '/profiles/epsmarg', LY.r_plt/R0);

    h5create(filename, '/profiles/depsdsmarg', 1);
    h5write(filename, '/profiles/depsdsmarg', 1/R0);

    fprintf('Done :)')

    %% POST PROCESSING
    
%     h5create(filename, '/postprocessing/S2bc', 1);
%     h5write(filename, '/postprocessing/S2bc', LX.Sbc(1));
% 
%     h5create(filename, '/postprocessing/S3bc', 1);
%     h5write(filename, '/postprocessing/S3bc', LX.Sbc(2));

%     h5create(filename, '/postprocessing/beta_poloidal', numel(beta_poloidal));
%     h5write(filename, '/postprocessing/beta_poloidal', beta_poloidal);
%     
%     h5create(filename, '/postprocessing/betap_rs', 1);
%     h5write(filename, '/postprocessing/betap_rs', beta_rs);

%     h5create(filename, '/postprocessing/betap_rsana', 1);
%     h5write(filename, '/postprocessing/betap_rsana', beta_rs_ana);

%     h5create(filename, '/postprocessing/kappa', 1);
%     h5write(filename, '/postprocessing/kappa', LY.kappa(end));

%     h5create(filename, '/postprocessing/kappa1', 1);
%     h5write(filename, '/postprocessing/kappa1', LY.kappa(Ii));
% 
%     h5create(filename, '/postprocessing/delta1', 1);
%     h5write(filename, '/postprocessing/delta1', LY.deltatrig(Ii));
% 
%     h5create(filename, '/postprocessing/deltatrig', 1);
%     h5write(filename, '/postprocessing/deltatrig', LY.deltatrig(end));
% 
%     h5create(filename, '/postprocessing/betap', Nr);
%     h5write(filename, '/postprocessing/betap', LX.kinetic_profiles.betap(LY.r_plt));
%     
%     h5create(filename, '/postprocessing/psiN', numel(LY.psiN));
%     h5write(filename, '/postprocessing/psiN', LY.psiN);
% 
%     h5create(filename, '/postprocessing/r_fine', numel(LY.r_fine));
%     h5write(filename, '/postprocessing/r_fine', LY.r_fine);


end