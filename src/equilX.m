function LX = equilX(L)
 % default inputs
  LX.eps_val = 0.3;

  LX.Sbc = zeros(L.P.Ns,1);
  LX.S2bc = -1;
  LX.Sbc(1) = LX.S2bc;
  LX.S3bc = 0.3;
  if L.P.Ns >1
    LX.Sbc(2) = LX.S3bc;
  end
  % %% Add a pedestal:
  r_ped = 0.6;
  beta_ped = 0;
  betapedfun = @(r) 0.5 * beta_ped * (1-tanh(4 * pi * (r - r_ped))); % this is wrong for now, but does not appear in isotropic
  betappedfun = @(r) - 4 * pi * beta_ped ./ (1+ cosh(8*pi*(r-r_ped)));

  LX.qfun = @(r) L.P.q0 * (1+0.5*L.P.s0*r.^2);
  LX.qpfun = @(r) L.P.q0*L.P.s0*r;
  
  
  betafun = @(r) L.P.beta * (1 - r.^2) + betapedfun(r);
  betapfun = @(r) -2 * L.P.beta * r + betappedfun(r);

% % test anisotropic, constant beta profile
%betafun = @(r) beta * ones(size(r));
%betapfun = @(r) zeros(size(r));

  LX.q_vec = LX.qfun(L.r_q);
  LX.qp_vec = LX.qpfun(L.r_q);
  %LX.beta_vec = LX.betafun(L.r_q);
  %LX.betap_vec = LX.betapfun(L.r_q);
  %LX.Ah_vec = LX.Ahfun(L.r_q);
  %LX.Ahp_vec = LX.Ahpfun(L.r_q);
  switch func2str(L.P.equation_of_state)
      case 'isotropic'
          LX.kinetic_profiles = struct('beta', betafun,'betap', betapfun);
      case 'parallel_rotating' 
          mach2fun = @(r) L.P.mach20;
          mach2pfun = @(r) 0;
          LX.kinetic_profiles = struct('beta', betafun,'betap', betapfun,...
              'mach2',mach2fun, 'mach2p',mach2pfun);
      case 'biMaxwellian'
          Bcfun = @(r) L.P.Bc0*ones(size(r));
          Bcpfun = @(r) zeros(size(r));
          Ahfun = @(r) L.P.A0*ones(size(r));
          Ahpfun = @(r) zeros(size(r));
          LX.kinetic_profiles = struct('beta', betafun,'betap', betapfun,...
              'Ah',Ahfun, 'Ahp',Ahpfun, 'Bc', Bcfun, 'Bcp',Bcpfun);
      otherwise
          fprintf('Equation of state not implemented!\n')
  end
end
