function LY = equil_ana(L, LX)
  r_fine = linspace(0,1,0.5*1e4); 
  r_fine = r_fine(2:end);
  [~, dbetapardr, ~, ~, ~, ...
d2betapardrdB, ~, d2betapardrdR, ~, d3betapardrdB2, ...
~, ~, ~, d3betapardrdRdB, d3betapardrdR2, ~, ...
~] = L.P.equation_of_state(LX.kinetic_profiles,r_fine.',ones(numel(r_fine), numel(L.omega)),ones(numel(r_fine), numel(L.omega)));
  dbetapardr = mean(dbetapardr,2).';
  d2betapardrdB = mean(d2betapardrdB,2).';
  d3betapardrdB2 = mean(d3betapardrdB2,2).';
  d3betapardrdRdB = mean(d3betapardrdRdB,2).';
  d2betapardrdR = mean(d2betapardrdR,2).';
  d3betapardrdR2 = mean(d3betapardrdR2,2).';

  % made anisotropic
  LY.t2p_fine = - 2 * r_fine ./ LX.qfun(r_fine) .^ 2 + r_fine .^ 2 .* LX.qpfun(r_fine) ./ LX.qfun(r_fine) .^ 3 - dbetapardr;
  LY.t2_fine =  cumtrapz(r_fine, LY.t2p_fine);
  LY.t2_ana = interp1(r_fine, LY.t2_fine, L.r_q, 'spline');
  LY.t2p_ana = interp1(r_fine, LY.t2p_fine, L.r_q, 'spline');

  LY.deltap_fine = LX.qfun(r_fine).^2 ./ r_fine .^ 3 .* cumtrapz(r_fine, r_fine.^3 ./ LX.qfun(r_fine).^2 + r_fine.^2 .* (-2 * dbetapardr - d2betapardrdR + d2betapardrdB));
  LY.delta_fine = cumtrapz(r_fine, LY.deltap_fine);
  LY.delta_ana = interp1(r_fine, LY.delta_fine, L.r_q, 'spline');
  LY.deltap_ana = interp1(r_fine, LY.deltap_fine, L.r_q, 'spline');


  [LY.S2_fine, LY.S2p_fine] = solve_S_eq(r_fine, LX.Sbc(1),LX.qfun,LX.qpfun, 2);
  [LY.S3_fine, LY.S3p_fine] = solve_S_eq(r_fine, LX.Sbc(2),LX.qfun,LX.qpfun, 3);
  LY.P_fine = -r_fine.^3/8 + LY.S2_fine .^2 ./ (2* r_fine) - r_fine .* LY.delta_fine / 2;  % todo: add other shaping
  LY.Pp_fine = -3*r_fine.^2/8 + 2 * LY.S2_fine .* LY.S2p_fine  ./ (2* r_fine) - LY.S2_fine.^2 ./ (2* r_fine.^2) -  LY.delta_fine / 2 - r_fine .* LY.deltap_fine / 2;
  LY.S2_ana = interp1(r_fine, LY.S2_fine, L.r_q, 'spline');
  LY.S3_ana = interp1(r_fine, LY.S3_fine, L.r_q, 'spline');
  LY.P_ana = interp1(r_fine, LY.P_fine, L.r_q, 'spline');
  LY.Pp_ana = interp1(r_fine, LY.Pp_fine, L.r_q, 'spline');

  %NLO shift
  delta1_RHS = (LY.S2_fine.*(-3.*(-d2betapardrdB + d2betapardrdR + 2.*dbetapardr).*LX.qfun(r_fine).^3.*r_fine + 4.*LX.qpfun(r_fine).*r_fine.*LY.S3_fine - 2.*LX.qfun(r_fine).*(-3.*LY.deltap_fine.*r_fine + 14.*LY.S3_fine + 5.*r_fine.*LY.S3p_fine)) + r_fine.*LY.S2p_fine.*(-3.*(-d2betapardrdB + d2betapardrdR + 2.*dbetapardr).*LX.qfun(r_fine).^3.*r_fine + 6.*LX.qpfun(r_fine).*r_fine.^2.*(LY.deltap_fine - LY.S3p_fine) + 2.*LX.qfun(r_fine).*(-10.*LY.S3_fine + r_fine.*(-3.*LY.deltap_fine + r_fine + LY.S3p_fine))))./(2..*LX.qfun(r_fine).^3);
  LY.delta1p_ana = LX.qfun(r_fine).^2 ./ r_fine .^ 3 .* cumtrapz(r_fine, delta1_RHS);
  LY.delta1_ana = cumtrapz(r_fine, LY.delta1p_ana);
  
  %NLO elongation
  if L.P.do_shift_NLO
      S2_1_RHS = (LX.qfun(r_fine).^2./r_fine.^3) .* (r_fine.*(6.*LY.deltap_fine.*LX.qpfun(r_fine).*r_fine.^2.*(LY.deltap_fine - 2.*LY.S3p_fine) + 2.*LX.qfun(r_fine).*(-6.*LY.deltap_fine.^2.*r_fine + r_fine.^3 + 2.*LY.deltap_fine.*(r_fine.^2 - 8.*LY.S3_fine) - 2.*r_fine.*(-2.*LY.deltap_fine + r_fine).*LY.S3p_fine) + LX.qfun(r_fine).^3.*(6.*d2betapardrdB.*LY.deltap_fine.*r_fine - 2.*d2betapardrdB.*r_fine.^2 + d3betapardrdR2.*r_fine.^2 + d3betapardrdB2.*r_fine.^2 - 2.*d3betapardrdRdB.*r_fine.^2 - 8.*d2betapardrdB.*LY.S3_fine - 6.*d2betapardrdB.*r_fine.*LY.S3p_fine + d2betapardrdR.*(-6.*LY.deltap_fine.*r_fine + 4.*r_fine.^2 + 8.*LY.S3_fine + 6.*r_fine.*LY.S3p_fine) + 2.*dbetapardr.*(8.*LY.S3_fine + r_fine.*(-6.*LY.deltap_fine + r_fine + 6.*LY.S3p_fine)))))./(4..*LX.qfun(r_fine).^3);
      S2_1_RHS_fun = @(r) interp1(r_fine, S2_1_RHS, r, 'spline');
      [LY.S2_1_fine, LY.S2_1p_fine] = solve_S_eq(r_fine, 0,LX.qfun,LX.qpfun, 2, S2_1_RHS_fun);
  end
  LY.r_fine = r_fine;
  LY.B1_ana =  - L.r_q;
end