function LY = equil_SFL(L, LX, LY, r_plt, RR, ZZ, delta, deltap, deltapp, P_hat, Pp, Ppp, S, Sp, Spp, ms)
  LY.theta_SFL = LX.eps_val^(-2) * cumtrapz(LY.omega_plt, LY.JoverR ./RR,2) ./ r_plt;
  LY.theta_SFL(1,:) = LY.omega_plt;
  [LY.RR_sfl, LY.ZZ_sfl] = regrid_to_SFL(RR, ZZ, LY.theta_SFL, r_plt, LY.omega_plt);

  % compute N = gtt/Jac in SFL coordinates
  goo = (LX.eps_val.^3.*P_hat.*cos(LY.omega_plt) + LX.eps_val.*r_plt.*cos(LY.omega_plt) - LX.eps_val.^2.*sum((-1 + ms).*cos((-1 + ms).*LY.omega_plt).*S,3)).^2 + (-(LX.eps_val.^3.*P_hat.*sin(LY.omega_plt)) - LX.eps_val.*r_plt.*sin(LY.omega_plt) + LX.eps_val.^2.*sum(-((-1 + ms).*S.*sin((-1 + ms).*LY.omega_plt)),3)).^2;
  J_SFL = LX.eps_val * r_plt .* RR .^2; % to be multiplied by R0^2
  dthetaSFLdomega = LY.JoverR ./ RR./ r_plt * LX.eps_val^(-2);
  LY.dthetaSFLdomega= dthetaSFLdomega;
  gtt = (dthetaSFLdomega).^(-2) .* goo; % to be multiplied by R0^2
  LY.N_SFL = gtt./J_SFL; % currently, on the (r,\omega) grid, to be regridded if we want to pass it to VENUS-MHD
  LY.N0 = mean(LY.N_SFL(:,1:end-1)  .* dthetaSFLdomega(:,1:end-1), 2);
  LY.Nm1 = mean( exp(1i * LY.theta_SFL(:,1:end-1)).* LY.N_SFL(:,1:end-1)  .* dthetaSFLdomega(:,1:end-1), 2);
  % compute M = grt/Jac in SFL coordinates
  gro = (-(LX.eps_val.^3.*P_hat.*sin(LY.omega_plt)) - LX.eps_val.*r_plt.*sin(LY.omega_plt) + LX.eps_val.^2.*sum(-((-1 + ms).*S.*sin((-1 + ms).*LY.omega_plt)),3)).*(-(deltap.*LX.eps_val.^2) + LX.eps_val.*cos(LY.omega_plt) + LX.eps_val.^3.*Pp.*cos(LY.omega_plt) + LX.eps_val.^2.*sum(cos((-1 + ms).*LY.omega_plt).*Sp,3)) + (LX.eps_val.^3.*P_hat.*cos(LY.omega_plt) + LX.eps_val.*r_plt.*cos(LY.omega_plt) - LX.eps_val.^2.*sum((-1 + ms).*cos((-1 + ms).*LY.omega_plt).*S,3)).*(LX.eps_val.*sin(LY.omega_plt) + LX.eps_val.^3.*Pp.*sin(LY.omega_plt) - LX.eps_val.^2.*sum(sin((-1 + ms).*LY.omega_plt).*Sp,3));
  dthetaSFLdr = cumtrapz(LY.omega_plt,(RR.*(-((LX.eps_val.^2.*P_hat.*cos(LY.omega_plt) + r_plt.*cos(LY.omega_plt) - LX.eps_val.*sum((-1 + ms).*cos((-1 + ms).*LY.omega_plt).*S,3)).*(-(deltap.*LX.eps_val) + cos(LY.omega_plt) + LX.eps_val.^2.*Pp.*cos(LY.omega_plt) + LX.eps_val.*sum(cos((-1 + ms).*LY.omega_plt).*Sp,3))) - (LX.eps_val.^2.*P_hat.*sin(LY.omega_plt) + r_plt.*sin(LY.omega_plt) - LX.eps_val.*sum(-((-1 + ms).*S.*sin((-1 + ms).*LY.omega_plt)),3)).*(sin(LY.omega_plt) + LX.eps_val.^2.*Pp.*sin(LY.omega_plt) - LX.eps_val.*sum(sin((-1 + ms).*LY.omega_plt).*Sp,3))) - LX.eps_val.*r_plt.*(-(deltap.*LX.eps_val) + cos(LY.omega_plt) + LX.eps_val.^2.*Pp.*cos(LY.omega_plt) + LX.eps_val.*sum(cos((-1 + ms).*LY.omega_plt).*Sp,3)).*((LX.eps_val.^2.*P_hat.*cos(LY.omega_plt) + r_plt.*cos(LY.omega_plt) - LX.eps_val.*sum((-1 + ms).*cos((-1 + ms).*LY.omega_plt).*S,3)).*(-(deltap.*LX.eps_val) + cos(LY.omega_plt) + LX.eps_val.^2.*Pp.*cos(LY.omega_plt) + LX.eps_val.*sum(cos((-1 + ms).*LY.omega_plt).*Sp,3)) + (LX.eps_val.^2.*P_hat.*sin(LY.omega_plt) + r_plt.*sin(LY.omega_plt) - LX.eps_val.*sum(-((-1 + ms).*S.*sin((-1 + ms).*LY.omega_plt)),3)).*(sin(LY.omega_plt) + LX.eps_val.^2.*Pp.*sin(LY.omega_plt) - LX.eps_val.*sum(sin((-1 + ms).*LY.omega_plt).*Sp,3))) + r_plt.*RR.*((-(deltap.*LX.eps_val) + cos(LY.omega_plt) + LX.eps_val.^2.*Pp.*cos(LY.omega_plt) + LX.eps_val.*sum(cos((-1 + ms).*LY.omega_plt).*Sp,3)).*(cos(LY.omega_plt) + LX.eps_val.^2.*Pp.*cos(LY.omega_plt) - LX.eps_val.*sum((-1 + ms).*cos((-1 + ms).*LY.omega_plt).*Sp,3)) + (sin(LY.omega_plt) + LX.eps_val.^2.*Pp.*sin(LY.omega_plt) - LX.eps_val.*sum(sin((-1 + ms).*LY.omega_plt).*Sp,3)).*(sin(LY.omega_plt) + LX.eps_val.^2.*Pp.*sin(LY.omega_plt) - LX.eps_val.*sum(-((-1 + ms).*sin((-1 + ms).*LY.omega_plt).*Sp),3)) + LX.eps_val.*(LX.eps_val.^2.*P_hat.*cos(LY.omega_plt) + r_plt.*cos(LY.omega_plt) - LX.eps_val.*sum((-1 + ms).*cos((-1 + ms).*LY.omega_plt).*S,3)).*(-deltapp + LX.eps_val.*Ppp.*cos(LY.omega_plt) + sum(cos((-1 + ms).*LY.omega_plt).*Spp,3)) + LX.eps_val.*(LX.eps_val.^2.*P_hat.*sin(LY.omega_plt) + r_plt.*sin(LY.omega_plt) - LX.eps_val.*sum(-((-1 + ms).*S.*sin((-1 + ms).*LY.omega_plt)),3)).*(LX.eps_val.*Ppp.*sin(LY.omega_plt) - sum(sin((-1 + ms).*LY.omega_plt).*Spp,3))))./(r_plt.^2.*RR.^2),2);
  grt = dthetaSFLdomega.^(-1) .* (gro - goo .* dthetaSFLdr ./ dthetaSFLdomega) / LX.eps_val; % to be multiplied by R0


  LY.M_SFL = grt./J_SFL; % this is Mhat = M * R0;
  LY.Mm1 = mean( exp(1i * LY.theta_SFL(:,1:end-1)).* LY.M_SFL(:,1:end-1)  .* dthetaSFLdomega(:,1:end-1), 2);
  LY.Mm2 = mean( exp(2i * LY.theta_SFL(:,1:end-1)).* LY.M_SFL(:,1:end-1)  .* dthetaSFLdomega(:,1:end-1), 2);
  % Need NS > 2
  if numel(ms) > 2
      LY.Mm2_ana = - LX.eps_val * (3 * S(:,:,1) + r_plt .* (Sp(:,:,1) + r_plt .* Spp(:,:,1))) ./ (4 * r_plt) + ...
          LX.eps_val^2 * (1 ./ (24 .* r_plt)) .* (-12 * P_hat + 12 * (r_plt.^3 + S(:,:,1) .* Sp(:,:,1)) - 8 * S(:,:,2) .* ...
          (-8 * r_plt - 7 * deltap + r_plt .* deltapp - 3* Sp(:,:,2) + 3 * r_plt .* Spp(:,:,2)) + r_plt .* ...
          (12 * Pp + 12 * deltap .^2 + 12 *r_plt .* Ppp + 9 * r_plt.^2 .* deltapp + 25 * r_plt .* Sp(:,:,2) + 10 * r_plt .* deltapp .* Sp(:,:,2) - 24* Sp(:,:,2).^2 - 12 * ...
          (Sp(:,:,1).^2 + S(:,:,1) .* Spp(:,:,1)) + 9 * r_plt.^2 .* Spp(:,:,2) + r_plt .* deltap .* (33 + 18 *deltapp + 2* Spp(:,:,2))));
      LY.Mm2_ana_LO = - LX.eps_val * (3 * S(:,:,1) + r_plt .* (Sp(:,:,1) + r_plt .* Spp(:,:,1))) ./ (4 * r_plt);
      LY.Mm2_LO = - LX.eps_val * (3 * S(:,:,1) + r_plt .* (Sp(:,:,1) + r_plt .* Spp(:,:,1))) ./ (4 * r_plt);
  else
      LY.Mm2_ana = - LX.eps_val * (3 * S(:,:,1) + r_plt .* (Sp(:,:,1) + r_plt .* Spp(:,:,1))) ./ (4 * r_plt);
  end
  grr = LX.eps_val.^2.*((-(deltap.*LX.eps_val) + cos(LY.omega_plt) + LX.eps_val.^2.*Pp.*cos(LY.omega_plt) + LX.eps_val.*sum(cos((-1 + ms).*LY.omega_plt).*Sp,3)).^2 + (sin(LY.omega_plt) + LX.eps_val.^2.*Pp.*sin(LY.omega_plt) - LX.eps_val.*sum(sin((-1 + ms).*LY.omega_plt).*Sp,3)).^2);
  dw_dr = -dthetaSFLdr ./ dthetaSFLdomega;  % (∂ω/∂r̂)_ϑ, dimensionless
  grr_SFL = grr + 2 * gro/LX.eps_val .* dw_dr + goo .* dw_dr.^2;
  LY.L_SFL = grr_SFL ./ J_SFL * LX.eps_val^(-2);
  L1 = mean( exp(-1i * LY.theta_SFL(:,1:end-1)).* LY.L_SFL(:,1:end-1)  .* dthetaSFLdomega(:,1:end-1), 2);
  LY.L1 = L1;
  LY.L0 =  mean(LY.L_SFL(:,1:end-1)  .* dthetaSFLdomega(:,1:end-1), 2);
  % compute (j^\phi / B^\phi)_{-1}
  LY.JoverBm1 = mean( exp(-1i * LY.theta_SFL(:,1:end-1)).* LY.jphi(:,1:end-1) ./ sqrt(LY.BBt2(:,1:end-1)) .* dthetaSFLdomega(:,1:end-1), 2);
  LY.JoverBm1_2 = -LX.eps_val * LX.qfun(r_plt) .* LX.kinetic_profiles.betap(r_plt);
  LY.prefacY0 = mean(J_SFL(:,1:end-1).*sqrt(LY.BBt2(:,1:end-1))./ RR(:,1:end-1) .* dthetaSFLdomega(:,1:end-1), 2) ./ r_plt ./ LY.N0;

  l1 = r_plt .* (1 ./ LX.qfun(r_plt)-1);
  l1p = -1 + (1./LX.qfun(r_plt)) - r_plt .* LX.qpfun(r_plt) ./ LX.qfun(r_plt).^2;
  
  [~ ,iq1] = min(abs(LX.qfun(r_plt)-1));
  rs = r_plt(iq1);
  %Z1 = 1i *l1;
  %Y1 = -l1p;
  %LY.Y0 = - (r_plt < rs) ./LY.N0 .* ( LX.eps_val .* (1i .* LY.JoverBm1 + LY.Mm1) .* Z1 + LY.Nm1 .* Y1);
  %LY.Y0_LO = real((1i * LX.qfun(r_plt) .* LX.kinetic_profiles.betap(r_plt) .* Z1 ./ r_plt - (1i/2) * (deltapp + deltap ./ r_plt + 1) .* Z1 - Y1 .* deltap)) .* (r_plt < rs);
  
  % Updated Y0
  LY.Y0 = (r_plt < rs) .* (2.*LX.qfun(r_plt).^2.*r_plt - 2.*deltap.*LX.qpfun(r_plt).*r_plt + LX.qfun(r_plt).*(3.*deltap + (-1 + deltapp).*r_plt) - LX.qfun(r_plt).^3.*(3.*deltap + r_plt + deltapp.*r_plt))./(2..*LX.qfun(r_plt).^3);
  % THIS COULD BE MOVED TO EQUIL_ANA
  %LY.Y0_ana = (L.r_q < rs) .* (L.r_q + LX.qfun(L.r_q).^3 .* LX.kinetic_profiles.betap(L.r_q) - LX.qfun(L.r_q) .* (L.r_q + LX.kinetic_profiles.betap(L.r_q)) ...
  %   - L.r_q .* LX.qpfun(L.r_q) .* LY.deltap_ana) ./ LX.qfun(L.r_q);
  LY.Y0_ana = (L.r_q < rs) .* (LX.kinetic_profiles.betap(L.r_q).*LX.qfun(L.r_q).^3 + L.r_q - LY.deltap_ana.*LX.qpfun(L.r_q).*L.r_q - LX.qfun(L.r_q).*(LX.kinetic_profiles.betap(L.r_q) + L.r_q))./LX.qfun(L.r_q);
  %LY.Y0_LO = (r_plt < rs) .* (-r_plt +(LX.qfun(r_plt).^2-1).* LX.kinetic_profiles.betap(r_plt)  + 2 * deltap .* (1 + 2 * r_plt .* LX.qpfun(r_plt) ./ LX.qfun(r_plt).^2 ...
  %    ) + deltap .* (r_plt - (2 + r_plt .* LX.qpfun(r_plt)) ./ LX.qfun(r_plt) ) );
  % jac B^\phi for Z1 computation

  % Only for quadratic profiles!!
  betapp = -2*L.P.beta;
  qpp = L.P.q0 * L.P.s0;
  LY.dW0 = - trapz(L.r_q, (L.r_q < rs) .*(pi.*L.r_q.*(LX.kinetic_profiles.betap(L.r_q).^2.*LX.qfun(L.r_q).^6 + LY.deltap_ana.*LX.qpfun(L.r_q).*(1 + ...
LY.deltap_ana.*LX.qpfun(L.r_q)).*L.r_q.^2 + LX.qfun(L.r_q).^3.*L.r_q.*(-((1 + ...
betapp).*LY.deltap_ana) + LX.kinetic_profiles.betap(L.r_q).*(2 - 5.*LY.deltap_ana.*LX.qpfun(L.r_q)) + L.r_q) + ...
LX.qfun(L.r_q).^5.*(-(LX.kinetic_profiles.betap(L.r_q).*L.r_q) + betapp.*LY.deltap_ana.*L.r_q) - ...
LX.qfun(L.r_q).^4.*(LX.kinetic_profiles.betap(L.r_q).^2 + betapp.*LY.deltap_ana.*L.r_q - ...
4.*LX.kinetic_profiles.betap(L.r_q).*LY.deltap_ana.*LX.qpfun(L.r_q).*L.r_q) - LX.qfun(L.r_q).*L.r_q.*(-L.r_q + ...
LY.deltap_ana.*(1 + LX.qpfun(L.r_q).*L.r_q) + LY.deltap_ana.^2.*(2.*LX.qpfun(L.r_q) + LX.qpfun(L.r_q).^2.*L.r_q - ...
qpp.*L.r_q)) - LX.qfun(L.r_q).^2.*L.r_q.*(LX.kinetic_profiles.betap(L.r_q) - (2 + ...
betapp).*LY.deltap_ana + 2.*L.r_q + LY.deltap_ana.^2.*(-2.*LX.qpfun(L.r_q) + ...
qpp.*L.r_q))))./LX.qfun(L.r_q).^3);

  % V2 potential for xi2
  LY.V2 = (L.r_q < rs) .*((-2.*(-2 + LX.qfun(L.r_q)).*(4.*LX.kinetic_profiles.betap(L.r_q).*LX.qfun(L.r_q).^3.*LX.qpfun(L.r_q).*L.r_q ...
+ 2.*LY.deltap_ana.*LX.qpfun(L.r_q).^2.*L.r_q.^2 + LX.qfun(L.r_q).*(LX.qpfun(L.r_q) + LY.deltap_ana.*LX.qpfun(L.r_q).^2 - ...
2.*LY.deltap_ana.*qpp).*L.r_q.^2 + LX.qfun(L.r_q).^5.*(LX.kinetic_profiles.betap(L.r_q) - ...
betapp.*L.r_q) + LX.qfun(L.r_q).^2.*L.r_q.*(LY.deltap_ana.*qpp.*L.r_q + ...
LX.qpfun(L.r_q).*(-3.*LY.deltap_ana + L.r_q)) - 2.*LX.qfun(L.r_q).^4.*(-(betapp.*L.r_q) + ...
LX.kinetic_profiles.betap(L.r_q).*(1 + 2.*LX.qpfun(L.r_q).*L.r_q))))./LX.qfun(L.r_q).^4);
  
  LY.jacBphiup =  mean(sqrt(LY.BBt2(:,1:end-1)) ./ RR(:,1:end-1) .* J_SFL(:,1:end-1) .* dthetaSFLdomega(:,1:end-1), 2);

  % %% Compute dMdr, dNdr, d(1/B_0^\phi)/dr, d(J_0^\phi/B_0^\phi)dr
  % % Let's try finite differences for now, even though that is not very
  % % clean
  % M1 = mean( exp(-1i * LY.theta_SFL(:,1:end-1)).* LY.M_SFL(:,1:end-1)  .* dthetaSFLdomega(:,1:end-1), 2);
  % M1p = gradient(M1, r_plt);
  % N1 = mean( exp(-1i * LY.theta_SFL(:,1:end-1)).* LY.N_SFL(:,1:end-1)  .* dthetaSFLdomega(:,1:end-1), 2);
  % N1p = gradient(N1, r_plt);
  % JoverB1 = mean( exp(-1i * LY.theta_SFL(:,1:end-1)).* LY.jphi(:,1:end-1) ./ sqrt(LY.BBt2(:,1:end-1)) .* dthetaSFLdomega(:,1:end-1), 2);
  % JoverB1p = gradient(JoverB1, r_plt);
  % oneoverB1 = mean( exp(-1i * LY.theta_SFL(:,1:end-1)).* RR(:,1:end-1) ./ sqrt(LY.BBt2(:,1:end-1)) .* dthetaSFLdomega(:,1:end-1), 2);
  % oneoverB1p = gradient(oneoverB1, r_plt);
  % 
  % % Actually need q'' and beta'' ...
  % % might as well sin all the way.
  % l1pp = gradient(l1p, r_plt);
  % betapp = gradient(LX.kinetic_profiles.betap(r_plt), r_plt);
  % LY.V2 = real(- 2 * 1i * r_plt .* (- (2 ./ LX.qfun(r_plt) - 1) .* ((M1p / LX.eps_val .* l1 + M1 / LX.eps_val .* l1p + ...
  %     1i * N1p / LX.eps_val^2 .* l1p + 1i .* N1 / LX.eps_val^2 .* l1pp)) - 2 * 1i * ...
  %     (L1 .* l1 + 1i * M1 / LX.eps_val .* l1p) + ...
  %     1i *( JoverB1p /LX.eps_val .* l1 - JoverB1/LX.eps_val .* l1p) - ...
  %     1i * (-1/LX.eps_val * oneoverB1p .* LX.kinetic_profiles.betap(r_plt) + 1/LX.eps_val * oneoverB1p .* betapp)));
  % LY.V2(1:2) = 0;
  % LY.dxidrjump = real(2/ rs * (N1(iq1) / LX.eps_val^2 .* LX.qpfun(rs) - 1/ rs / LX.eps_val * oneoverB1(iq1) .* LX.kinetic_profiles.betap(rs))); 
end


function [R_sfl, Z_sfl] = regrid_to_SFL(RR, ZZ, theta_SFL, r_plt, omega_plt)
    Nr     = numel(r_plt);
    Ntheta = numel(omega_plt);
    theta_uniform = linspace(0, 2*pi, Ntheta);

    R_sfl = ones(Nr, Ntheta);
    Z_sfl = zeros(Nr, Ntheta);

    for i = 2 : Nr
        th_ext = [theta_SFL(i,end)-2*pi, theta_SFL(i,:), theta_SFL(i,1)+2*pi];
        R_ext  = [RR(i,end),             RR(i,:),        RR(i,1)];
        Z_ext  = [ZZ(i,end),             ZZ(i,:),        ZZ(i,1)];

        [th_ext, idx] = unique(th_ext);
        R_ext = R_ext(idx);
        Z_ext = Z_ext(idx);

        R_sfl(i,:) = interp1(th_ext, R_ext, theta_uniform, 'spline');
        Z_sfl(i,:) = interp1(th_ext, Z_ext, theta_uniform, 'spline');
    end
end
