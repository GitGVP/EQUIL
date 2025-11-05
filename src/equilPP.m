function LY = equilPP(L, LX, LY)
  % unpack solution vector
  idx = 0;
  t2_c    = LY.x(1:L.dof_count); idx = idx + L.dof_count;
  delta_c = LY.x(idx+1:idx+L.dof_count); idx = idx + L.dof_count;
  P_c     = LY.x(idx+1:idx+L.dof_count); idx = idx + L.dof_count;
  Bs_c    = reshape(LY.x(idx+1:idx+L.P.Nb*L.dof_count), [L.dof_count,1,L.P.Nb]); idx = idx + L.P.Nb*L.dof_count;
  S_unkn  = reshape(LY.x(idx+1:end), [L.dof_count-1,1,L.P.Ns]);

  % build S and B including boundary components
  S_c = cat(1, S_unkn, reshape(LX.Sbc, [1,1,L.P.Ns]));
  Bs_flat = reshape(Bs_c, size(Bs_c,1), []);  % dof_count x Nb  
  S_flat  = reshape(S_c,  size(S_c,1),  []);  % dof_count x Ns

  t2      = L.P0 * t2_c;
  t2p     = L.P1 * t2_c;
  delta   = L.P0 * delta_c;
  deltap  = L.P1 * delta_c;
  deltapp = L.P2 * delta_c;
  P_hat   = L.P0 * P_c;
  Pp      = L.P1 * P_c;
  Ppp     = L.P2 * P_c;

  Bs      = reshape(L.P0 * Bs_flat, L.Nq, 1, []);
  Bsp     = reshape(L.P1 * Bs_flat, L.Nq, 1, []);
  S       = reshape(L.P0 * S_flat, L.Nq, 1, []);
  Sp      = reshape(L.P1 * S_flat, L.Nq, 1, []);
  Spp     = reshape(L.P2 * S_flat, L.Nq, 1, []);
  ms = reshape(linspace(2,L.P.Ns+1,L.P.Ns), [1,1,L.P.Ns]);
  mb = reshape(linspace(1,L.P.Nb,L.P.Nb),   [1,1,L.P.Nb]);

  %% bp and li (MEQ conventions)
  RR = 1 - delta.*LX.eps_val.^2 + LX.eps_val.^3.*P_hat.*cos(L.omega) + LX.eps_val.*L.r_q.*cos(L.omega) + LX.eps_val.^2.*sum(cos((-1 + ms).*L.omega).*S,3);
  BB = 1 + LX.eps_val.*sum(Bs.*cos(mb.*L.omega),3);
  [dbetapardB0, ~, ~, ~, ~, ...
    ~, ~, ~, ~, ~, ...
    ~, ~, ~, ~, betapar, ...
    betaperp] = L.P.equation_of_state(LX.kinetic_profiles,L.r_q,RR,BB);
  % this is the averaged (p_par + p_perp)/2
  
  Intp = trapz(L.r_q,2*pi*mean(((betapar + betaperp).*RR.*((LX.eps_val.^2.*P_hat.*cos(L.omega) + L.r_q.*cos(L.omega) - LX.eps_val.*sum((-1 + ms).*cos((-1 + ms).*L.omega).*S,3)).*(-(deltap.*LX.eps_val) + cos(L.omega) + LX.eps_val.^2.*Pp.*cos(L.omega) + LX.eps_val.*sum(cos((-1 + ms).*L.omega).*Sp,3)) + (LX.eps_val.^2.*P_hat.*sin(L.omega) + L.r_q.*sin(L.omega) - LX.eps_val.*sum(-((-1 + ms).*S.*sin((-1 + ms).*L.omega)),3)).*(sin(L.omega) + LX.eps_val.^2.*Pp.*sin(L.omega) - LX.eps_val.*sum(sin((-1 + ms).*L.omega).*Sp,3))))./2.,2));  
  Intphir = 2*pi*mean((L.r_q.*(1 - dbetapardB0.*LX.eps_val.^2 + LX.eps_val.^2.*t2).*((LX.eps_val.^3.*P_hat.*cos(L.omega) + LX.eps_val.*L.r_q.*cos(L.omega) - LX.eps_val.^2.*sum((-1 + ms).*cos((-1 + ms).*L.omega).*S,3)).^2 + (LX.eps_val.^3.*P_hat.*sin(L.omega) + LX.eps_val.*L.r_q.*sin(L.omega) - LX.eps_val.^2.*sum(-((-1 + ms).*S.*sin((-1 + ms).*L.omega)),3)).^2))./(LX.eps_val.^2.*LX.qfun(L.r_q).*RR.*((LX.eps_val.^2.*P_hat.*cos(L.omega) + L.r_q.*cos(L.omega) - LX.eps_val.*sum((-1 + ms).*cos((-1 + ms).*L.omega).*S,3)).*(-(deltap.*LX.eps_val) + cos(L.omega) + LX.eps_val.^2.*Pp.*cos(L.omega) + LX.eps_val.*sum(cos((-1 + ms).*L.omega).*Sp,3)) + (LX.eps_val.^2.*P_hat.*sin(L.omega) + L.r_q.*sin(L.omega) - LX.eps_val.*sum(-((-1 + ms).*S.*sin((-1 + ms).*L.omega)),3)).*(sin(L.omega) + LX.eps_val.^2.*Pp.*sin(L.omega) - LX.eps_val.*sum(sin((-1 + ms).*L.omega).*Sp,3)))),2);  
  Intphi = Intphir(end);
  Intq = trapz(L.r_q,2*pi*mean((L.r_q.^2.*(1 - dbetapardB0.*LX.eps_val.^2 + LX.eps_val.^2.*t2).^2.*((LX.eps_val.^3.*P_hat.*cos(L.omega) + LX.eps_val.*L.r_q.*cos(L.omega) - LX.eps_val.^2.*sum((-1 + ms).*cos((-1 + ms).*L.omega).*S,3)).^2 + (LX.eps_val.^3.*P_hat.*sin(L.omega) + LX.eps_val.*L.r_q.*sin(L.omega) - LX.eps_val.^2.*sum(-((-1 + ms).*S.*sin((-1 + ms).*L.omega)),3)).^2))./(LX.eps_val.^2.*LX.qfun(L.r_q).^2.*RR.*((LX.eps_val.^2.*P_hat.*cos(L.omega) + L.r_q.*cos(L.omega) - LX.eps_val.*sum((-1 + ms).*cos((-1 + ms).*L.omega).*S,3)).*(-(deltap.*LX.eps_val) + cos(L.omega) + LX.eps_val.^2.*Pp.*cos(L.omega) + LX.eps_val.*sum(cos((-1 + ms).*L.omega).*Sp,3)) + (LX.eps_val.^2.*P_hat.*sin(L.omega) + L.r_q.*sin(L.omega) - LX.eps_val.*sum(-((-1 + ms).*S.*sin((-1 + ms).*L.omega)),3)).*(sin(L.omega) + LX.eps_val.^2.*Pp.*sin(L.omega) - LX.eps_val.*sum(sin((-1 + ms).*L.omega).*Sp,3)))),2));
  LY.intp = Intp;
  LY.intphi = Intphi;
  LY.intq = Intq;

  LY.bp_liu = 8 * pi * Intp / Intphi^2;
  LY.li_liu = 8 * pi * Intq / Intphi^2; % unsure about this 8, maybe 4.
  
  % augment with r=0 point where required
  r_plt = [0; L.r_q];
  delta = [0; delta];
  P_hat = [0; P_hat];

  % psiN
  %% ! here missing the anisotropy case!!
  psi = cumtrapz(L.r_q, (1 + LX.eps_val^2 * t2)./LX.q_vec .* L.r_q); % here factor epsilon might be missing
  psiN = (psi - psi(1)) / (psi(end) - psi(1)); % but it does not matter here

  % pad S and Bs with zero at index 0 to match r_plt
  S = cat(1, zeros(1,1,L.P.Ns), S);
  Bs = cat(1, zeros(1,1,L.P.Nb), Bs);
  Bsp = cat(1, zeros(1,1,L.P.Nb), Bsp);

  BB = 1 + LX.eps_val * sum(Bs .* cos(mb .* L.omega), 3);

  % ---- Geometry ----
  RR = 1 + LX.eps_val .* r_plt .* cos(L.omega) - LX.eps_val^2 .* delta + ...
    LX.eps_val^2 .* sum(S .* cos((ms-1) .* L.omega), 3) + ...
    LX.eps_val^3 .* P_hat .* cos(L.omega);
  ZZ = LX.eps_val .* r_plt .* sin(L.omega) - ...
    LX.eps_val^2 .* sum(S .* sin((ms-1) .* L.omega), 3) + ...
    LX.eps_val^3 .* P_hat .* sin(L.omega);

  % extrema and geometric measures
  [rmin, irmin] = min(RR, [], 2);
  [rmax, irmax] = max(RR, [], 2);
  [zmin, izmin] = min(ZZ, [], 2);
  [zmax, izmax] = max(ZZ, [], 2);

  Nr = numel(r_plt);
  zrmin = ZZ(sub2ind(size(ZZ), (1:Nr).', irmin));
  zrmax = ZZ(sub2ind(size(ZZ), (1:Nr).', irmax));
  rzmin = RR(sub2ind(size(RR), (1:Nr).', izmin));
  rzmax = RR(sub2ind(size(RR), (1:Nr).', izmax));

  kappa = (zmax - zmin) ./ (rmax - rmin);
  rgeom = (rmax + rmin)/2;
  zgeom = (zmax + zmin)/2;
  aminor = (rmax - rmin)/2;
  deltal = (rgeom - rzmin) ./ aminor;
  deltau = (rgeom - rzmax) ./ aminor;
  deltatrig = (deltal + deltau)/2;


  
  % ---- Package results ----
  LY.r_plt = r_plt;
  LY.t2 = t2;
  LY.t2p = t2p;
  LY.delta = delta;
  LY.deltap = deltap;
  LY.deltapp = deltapp;
  LY.P = P_hat;
  LY.Pp = Pp;
  LY.Bs = Bs;
  LY.Bsp = Bsp;
  LY.S = S;
  LY.Sp = Sp;
  LY.Spp = Spp;
  LY.BB = BB;
  LY.psiN = psiN;
  LY.RR = RR;
  LY.ZZ = ZZ;
  LY.kappa = kappa;
  LY.deltatrig = deltatrig;
  LY.aminor = aminor;
end