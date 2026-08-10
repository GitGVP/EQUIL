function [r_q, ...
          P0, P1, P2, P3, ...
          M_profiles, M_extended, P_templates, P_extended, ...
          A_global, profile_lengths, profile_starts, P0_end, ...
          P0_axis, P1_axis, P0_edge, P1_edge, w_q, ...
          axis_jet_scales] = ...
          assemble_FE_matrices_polynomial_axis_regular( ...
              r_nodes, Nb, Ns, nq, p, max_power, Nh, parity_basis)
%ASSEMBLE_FE_MATRICES_POLYNOMIAL_AXIS_REGULAR Axis-regular polynomial basis.
%
% This is an opt-in experimental alternative to
% assemble_FE_matrices_bspline_neumann.  It uses global polynomials with
% the desired magnetic-axis powers built into the trial space:
%
%   t2, delta : r^2 v(r^2)       (even powers only)
%   P         : r^3 v(r^2)
%   B_m       : r   v(r)
%   S_m, V_m  : r^(m-1) v(r^2),   m = 2,...
%
% The maximum radial power is m+p, matching the number of basis functions
% in the standard open-knot spline space.  Shifted Legendre polynomials are
% used instead of raw monomials, followed by a weighted QR change of basis.
%
% For S_m and V_m, a distinguished axis-jet column is followed by bubbles
% with zero leading Frobenius coefficient. The final fixed column supplies
% the exact edge value. Set parity_basis=false for the legacy all-power
% space used by the standard solver.

m = (numel(r_nodes)-1)/2;
if m ~= round(m) || m < 1
    error('r_nodes must contain 2*m+1 entries.');
end

% Use the same number of radial quadrature points as the standard path.
[gauss_pts, gauss_wts] = lgwt(nq, 0, 1);
gauss_pts = flip(gauss_pts(:));
gauss_wts = flip(gauss_wts(:));

Nq = nq*m;
r_q = zeros(Nq,1);
w_q = zeros(Nq,1);
for e = 1:m
    r0 = (e-1)/m;
    r1 = e/m;
    rows = (e-1)*nq + (1:nq);
    r_q(rows) = r0 + (r1-r0)*gauss_pts;
    w_q(rows) = (r1-r0)*gauss_wts;
end

if nargin < 6 || isempty(max_power)
    max_power = m+p;
end
if nargin < 7 || isempty(Nh)
    Nh = 0;
end
if nargin < 8 || isempty(parity_basis)
    parity_basis = false;
end
minimum_power = max([3,Ns+1,Nh+1]);
if parity_basis
    minimum_power = max([3,Ns+2,Nh+2]);
end
if max_power ~= round(max_power) || max_power < minimum_power
    error(['axis_regular_max_power must be an integer >= ', ...
           'the largest required axis/edge power.']);
end
n_even = floor(max_power/2);       % r^2,r^4,...
n_P = max_power-2;
n_B = max_power;                   % r^1,...,r^max_power

[Nr0, Nr1, Nr2, Nr3] = shifted_legendre_all(r_q,n_B);
[Ne0, Ne1, Ne2, Ne3] = shifted_legendre_all(r_q.^2,n_even);

if parity_basis
    n_P = floor((max_power-1)/2);
    [Pe0, Pe1, Pe2, Pe3] = parity_axis_space( ...
        Ne0,Ne1,Ne2,Ne3,r_q,2,n_even);
    [Pp0, Pp1, Pp2, Pp3] = parity_axis_space( ...
        Ne0,Ne1,Ne2,Ne3,r_q,3,n_P);
else
    [Pe0, Pe1, Pe2, Pe3] = even_r2_basis( ...
        Ne0(:,1:n_even),Ne1(:,1:n_even), ...
        Ne2(:,1:n_even),Ne3(:,1:n_even),r_q);
    [Pp0, Pp1, Pp2, Pp3] = radial_power_basis( ...
        Nr0(:,1:n_P),Nr1(:,1:n_P),Nr2(:,1:n_P),Nr3(:,1:n_P),r_q,3);
end
[Pb0, Pb1, Pb2, Pb3] = radial_power_basis( ...
    Nr0(:,1:n_B),Nr1(:,1:n_B),Nr2(:,1:n_B),Nr3(:,1:n_B),r_q,1);

nProfiles = 3 + Nb + Ns + Nh;
axis_jet_scales = nan(nProfiles,1);
P0 = cell(1,nProfiles);
P1 = cell(1,nProfiles);
P2 = cell(1,nProfiles);
P3 = cell(1,nProfiles);

% [t2 | delta | P]
P0(1:3) = {sparse(Pe0), sparse(Pe0), sparse(Pp0)};
P1(1:3) = {sparse(Pe1), sparse(Pe1), sparse(Pp1)};
P2(1:3) = {sparse(Pe2), sparse(Pe2), sparse(Pp2)};
P3(1:3) = {sparse(Pe3), sparse(Pe3), sparse(Pp3)};

% B_0,...,B_(Nb-1): all start at r^1.
for ib = 1:Nb
    alpha = 3 + ib;
    P0{alpha} = sparse(Pb0);
    P1{alpha} = sparse(Pb1);
    P2{alpha} = sparse(Pb2);
    P3{alpha} = sparse(Pb3);
end

% S(:,:,is) corresponds to m=is+1 and therefore starts at r^is.  Unknown
% columns are r^is*(1-r)*polynomial and hence vanish at the outer edge;
% the final fixed column r^is carries Sbc.
for is = 1:Ns
    alpha = 3 + Nb + is;
    if parity_basis
        [P0{alpha},P1{alpha},P2{alpha},P3{alpha}] = ...
            parity_shape_space(Ne0,Ne1,Ne2,Ne3,r_q,is,max_power);
    else
        n_S_unknown = max_power-is;
        cols = 1:n_S_unknown;
        [Sk0,Sk1,Sk2,Sk3] = radial_power_basis( ...
            Nr0(:,cols),Nr1(:,cols),Nr2(:,cols),Nr3(:,cols),r_q,is);
        [Skp0,Skp1,Skp2,Skp3] = radial_power_basis( ...
            Nr0(:,cols),Nr1(:,cols),Nr2(:,cols),Nr3(:,cols),r_q,is+1);
        [F0,F1,F2,F3] = pure_power_column(r_q,is);
        P0{alpha} = sparse([Sk0-Skp0,F0]);
        P1{alpha} = sparse([Sk1-Skp1,F1]);
        P2{alpha} = sparse([Sk2-Skp2,F2]);
        P3{alpha} = sparse([Sk3-Skp3,F3]);
    end
end

% V(:,:,iv) corresponds to Fitzpatrick V_(iv+1) and has the same radial
% regularity and fixed-edge construction as S_(iv+1).
for iv = 1:Nh
    alpha = 3 + Nb + Ns + iv;
    if parity_basis
        [P0{alpha},P1{alpha},P2{alpha},P3{alpha}] = ...
            parity_shape_space(Ne0,Ne1,Ne2,Ne3,r_q,iv,max_power);
    else
        n_V_unknown = max_power-iv;
        cols = 1:n_V_unknown;
        [Vk0,Vk1,Vk2,Vk3] = radial_power_basis( ...
            Nr0(:,cols),Nr1(:,cols),Nr2(:,cols),Nr3(:,cols),r_q,iv);
        [Vkp0,Vkp1,Vkp2,Vkp3] = radial_power_basis( ...
            Nr0(:,cols),Nr1(:,cols),Nr2(:,cols),Nr3(:,cols),r_q,iv+1);
        [F0,F1,F2,F3] = pure_power_column(r_q,iv);
        P0{alpha} = sparse([Vk0-Vkp0,F0]);
        P1{alpha} = sparse([Vk1-Vkp1,F1]);
        P2{alpha} = sparse([Vk2-Vkp2,F2]);
        P3{alpha} = sparse([Vk3-Vkp3,F3]);
    end
end

% Exact endpoint evaluation before the coefficient-space QR transform.
rb = [0;1];
[Nr0b, Nr1b, Nr2b, Nr3b] = shifted_legendre_all(rb,n_B);
[Ne0b, Ne1b, Ne2b, Ne3b] = shifted_legendre_all(rb.^2,n_even);
if parity_basis
    [Pe0b, Pe1b] = parity_axis_space( ...
        Ne0b,Ne1b,Ne2b,Ne3b,rb,2,n_even);
    [Pp0b, Pp1b] = parity_axis_space( ...
        Ne0b,Ne1b,Ne2b,Ne3b,rb,3,n_P);
else
    [Pe0b, Pe1b] = even_r2_basis( ...
        Ne0b(:,1:n_even),Ne1b(:,1:n_even), ...
        Ne2b(:,1:n_even),Ne3b(:,1:n_even),rb);
    [Pp0b, Pp1b] = radial_power_basis( ...
        Nr0b(:,1:n_P),Nr1b(:,1:n_P),Nr2b(:,1:n_P),Nr3b(:,1:n_P),rb,3);
end
[Pb0b, Pb1b] = radial_power_basis( ...
    Nr0b(:,1:n_B),Nr1b(:,1:n_B),Nr2b(:,1:n_B),Nr3b(:,1:n_B),rb,1);

P0b = cell(1,nProfiles);
P1b = cell(1,nProfiles);
P0b(1:3) = {Pe0b,Pe0b,Pp0b};
P1b(1:3) = {Pe1b,Pe1b,Pp1b};
for ib = 1:Nb
    alpha = 3+ib;
    P0b{alpha} = Pb0b;
    P1b{alpha} = Pb1b;
end
for is = 1:Ns
    alpha = 3+Nb+is;
    if parity_basis
        [P0b{alpha},P1b{alpha}] = ...
            parity_shape_space(Ne0b,Ne1b,Ne2b,Ne3b,rb,is,max_power);
    else
        n_S_unknown = max_power-is;
        cols = 1:n_S_unknown;
        [Sk0,Sk1] = radial_power_basis( ...
            Nr0b(:,cols),Nr1b(:,cols),Nr2b(:,cols),Nr3b(:,cols),rb,is);
        [Skp0,Skp1] = radial_power_basis( ...
            Nr0b(:,cols),Nr1b(:,cols),Nr2b(:,cols),Nr3b(:,cols),rb,is+1);
        [F0,F1] = pure_power_column(rb,is);
        P0b{alpha} = [Sk0-Skp0,F0];
        P1b{alpha} = [Sk1-Skp1,F1];
    end
end
for iv = 1:Nh
    alpha = 3+Nb+Ns+iv;
    if parity_basis
        [P0b{alpha},P1b{alpha}] = ...
            parity_shape_space(Ne0b,Ne1b,Ne2b,Ne3b,rb,iv,max_power);
    else
        n_V_unknown = max_power-iv;
        cols = 1:n_V_unknown;
        [Vk0,Vk1] = radial_power_basis( ...
            Nr0b(:,cols),Nr1b(:,cols),Nr2b(:,cols),Nr3b(:,cols),rb,iv);
        [Vkp0,Vkp1] = radial_power_basis( ...
            Nr0b(:,cols),Nr1b(:,cols),Nr2b(:,cols),Nr3b(:,cols),rb,iv+1);
        [F0,F1] = pure_power_column(rb,iv);
        P0b{alpha} = [Vk0-Vkp0,F0];
        P1b{alpha} = [Vk1-Vkp1,F1];
    end
end

% The fixed final S/V column is not part of the nonlinear state.
profile_lengths = cellfun(@(A)size(A,2),P0);
profile_lengths((4+Nb):end) = profile_lengths((4+Nb):end)-1;
profile_starts = cumsum([1, profile_lengths(1:end-1)]);
totalDofs = sum(profile_lengths);

% In the parity basis, only the first column has a nonzero physical axis
% jet. Higher bubbles are QR-scaled, the first column is made orthogonal
% using those bubbles, and its normalization is recorded explicitly.
for alpha = 1:nProfiles
    n_unknown = profile_lengths(alpha);
    preserve_axis_jet = parity_basis ...
        && (alpha <= 3 || alpha > 3+Nb);
    if preserve_axis_jet
        columns = 2:n_unknown;
    else
        columns = 1:n_unknown;
    end
    if ~isempty(columns)
        weighted_P0 = sqrt(w_q).*full(P0{alpha}(:,columns));
        [~,Rscale] = qr(weighted_P0,0);
        if rcond(Rscale) < 10*eps
            error(['Rank-deficient axis-regular polynomial space ', ...
                   'in profile %d.'],alpha);
        end
        Tscale = Rscale\eye(numel(columns));
        P0{alpha}(:,columns) = P0{alpha}(:,columns)*Tscale;
        P1{alpha}(:,columns) = P1{alpha}(:,columns)*Tscale;
        P2{alpha}(:,columns) = P2{alpha}(:,columns)*Tscale;
        P3{alpha}(:,columns) = P3{alpha}(:,columns)*Tscale;
        P0b{alpha}(:,columns) = P0b{alpha}(:,columns)*Tscale;
        P1b{alpha}(:,columns) = P1b{alpha}(:,columns)*Tscale;
    end
    if preserve_axis_jet
        if ~isempty(columns)
            projection = full(P0{alpha}(:,columns)' ...
                *(w_q.*P0{alpha}(:,1)));
            P0{alpha}(:,1) = P0{alpha}(:,1) ...
                -P0{alpha}(:,columns)*projection;
            P1{alpha}(:,1) = P1{alpha}(:,1) ...
                -P1{alpha}(:,columns)*projection;
            P2{alpha}(:,1) = P2{alpha}(:,1) ...
                -P2{alpha}(:,columns)*projection;
            P3{alpha}(:,1) = P3{alpha}(:,1) ...
                -P3{alpha}(:,columns)*projection;
            P0b{alpha}(:,1) = P0b{alpha}(:,1) ...
                -P0b{alpha}(:,columns)*projection;
            P1b{alpha}(:,1) = P1b{alpha}(:,1) ...
                -P1b{alpha}(:,columns)*projection;
        end
        axis_norm = sqrt(full(P0{alpha}(:,1)' ...
            *(w_q.*P0{alpha}(:,1))));
        axis_scale = 1/axis_norm;
        P0{alpha}(:,1) = P0{alpha}(:,1)*axis_scale;
        P1{alpha}(:,1) = P1{alpha}(:,1)*axis_scale;
        P2{alpha}(:,1) = P2{alpha}(:,1)*axis_scale;
        P3{alpha}(:,1) = P3{alpha}(:,1)*axis_scale;
        P0b{alpha}(:,1) = P0b{alpha}(:,1)*axis_scale;
        P1b{alpha}(:,1) = P1b{alpha}(:,1)*axis_scale;
        axis_jet_scales(alpha) = axis_scale;
    end
end

P0_end = zeros(Nq,Ns+Nh);
for ishape = 1:(Ns+Nh)
    alpha = 3 + Nb + ishape;
    P0_end(:,ishape) = P0{alpha}(:,end);
end

P0_axis = cellfun(@(A)sparse(A(1,:)),P0b,'UniformOutput',false);
P1_axis = cellfun(@(A)sparse(A(1,:)),P1b,'UniformOutput',false);
P0_edge = cellfun(@(A)sparse(A(2,:)),P0b,'UniformOutput',false);
P1_edge = cellfun(@(A)sparse(A(2,:)),P1b,'UniformOutput',false);

% Unknown-only trial matrices. The full S/V matrices keep their final,
% prescribed column for state reconstruction, but that column is absent
% from Newton derivatives and test functions.
Pu = cell(4,nProfiles);
Pall = {P0,P1,P2,P3};
for d = 1:4
    for alpha = 1:nProfiles
        Pu{d,alpha} = Pall{d}{alpha}(:,1:profile_lengths(alpha));
    end
end

W = spdiags(w_q,0,Nq,Nq);
M_blocks = cell(1,nProfiles);
A_blocks = cell(1,nProfiles);
for alpha = 1:nProfiles
    M_blocks{alpha} = Pu{1,alpha}'*W;
    A_blocks{alpha} = M_blocks{alpha}*Pu{1,alpha};
end
M_profiles = blkdiag(M_blocks{:});
A_global = blkdiag(A_blocks{:});

% M_extended acts on blocks ordered as (residual alpha, variable gamma,
% quadrature q), which is the ordering used by the current Jacobians.
nnz_M = nProfiles*sum(cellfun(@nnz,M_blocks));
Mi = zeros(nnz_M,1);
Mj = zeros(nnz_M,1);
Mv = zeros(nnz_M,1);
ptr = 1;
for alpha = 1:nProfiles
    [ii,jj,vv] = find(M_blocks{alpha});
    block_nnz = numel(vv);
    row0 = profile_starts(alpha)-1;
    for gamma = 1:nProfiles
        take = ptr:(ptr+block_nnz-1);
        col0 = ((alpha-1)*nProfiles + (gamma-1))*Nq;
        Mi(take) = row0 + ii;
        Mj(take) = col0 + jj;
        Mv(take) = vv;
        ptr = ptr + block_nnz;
    end
end
M_extended = sparse(Mi,Mj,Mv,totalDofs,nProfiles^2*Nq);

% Derivative projection templates.  Only P0/P1/P2 are used by the current
% residual Jacobians; P3 is retained for postprocessing.
P_templates = cell(1,3);
P_extended = cell(1,3);
for d = 1:3
    nnz_P = nProfiles*sum(cellfun(@nnz,Pu(d,:)));
    Pi = zeros(nnz_P,1);
    Pj = zeros(nnz_P,1);
    Pv = zeros(nnz_P,1);
    ptr = 1;
    for gamma = 1:nProfiles
        [ii,jj,vv] = find(Pu{d,gamma});
        block_nnz = numel(vv);
        col0 = profile_starts(gamma)-1;
        for alpha = 1:nProfiles
            take = ptr:(ptr+block_nnz-1);
            row0 = ((alpha-1)*nProfiles + (gamma-1))*Nq;
            Pi(take) = row0 + ii;
            Pj(take) = col0 + jj;
            Pv(take) = vv;
            ptr = ptr + block_nnz;
        end
    end
    P_extended{d} = sparse(Pi,Pj,Pv,nProfiles^2*Nq,totalDofs);
    P_templates{d} = struct('i',Pi,'j',Pj,'v_template',Pv);
end

end


function [F0,F1,F2,F3] = even_r2_basis(N0,N1,N2,N3,r)
% f(r)=r^2*N(r^2), where the N derivatives are with respect to s=r^2.
r = r(:);
F0 = r.^2.*N0;
F1 = 2*r.*N0 + 2*r.^3.*N1;
F2 = 2*N0 + 10*r.^2.*N1 + 4*r.^4.*N2;
F3 = 24*r.*N1 + 36*r.^3.*N2 + 8*r.^5.*N3;
end


function [F0,F1,F2,F3] = parity_axis_space(N0,N1,N2,N3,r,k,ncols)
% r^k*A plus parity-correct bubbles r^(k+2)*Phi_j(r^2).
[A0,A1,A2,A3] = pure_power_column(r,k);
if ncols <= 1
    F0 = sparse(A0); F1 = sparse(A1);
    F2 = sparse(A2); F3 = sparse(A3);
    return
end
cols = 1:(ncols-1);
[B0,B1,B2,B3] = radial_power_s_basis( ...
    N0(:,cols),N1(:,cols),N2(:,cols),N3(:,cols),r,k+2);
F0 = sparse([A0,B0]);
F1 = sparse([A1,B1]);
F2 = sparse([A2,B2]);
F3 = sparse([A3,B3]);
end


function [F0,F1,F2,F3] = parity_shape_space( ...
        N0,N1,N2,N3,r,n,max_power)
% r^n[A*(1-r^2)+bc*r^2+r^2*(1-r^2)*sum c_j Phi_j(r^2)].
[An0,An1,An2,An3] = pure_power_column(r,n);
[E0,E1,E2,E3] = pure_power_column(r,n+2);
axis0 = An0-E0; axis1 = An1-E1;
axis2 = An2-E2; axis3 = An3-E3;
n_bubbles = max(0,floor((max_power-n-2)/2));
if n_bubbles > 0
    cols = 1:n_bubbles;
    [B0,B1,B2,B3] = radial_power_s_basis( ...
        N0(:,cols),N1(:,cols),N2(:,cols),N3(:,cols),r,n+2);
    [C0,C1,C2,C3] = radial_power_s_basis( ...
        N0(:,cols),N1(:,cols),N2(:,cols),N3(:,cols),r,n+4);
    B0 = B0-C0; B1 = B1-C1;
    B2 = B2-C2; B3 = B3-C3;
else
    B0 = zeros(numel(r),0); B1 = B0;
    B2 = B0; B3 = B0;
end
F0 = sparse([axis0,B0,E0]);
F1 = sparse([axis1,B1,E1]);
F2 = sparse([axis2,B2,E2]);
F3 = sparse([axis3,B3,E3]);
end


function [F0,F1,F2,F3] = radial_power_s_basis(N0,N1,N2,N3,r,k)
% f(r)=r^k*N(r^2); N1,N2,N3 are derivatives with respect to s=r^2.
r = r(:);
F0 = r.^k.*N0;
F1 = 2*r.^(k+1).*N1;
F2 = (4*k+2)*r.^k.*N1+4*r.^(k+2).*N2;
F3 = 12*(k+1)*r.^(k+1).*N2+8*r.^(k+3).*N3;
if k >= 1
    F1 = F1+k*r.^(k-1).*N0;
    F3 = F3+6*k^2*r.^(k-1).*N1;
end
if k >= 2
    F2 = F2+k*(k-1)*r.^(k-2).*N0;
end
if k >= 3
    F3 = F3+k*(k-1)*(k-2)*r.^(k-3).*N0;
end
end


function [F0,F1,F2,F3] = radial_power_basis(N0,N1,N2,N3,r,k)
% f(r)=r^k*N(r), with derivatives through third order.
r = r(:);

F0 = r.^k.*N0;
F1 = r.^k.*N1;
F2 = r.^k.*N2;
F3 = r.^k.*N3;

if k >= 1
    F1 = F1 + k*r.^(k-1).*N0;
    F2 = F2 + 2*k*r.^(k-1).*N1;
    F3 = F3 + 3*k*r.^(k-1).*N2;
end
if k >= 2
    F2 = F2 + k*(k-1)*r.^(k-2).*N0;
    F3 = F3 + 3*k*(k-1)*r.^(k-2).*N1;
end
if k >= 3
    F3 = F3 + k*(k-1)*(k-2)*r.^(k-3).*N0;
end
end


function [N0,N1,N2,N3] = shifted_legendre_all(z,ncols)
% Shifted Legendre polynomials P_n(2*z-1) and z derivatives.
z = z(:);
x = 2*z-1;
Nq = numel(z);
N0 = zeros(Nq,ncols);
N1 = zeros(Nq,ncols);
N2 = zeros(Nq,ncols);
N3 = zeros(Nq,ncols);
N0(:,1) = 1;
if ncols == 1, return; end
N0(:,2) = x;
N1(:,2) = 2;

for n = 1:(ncols-2)
    a = 2*n+1;
    b = n;
    c = n+1;
    N0(:,n+2) = (a*x.*N0(:,n+1)-b*N0(:,n))/c;
    N1(:,n+2) = (a*(2*N0(:,n+1)+x.*N1(:,n+1))-b*N1(:,n))/c;
    N2(:,n+2) = (a*(4*N1(:,n+1)+x.*N2(:,n+1))-b*N2(:,n))/c;
    N3(:,n+2) = (a*(6*N2(:,n+1)+x.*N3(:,n+1))-b*N3(:,n))/c;
end
end


function [F0,F1,F2,F3] = pure_power_column(r,k)
r = r(:);
F0 = r.^k;
F1 = k*r.^(k-1);
F2 = zeros(size(r));
F3 = zeros(size(r));
if k >= 2
    F2 = k*(k-1)*r.^(k-2);
end
if k >= 3
    F3 = k*(k-1)*(k-2)*r.^(k-3);
end
end
