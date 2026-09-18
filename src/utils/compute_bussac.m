function [norm_growth,bpbar,sigma,tau1,b,c] = compute_bussac(L,LX,LY)
%COMPUTE_BUSSAC Leading-order anisotropic Bussac internal-kink estimate.
%   [NORM_GROWTH,BPBAR,SIGMA,TAU1,B,C] = COMPUTE_BUSSAC(L,LX,LY)
%   returns NORM_GROWTH = (gamma/omega_A)/epsilon^2.  L and LX must be
%   variational-EQUIL inputs, and LY must include the analytical output
%   produced with do_ana=true.  Only q, the pressure map, and the
%   analytical radial cutoff are used; the solved EQUIL state is not.

    if ~isfield(LY,'r_fine') || isempty(LY.r_fine)
        error('compute_bussac:MissingAnalyticalOutput', ...
            'LY must contain equil_ana output (run with do_ana=true).');
    end
    if ~isfield(LX,'qfun') || ~isfield(LX,'qpfun')
        error('compute_bussac:MissingQProfile', ...
            'LX.qfun and LX.qpfun are required.');
    end

    q = LX.qfun;
    qp = LX.qpfun;
    rs = fzero(@(r) q(r)-1,[0,1]);
    shear = rs*qp(rs)/q(rs);
    if shear == 0
        error('compute_bussac:ZeroShear','The q=1 shear must be nonzero.');
    end

    % A modest Gauss--Legendre rule is ample for the smooth input maps.
    [r,w] = lgwt(64,0,rs);
    [Pr,PrB,PrBB] = pressure_derivatives(L,LX,r);
    mean_Pr = Pr-0.5*PrB;
    bpbar = -2*sum(w.*r.^2.*mean_Pr)/rs^4;
    tau1_bpbar = -sum(w.*r.^2.*(PrB-PrBB))/(2*rs^4);
    sigma = sum(w.*r.^3.*(1./q(r).^2-1))/rs^4;

    if bpbar == 0
        tau1 = NaN;
    else
        tau1 = tau1_bpbar/bpbar;
    end

    % Regular homogeneous upper-sideband solutions inside q=1 and outside
    % it.  If q=2 lies in the plasma, regularity there is Neumann; if not,
    % the outer solution is set to zero at the plasma edge.
    cutoff = LY.r_fine(1);
    if q(1) >= 2
        r_outer = fzero(@(r) q(r)-2,[rs,1]);
        outer_endpoint = r_outer-cutoff;
        outer_neumann = true;
    else
        outer_endpoint = 1;
        outer_neumann = false;
    end
    if outer_endpoint <= rs
        error('compute_bussac:SidebandInterval', ...
            'The outer sideband interval is too short for the cutoff.');
    end
    [b,c] = sideband_slopes( ...
        q,qp,rs,cutoff,outer_endpoint,outer_neumann);

    beta_sigma = bpbar+sigma;
    W = tau1_bpbar ...
        +(32*(b-c)*sigma+9*(b-1)*(1-c))/(64*(b-c)) ...
        -3*(b-1)*(c+3)*beta_sigma/(8*(b-c)) ...
        -(b+3)*(c+3)*beta_sigma^2/(4*(b-c));
    norm_growth = -pi*rs^2*W/(sqrt(3)*shear);
end

function [b,c] = sideband_slopes( ...
        q,qp,rs,cutoff,outer_endpoint,outer_neumann)
% Integrate arbitrarily normalized regular solutions; their scale cancels.
    odefun = @(r,Y) sideband_ode(r,Y,q,qp);
    options = odeset('RelTol',1e-8,'AbsTol',1e-10);

    [~,inner] = ode113(odefun,[cutoff,rs],[cutoff;1],options);
    b = rs*inner(end,2)/inner(end,1);

    if outer_neumann
        outer_initial = [1;0];
    else
        outer_initial = [0;1];
    end
    [~,outer] = ode113( ...
        odefun,[outer_endpoint,rs],outer_initial,options);
    c = rs*outer(end,2)/outer(end,1);
end

function dY = sideband_ode(r,Y,q,qp)
    qr = q(r);
    qpr = qp(r);
    d2 = (-1+2/qr)^2*r^2;
    d1 = -((2-qr)*r*(-6*qr+3*qr^2+4*qpr*r))/qr^3;
    d0 = -3*(-1+2/qr)^2;
    dY = [Y(2);-(d1*Y(2)+d0*Y(1))/d2];
end

function [Pr,PrB,PrBB] = pressure_derivatives(L,LX,r)
% Return derivatives of P_parallel, where Pi_parallel=epsilon^2 P_parallel.
    nr = numel(r);
    if isfield(LX,'pressure_map') && ~isempty(LX.pressure_map)
        pressure = LX.pressure_map(r,ones(nr,1),ones(nr,1));
        required = {'Pr','PrB','PrBB'};
        for k = 1:numel(required)
            if ~isfield(pressure,required{k})
                error('compute_bussac:PressureDerivative', ...
                    'LX.pressure_map must return field %s.',required{k});
            end
        end
        scale = LX.eps_val^2;
        Pr = expand_column(pressure.Pr,nr)/scale;
        PrB = expand_column(pressure.PrB,nr)/scale;
        PrBB = expand_column(pressure.PrBB,nr)/scale;
        return
    end

    one = ones(nr,numel(L.omega));
    [~,Pr,~,~,~,PrB,~,~,~,PrBB] = ...
        L.P.equation_of_state(LX.kinetic_profiles,r,one,one);
    Pr = mean(Pr,2);
    PrB = mean(PrB,2);
    PrBB = mean(PrBB,2);
end

function value = expand_column(value,nr)
    if isscalar(value)
        value = value+zeros(nr,1);
    elseif numel(value) == nr
        value = value(:);
    else
        error('compute_bussac:PressureSize', ...
            'A pressure derivative has %d entries; expected %d.', ...
            numel(value),nr);
    end
end
