function [dbetapardB0, dbetapardr, dbetapardB, dbetapardR, d2betapardB2, ...
    d2betapardrdB, d2betapardRdB, d2betapardrdR, d2betapardR2, d3betapardrdB2, ...
    d3betapardB3, d3betapardBdR2, d3betapardB2dR, d3betapardrdRdB, ...
    d3betapardrdR2, betapar, betaperp] = biMaxwellian(kinetic_profiles,r,RR,BB,varargin)
% BIMAXWELLIAN Bi-Maxwellian pressure closure and analytical derivatives.
%
% The pressure is written in both branches as
%
%   betapar = beta(r) * x * f(y),
%   x = B/Bc(r),
%   y = Ah(r) * (1 - x).
%
% For B >= Bc,
%
%   f(y) = 1/(1-y),
%
% which is algebraically identical to the original rational expression.
% For B < Bc,
%
%   f(y) = (1 + y - 2*y^(5/2))/(1-y^2).
%
% All B and r derivatives below are derivatives of the same branch used
% for betapar. The radial derivative is the partial derivative at fixed
% (R,B). There is no explicit R dependence in this closure.
%
% The low-B branch is C^2 across B = Bc, while third derivatives involving
% sufficiently many B derivatives are singular as B -> Bc^- because of
% the y^(5/2) term. This is physical/mathematical behaviour of the chosen
% regularization and is not artificially removed here.

    beta  = kinetic_profiles.beta(r)  * ones(1,size(RR,2));
    betap = kinetic_profiles.betap(r) * ones(1,size(RR,2));
    Ah    = kinetic_profiles.Ah(r)    * ones(1,size(RR,2));
    Ahp   = kinetic_profiles.Ahp(r)   * ones(1,size(RR,2));
    Bc    = kinetic_profiles.Bc(r)    * ones(1,size(RR,2));
    Bcp   = kinetic_profiles.Bcp(r)   * ones(1,size(RR,2));

    % Compact variables.
    x = BB ./ Bc;
    y = Ah .* (1 - x);

    % Radial derivatives at fixed B.
    c  = Bcp ./ Bc;                           % Bc_r / Bc
    yr = Ahp .* (1 - x) + Ah .* x .* c;      % (partial y / partial r)_B

    % ------------------------------------------------------------------
    % f(y) and its y-derivatives.
    %
    % Start with the B >= Bc branch:
    %   f = 1/(1-y).
    % ------------------------------------------------------------------
    one_minus_y = 1 - y;
    f  = 1 ./ one_minus_y;
    f1 = 1 ./ one_minus_y.^2;
    f2 = 2 ./ one_minus_y.^3;

    if nargout > 7
        f3 = 6 ./ one_minus_y.^4;
    end

    % ------------------------------------------------------------------
    % Replace by the actual low-B branch wherever B < Bc.
    %
    %   f = n/d,
    %   n = 1 + y - 2 y^(5/2),
    %   d = 1 - y^2.
    %
    % Derivatives are evaluated recursively from d*f = n:
    %   d f'   + d' f                         = n'
    %   d f''  + 2 d' f' + d'' f            = n''
    %   d f''' + 3 d' f'' + 3 d'' f'        = n'''
    % ------------------------------------------------------------------
    idx = BB < Bc;

    if any(y(idx) < 0)
        error('biMaxwellian:InvalidLowBBranch', ...
            ['The low-B branch requires Ah*(1-B/Bc) >= 0 for the ', ...
             'fractional power y^(5/2) to remain real.']);
    end

    % At y = 0 the low-B branch and its derivatives required through
    % second order have the same limiting values as the high-B branch.
    % Only overwrite strictly positive y, which also avoids 0*Inf in the
    % third-derivative chain rule when Ah = 0.
    idx_low = idx & (y > 0);

    if any(idx_low(:))
        yl = y(idx_low);

        n   = 1 + yl - 2 .* yl.^2.5;
        n1  = 1 - 5 .* yl.^1.5;
        n2  = -(15/2) .* yl.^0.5;

        d   = 1 - yl.^2;
        d1  = -2 .* yl;
        d2  = -2 .* ones(size(yl));

        fl  = n ./ d;
        f1l = (n1 - d1 .* fl) ./ d;
        f2l = (n2 - 2 .* d1 .* f1l - d2 .* fl) ./ d;

        f(idx_low)  = fl;
        f1(idx_low) = f1l;
        f2(idx_low) = f2l;

        if nargout > 7
            % Singular as y -> 0+; this is the genuine third derivative
            % of the y^(5/2) low-B correction.
            n3  = -(15/4) .* yl.^(-0.5);
            f3l = (n3 - 3 .* d1 .* f2l - 3 .* d2 .* f1l) ./ d;
            f3(idx_low) = f3l;
        end
    end

    % ------------------------------------------------------------------
    % Pressure derivatives.
    %
    % P = beta * x * f(y)
    %
    % with
    %   x_B = 1/Bc,        y_B = -Ah/Bc,
    %   x_r = -x*c,        y_r = yr,
    %
    % where r derivatives are taken at fixed B.
    % ------------------------------------------------------------------

    % P_r
    dbetapardr = x .* ( ...
        (betap - beta .* c) .* f + beta .* f1 .* yr);

    % P_B
    H1 = f - Ah .* x .* f1;
    dbetapardB = (beta ./ Bc) .* H1;

    % No explicit R dependence.
    dbetapardR = zeros(size(RR));

    % P_BB
    H2 = -2 .* Ah .* f1 + Ah.^2 .* x .* f2;
    d2betapardB2 = (beta ./ Bc.^2) .* H2;

    % P_rB
    H1r = ...
        f1 .* yr ...
        - Ahp .* x .* f1 ...
        + Ah .* x .* c .* f1 ...
        - Ah .* x .* f2 .* yr;

    d2betapardrdB = ...
        ((betap - beta .* c) ./ Bc) .* H1 ...
        + (beta ./ Bc) .* H1r;

    % P_RB = 0.
    d2betapardRdB = zeros(size(RR));

    % Axis/reference value used by the equilibrium code.
    dbetapardB0 = dbetapardB(1,1);

    if nargin > 4
        betapar = beta .* x .* f;
        betaperp = betapar-BB.*dbetapardB;
        dbetapardB0 = variational_pressure_output( ...
            betapar,dbetapardr,dbetapardR,dbetapardB, ...
            d2betapardB2,d2betapardrdB,d2betapardRdB, ...
            zeros(size(RR)),zeros(size(RR)),betaperp);
        return
    end

    if nargout > 7
        % No explicit R dependence.
        d2betapardrdR = zeros(size(RR));
        d2betapardR2  = zeros(size(RR));

        % P_rBB
        H2r = ...
            -2 .* Ahp .* f1 ...
            -2 .* Ah .* f2 .* yr ...
            +2 .* Ah .* Ahp .* x .* f2 ...
            -Ah.^2 .* x .* c .* f2 ...
            +Ah.^2 .* x .* f3 .* yr;

        d3betapardrdB2 = ...
            ((betap - 2 .* beta .* c) ./ Bc.^2) .* H2 ...
            + (beta ./ Bc.^2) .* H2r;

        % P_BBB
        d3betapardB3 = (beta ./ Bc.^3) .* ( ...
            3 .* Ah.^2 .* f2 - Ah.^3 .* x .* f3);

        % All derivatives containing explicit R derivatives vanish.
        d3betapardBdR2   = zeros(size(RR));
        d3betapardB2dR   = zeros(size(RR));
        d3betapardrdRdB  = zeros(size(RR));
        d3betapardrdR2   = zeros(size(RR));

        if nargout > 14
            betapar = beta .* x .* f;

            % P_perp = P_parallel - B * (partial P_parallel / partial B).
            % This must use the branch-consistent PB computed above.
            betaperp = betapar - BB .* dbetapardB;
        end
    end
end
