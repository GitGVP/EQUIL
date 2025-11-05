% testing in a script
clear all;

% used in G_JFNK_FULL
nq = 2; % Number of quadrature points

syms xi_sym;
phi_sym = [2*(xi_sym-0.5)*(xi_sym-1);
          4*xi_sym*(1-xi_sym);
          2*xi_sym*(xi_sym-0.5)];
dphi_dxi_sym = diff(phi_sym, xi_sym);
d2phi_dxi2_sym = diff(dphi_dxi_sym, xi_sym);

dopp = false;
%ms = linspace(20,60,11);
%ms = [60];
ms = linspace(20,50,7);
L2s = []; H1s = []; r0s = [];
!rm last_solution.mat
for m = ms
    H1_elem = []; L2_elem = [];
    h = 1/(2*m);
    nq = 2; % Number of quadrature points
    G_JFNK_full
    !rm last_solution.mat
    if res_norm > 1e-2
        L2s = [L2s, NaN]; H1s = [H1s, NaN];
        continue;
    end
    x_full = [0; x(1:2*m); ... % t2
        0; x(2*m+1:4*m);... % delta
        0; x(4*m+1:6*m);... % P
        0; x(6*m+1:8*m);... % B1
        0; x(8*m+1:end);S2bc;... % S2
        ];        


    %gauss_pts = [(1 - sqrt(3/5))/2, 0.5, (1 + sqrt(3/5))/2];
    %gauss_wts = [5/18, 4/9, 5/18];
    nq = 10; % Number of quadrature points
    [gauss_pts, gauss_wts] = lgwt(nq, 0, 1);


    phi = zeros(3,3);
    dphi_dxi = zeros(3,3);
    d2phi_dxi2 = zeros(3,3);
    for i = 1:3
        for q = 1:nq
            xi_val = gauss_pts(q);
            phi(i,q) = double(subs(phi_sym(i), xi_sym, xi_val));
            dphi_dxi(i,q) = double(subs(dphi_dxi_sym(i), xi_sym, xi_val));
            d2phi_dxi2(i,q) = double(subs(d2phi_dxi2_sym(i), xi_sym, xi_val));
        end
    end
    dphi_dr = dphi_dxi / (2*h);
    d2phi_dr2 = d2phi_dxi2 / (2*h)^2;

    % --- initialize accumulators ---
    L2_err_sq   = 0.0;
    H1_sem_sq   = 0.0;

    % Loop over elements
    for e = 1:m
        % local DOF indices into U
        nodes = [2*e-1, 2*e, 2*e+1];     % left vertex, mid, right vertex
        r0 = r_nodes(nodes(1));          % left node coordinate
        r0s = [r0s, r0];
        % loop over quadrature points (xi in [0,1])
        for q = 1:nq
            q_global = (e-1)*nq + q;     % if you have r_q precomputed
            rq = r0 + 2*h*gauss_pts(q); % physical quadrature point on this element
            uh_q = 0; duh_dr_q = 0;
            for i = 1:3 
                node_global = nodes(i);
                if node_global == 1
                    continue  
                end
                uh_q = uh_q + x_full(nodes(i))*phi(i,q);
                duh_dr_q = duh_dr_q + x_full(nodes(i))*dphi_dr(i,q);
            end

            % exact solution at rq
            uex = interp1(r_fine, t2_fine, rq, 'spline');%t2_ana(q_global);
            duex = interp1(r_fine, t2p_fine, rq, 'spline');%t2p_ana(q_global);
            %duex = results.t2p(q_global);
            % pointwise differences
            diff   = uex - uh_q;
            diffp  = duex - duh_dr_q;
            % Jacobian for xi in [0,1] -> r in [r0, r0+2h] is dx/dxi = 2*h
            J = 2*h;

            % accumulate 
            L2_err_sq = L2_err_sq + J * gauss_wts(q) * (diff^2);
            H1_sem_sq = H1_sem_sq + J * gauss_wts(q) * (diffp^2);
        end
        L2_elem = [L2_elem, sqrt(L2_err_sq)];
        H1_elem = [H1_elem, sqrt(H1_sem_sq)];
    end

    % final errors
    L2_err   = sqrt(L2_err_sq);
    H1_semi  = sqrt(H1_sem_sq);

    fprintf('L2 error = %.6e\n', L2_err);
    fprintf('H1 seminorm error = %.6e\n', H1_semi);
    L2s = [L2s, L2_err]; H1s = [H1s, H1_semi];
end


hvals = 1./(2*ms);

figure('WindowStyle','docked');
tiledlayout(2,1, 'TileSpacing','compact','Padding','compact');

nexttile;
loglog(hvals, L2s, 'o-', 'LineWidth', 1.5); 
hold on;
p = polyfit(log(hvals), log(L2s), 1);
hfit = linspace(min(hvals), max(hvals), 100);
slope = p(1);
L2fit = exp(polyval(p, log(hfit)));
loglog(hfit, L2fit, '--', 'LineWidth', 1.5);
legend('L2 error', sprintf('Fit slope = %.2f', slope), 'Location','best');
ylabel('L2 error');
grid on;

nexttile;
loglog(hvals, H1s, 'o-', 'LineWidth', 1.5); 
hold on;
p = polyfit(log(hvals), log(H1s), 1);
slope = p(1);
hfit = linspace(min(hvals), max(hvals), 100);
L2fit = exp(polyval(p, log(hfit)));
loglog(hfit, L2fit, '--', 'LineWidth', 1.5);
legend('H1 error', sprintf('Fit slope = %.2f', slope), 'Location','best');
xlabel('h'); ylabel('H1 error');
grid on;


% %% Test interpolation error 
% 
% % --- pick a smooth exact solution ---
% eps_val = 0.1;
% u_exact = @(r) -r.^2 + 3*r.^4 * eps_val^2 / 32 + 113 * r.^6 * eps_val ^ 4 / 2048;
% du_exact = @(r) -2 *r + 3*4*r.^3 * eps_val^2 / 32 + 113 * 6 * r.^5 * eps_val ^ 4 / 2048;
% 
% % pick several meshes for the test
% mvals = linspace(10,100,10); 
% L2_test = zeros(size(mvals));
% H1_test = zeros(size(mvals));
% 
% for k = 1:numel(mvals)
%     mloc = mvals(k);
%     hloc = 1/(2*mloc);
%     r_nodes_loc = linspace(0,1,2*mloc+1);
% 
%     % build nodal DOFs by sampling exact solution at vertices and midpoints
%     Uloc = zeros(2*mloc+1,1);
%     for ii = 1:(2*mloc+1)
%         Uloc(ii) = u_exact(r_nodes_loc(ii));
%     end
%     % build phi, dphi_dr etc as you already do for each h (reuse your code)
%     % For brevity assume phi, dphi_dr and gauss_pts, gauss_wts already match hloc:
%     % (recompute them if needed using your symbolic setup with hloc)
%     % Now call your error loop replacing U by Uloc (use the code you already have)
%     % Suppose you have a function compute_errors(U, m, h, r_nodes, phi, dphi_dr, ...)
%     % [L2_test(k), H1_test(k)] = compute_errors(Uloc, mloc, ...);
% 
%     % --- If you don't have the function, simply inline your loop here using Uloc ---
%     % (Use the loop from previous messages to compute L2 and H1 with Uloc)
%         % --- initialize accumulators ---
%     phi = zeros(3,3);
%     dphi_dxi = zeros(3,3);
%     d2phi_dxi2 = zeros(3,3);
%     for i = 1:3
%         for q = 1:3
%             xi_val = gauss_pts(q);
%             phi(i,q) = double(subs(phi_sym(i), xi_sym, xi_val));
%             dphi_dxi(i,q) = double(subs(dphi_dxi_sym(i), xi_sym, xi_val));
%             d2phi_dxi2(i,q) = double(subs(d2phi_dxi2_sym(i), xi_sym, xi_val));
%         end
%     end
%     dphi_dr = dphi_dxi / (2*hloc);
%     d2phi_dr2 = d2phi_dxi2 / (2*hloc)^2;
%     L2_err_sq   = 0.0;
%     H1_sem_sq   = 0.0;
% 
%     % Loop over elements
%     for e = 1:mloc
%         % local DOF indices into U
%         nodes = [2*e-1, 2*e, 2*e+1];     % left vertex, mid, right vertex
%         r0 = r_nodes_loc(nodes(1));          % left node coordinate
%         % local DOF values
%         u_l = Uloc(nodes(1));
%         u_m = Uloc(nodes(2));
%         u_r = Uloc(nodes(3));
% 
%         % loop over quadrature points (xi in [0,1])
%         for q = 1:3
%             q_global = (e-1)*3 + q;     % if you have r_q precomputed
%             rq = r0 + 2*hloc*gauss_pts(q); % physical quadrature point on this element
% 
%             % finite-element approximation at rq (using local DOFs)
%             uh_q = u_l*phi(1,q) + u_m*phi(2,q) + u_r*phi(3,q);
% 
%             % derivative w.r.t r at rq (use your dphi_dr)
%             duh_dr_q = u_l*dphi_dr(1,q) + u_m*dphi_dr(2,q) + u_r*dphi_dr(3,q);
% 
%             % exact solution at rq
%             uex = u_exact(rq);
%             duex = du_exact(rq);
% 
%             % pointwise differences
%             diff   = uex - uh_q;
%             diffp  = duex - duh_dr_q;
% 
%             % Jacobian for xi in [0,1] -> r in [r0, r0+2h] is dx/dxi = 2*h
%             J = 2*hloc;
% 
%             % accumulate (using your gauss weights on [0,1])
%             L2_err_sq = L2_err_sq + J * gauss_wts(q) * (diff^2);
%             H1_sem_sq = H1_sem_sq + J * gauss_wts(q) * (diffp^2);
%         end
%     end
%     L2_test(k) = sqrt(L2_err_sq);
%     H1_test(k) = sqrt(H1_sem_sq);
% end
% 
% disp([mvals; L2_test; H1_test]');
% 
% L2s = L2_test;
% H1s = H1_test;
% hvals = 1./(2*mvals);
