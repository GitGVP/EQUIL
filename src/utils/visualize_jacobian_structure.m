function visualize_jacobian_structure(J, dof_count, Nb, Ns, jacTotal, M_extended, P_templates)
    % Enhanced visualization with per-derivative contributions

    % Create variable names
    n_vars = 3 + Nb + Ns;
    block_size = dof_count;
    block_size_S = dof_count - 1;
    var_names = cell(1,n_vars);
    var_names{1} = 't2';
    var_names{2} = 'delta';
    var_names{3} = 'P';
    for k = 1:Nb
        var_names{3+k} = sprintf('B%d', k-1);
    end
    for k = 1:Ns
        var_names{3+Nb+k} = sprintf('S%d', k+1);
    end

    % Create figure with more subplots
    figure('Position', [50, 50, 2000, 1200]);

    %% Row 1: Original plots
    % Sparsity pattern
    subplot(2, 4, 1);
    spy(J);
    title(sprintf('Full Jacobian Sparsity\n(nnz=%d, %.1f%%)', ...
        nnz(J), 100*nnz(J)/numel(J)));
    xlabel('Columns (dofs)'); ylabel('Rows (equations)');
    grid on;

    % Singular values
    subplot(2,4,2);
    s = svd(full(J));
    semilogy(s, 'o-', 'LineWidth', 1.5, 'MarkerSize', 4);
    title(sprintf('Singular Values\n(cond=%.2e)', cond(J)));
    xlabel('Index'); ylabel('\sigma_i');
    grid on;

    % Norm of each block
    subplot(2,4,3);
    block_norms = zeros(n_vars, n_vars);
    start_row = 1;
    for i = 1:n_vars
        var_size = block_size;
        if i > 3+Nb
            var_size = block_size_S;
        end
        start_col = 1;
        for j = 1:n_vars
            dest_size = block_size;
            if j > 3+Nb
                dest_size = block_size_S;
            end
            block = J(start_row:start_row+var_size-1, start_col:start_col+dest_size-1);
            block_norms(i,j) = norm(block, 'fro');
            start_col = start_col + dest_size;
        end
        start_row = start_row + var_size;
    end
    imagesc(log10(block_norms + 1e-16));
    colorbar;
    title('log10(Frobenius norm)\nof Blocks');
    set(gca, 'XTick', 1:n_vars, 'XTickLabel', var_names, ...
        'XTickLabelRotation', 45, 'FontSize', 7);
    set(gca, 'YTick', 1:n_vars, 'YTickLabel', var_names, 'FontSize', 7);
    axis tight;

    % Per-derivative contribution comparison
    subplot(2,4,4);
    totalDofs = size(J,1);
    nCols = size(jacTotal, 4) * (n_vars)^2;
    nnz_per_d = zeros(3,1);
    norm_per_d = zeros(3,1);

    for d = 1:3
        jac_d = squeeze(jacTotal(d,:,:,:));
        jacTotal_flat = reshape(permute(jac_d, [3, 1, 2]), [], 1);
        vP_scaled = P_templates{d}.v_template .* jacTotal_flat(P_templates{d}.i);
        P_scaled = sparse(P_templates{d}.i, P_templates{d}.j, vP_scaled, nCols, totalDofs);
        J_d = M_extended * P_scaled;
        nnz_per_d(d) = nnz(J_d);
        norm_per_d(d) = norm(J_d, 'fro');
    end

    yyaxis left
    bar(1:3, nnz_per_d);
    ylabel('Number of nonzeros');
    yyaxis right
    plot(1:3, norm_per_d, 'ro-', 'LineWidth', 2, 'MarkerSize', 8);
    ylabel('Frobenius norm');
    xlabel('Derivative order d');
    title(sprintf('Contribution per d\nTotal nnz=%d', nnz(J)));
    set(gca, 'XTick', 1:3, 'XTickLabel', {'d=1', 'd=2', 'd=3'});
    grid on;

    %% Row 2: Per-derivative sparsity patterns
    for d = 1:3
        subplot(2, 4, 4+d);

        % Reconstruct J contribution from d
        jac_d = squeeze(jacTotal(d,:,:,:));
        jacTotal_flat = reshape(permute(jac_d, [3, 1, 2]), [], 1);
        vP_scaled = P_templates{d}.v_template .* jacTotal_flat(P_templates{d}.i);
        P_scaled = sparse(P_templates{d}.i, P_templates{d}.j, vP_scaled, nCols, totalDofs);
        J_d = M_extended * P_scaled;

        spy(J_d);
        title(sprintf('d=%d Sparsity\n(nnz=%d, norm=%.2e)', ...
            d, nnz(J_d), norm(J_d, 'fro')));
        xlabel('Columns'); ylabel('Rows');
    end

    %% Last subplot: jacTotal statistics across radial points
    subplot(2, 4, 8);
    Nq = size(jacTotal, 4);

    % Compute max absolute value for each (d, equation, profile) across radial points
    max_vals = squeeze(max(abs(jacTotal), [], 4));  % 3 x n_vars x n_vars

    % Show heatmap of max values averaged over derivatives
    avg_max = squeeze(mean(max_vals, 1));  % n_vars x n_vars
    imagesc(log10(avg_max + 1e-16));
    colorbar;
    title(sprintf('log10(max|jacTotal|)\naveraged over d, max over q'));
    set(gca, 'XTick', 1:n_vars, 'XTickLabel', var_names, ...
        'XTickLabelRotation', 45, 'FontSize', 7);
    set(gca, 'YTick', 1:n_vars, 'YTickLabel', var_names, 'FontSize', 7);
    xlabel('Profile (gamma)');
    ylabel('Equation (alpha)');
    axis tight;

    sgtitle(sprintf('Jacobian Analysis (dof=%d, Nb=%d, Ns=%d, Nq=%d)', ...
        dof_count, Nb, Ns, Nq), 'FontSize', 14, 'FontWeight', 'bold');
end