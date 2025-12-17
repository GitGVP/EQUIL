function check_projection_matrices(P0, P1, P2, r, basis_name)
    % Check basic properties of discretization matrices
    fprintf('\n=== Discretization Validation (%s) ===\n', basis_name);
    fprintf('P0 size: %dx%d\n', size(P0));
    fprintf('P1 size: %dx%d\n', size(P1));
    fprintf('P2 size: %dx%d\n', size(P2));
    
    % Test derivative consistency on polynomials
    dof = size(P0, 2);
    test_coeffs = randn(dof, 3);
    
    % Create figure for visualization
    figure('Position', [100, 100, 1200, 900]);
    
    for i = 1:3
        f_c = test_coeffs(:, i);
        f = P0 * f_c;
        f_prime = P1 * f_c;
        f_dprime = P2 * f_c;
        
        % Numerical derivatives for comparison
        f_prime_num = gradient(f, r);
        f_dprime_num = gradient(f_prime_num, r);
        
        fprintf('\nTest function %d:\n', i);
        rel_error1 = norm(f_prime - f_prime_num) / norm(f_prime_num);
        rel_error2 = norm(f_dprime - f_dprime_num) / norm(f_dprime_num);
        fprintf('  ||P1*f - numerical f''||/||f''||: %.2e\n', rel_error1);
        fprintf('  ||P2*f - numerical f"||/||f"||: %.2e\n', rel_error2);
        
        % Plot f and its derivative
        subplot(3, 3, (i-1)*3 + 1);
        plot(r, f, 'b-', 'LineWidth', 2);
        title(sprintf('Test function %d: f(r)', i));
        xlabel('r'); ylabel('f(r)');
        grid on;
        
        % Plot first derivative comparison
        subplot(3, 3, (i-1)*3 + 2);
        plot(r, f_prime, 'r-', 'LineWidth', 2, 'DisplayName', 'P1*f');
        hold on;
        plot(r, f_prime_num, 'k--', 'LineWidth', 1.5, 'DisplayName', 'Numerical f''');
        plot(r, f_prime - f_prime_num, 'g:', 'LineWidth', 1, 'DisplayName', 'Difference');
        hold off;
        title(sprintf('First derivative (error: %.1e)', rel_error1));
        xlabel('r'); ylabel('f''(r)');
        legend('Location', 'best');
        grid on;
        
        % Plot second derivative comparison
        subplot(3, 3, (i-1)*3 + 3);
        plot(r, f_dprime, 'r-', 'LineWidth', 2, 'DisplayName', 'P2*f');
        hold on;
        plot(r, f_dprime_num, 'k--', 'LineWidth', 1.5, 'DisplayName', 'Numerical f"');
        plot(r, f_dprime - f_dprime_num, 'g:', 'LineWidth', 1, 'DisplayName', 'Difference');
        hold off;
        title(sprintf('Second derivative (error: %.1e)', rel_error2));
        xlabel('r'); ylabel('f"(r)');
        legend('Location', 'best');
        grid on;
    end
    
    % Add overall title
    sgtitle(sprintf('Discretization Validation: %s Basis', basis_name));
    
    % Check null spaces
    fprintf('\nNull space checks:\n');
    fprintf('  rank(P0): %d (should be %d)\n', rank(full(P0)), min(size(P0)));
    fprintf('  cond(P0): %.2e\n', cond(full(P0)));
    
    % Additional diagnostics
    fprintf('\nAdditional diagnostics:\n');
    fprintf('  ||P1 - gradient(P0)||/||P1||: %.2e\n', ...
        norm(P1 - gradient_matrix(r) * P0, 'fro') / norm(P1, 'fro'));
    
    % Visualize the matrices
    figure('Position', [100, 100, 1400, 400]);
    
    subplot(1, 3, 1);
    spy(P0);
    title(sprintf('P0 Sparsity (nnz=%d, %.1f%%)', nnz(P0), 100*nnz(P0)/numel(P0)));
    xlabel('Coefficients'); ylabel('Radial Points');
    
    subplot(1, 3, 2);
    spy(P1);
    title(sprintf('P1 Sparsity (nnz=%d, %.1f%%)', nnz(P1), 100*nnz(P1)/numel(P1)));
    xlabel('Coefficients'); ylabel('Radial Points');
    
    subplot(1, 3, 3);
    spy(P2);
    title(sprintf('P2 Sparsity (nnz=%d, %.1f%%)', nnz(P2), 100*nnz(P2)/numel(P2)));
    xlabel('Coefficients'); ylabel('Radial Points');
    
    sgtitle(sprintf('%s Basis Matrices', basis_name));
end

function G = gradient_matrix(r)
    % Create a simple finite difference matrix for comparison
    n = length(r);
    G = zeros(n, n);
    dr = r(2) - r(1);
    
    for i = 2:n-1
        G(i, i-1) = -1/(2*dr);
        G(i, i+1) = 1/(2*dr);
    end
    G(1, 1) = -1/dr; G(1, 2) = 1/dr;
    G(n, n-1) = -1/dr; G(n, n) = 1/dr;
end