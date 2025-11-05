function plot_profile_panel(r, ana_vec, num_mat, cmap, col_inds, title_str, ylabel_str)
    % ana_vec : 1 x Nq (analytical, same for all eps)
    % num_mat : nsim x Nq (numerical profiles for each eps)
    hold on;
    plot(r, ana_vec, 'k--', 'LineWidth', 3);
    ncurves = size(num_mat,1);
    for ii = 1:ncurves
        color_i = cmap(col_inds(ii), :);
        plot(r, num_mat(ii,:), '-', ...
                'Color',color_i);
    end
    grid on;
    xlabel('$\hat r$','Interpreter','latex','FontSize',12);
    ylabel(ylabel_str,'Interpreter','latex','FontSize',12);
    title(title_str, 'Interpreter','latex','FontSize',12);
    hold off;
end