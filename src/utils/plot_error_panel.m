
function plot_error_panel(err_all, values_eps, outerloop_pvals, cmap, col_inds, marker_list, title_str, outerloop_legend_name, skip_count)    % err_all : n_om x n_eps
    n_om = size(err_all,1);
    nsim = size(err_all,2);

    % compute per-om fits
    pcell = cell(n_om,1);
    for k = 1:n_om
        x = log10(values_eps(:));
        y = log10((err_all(k,:)).');
        valid = isfinite(x) & isfinite(y) & ((err_all(k,:)).' > 0);
        if skip_count > 0
            nskip = min(skip_count, numel(valid));
            % skip the first nskip entries in the current ordering of values_eps
            valid(1:nskip) = false;
        end
        if any(valid)
            pcell{k} = polyfit(x(valid), y(valid), 1);
        else
            pcell{k} = [NaN NaN];
        end
    end

    % plot
    hold on;
    set(gca,'XScale','log','YScale','log');
    grid on;


    for k = 1:n_om
        x_fit = linspace(min(log10(values_eps)), max(log10(values_eps)), 200);
        y_fit = polyval(pcell{k}, x_fit);
        loglog(10.^x_fit, 10.^y_fit, 'k-', 'LineWidth', 1.0);
        marker_k = marker_list{mod(k-1, numel(marker_list))+1};
        for i = 1:numel(values_eps)
            color_i = cmap(col_inds(i), :);
            v = err_all(k,i);
            if ~isfinite(v) || v<=0
                continue
            end
            loglog(values_eps(i), v, marker_k, 'MarkerFaceColor', color_i, ...
                   'MarkerEdgeColor', 'k', 'MarkerSize', 6);
        end
    end

    xlabel('$\epsilon$','Interpreter','latex','FontSize',12);
    ylabel(title_str,'Interpreter','latex','FontSize',12);

    % legend with same markers and slope
    h = gobjects(n_om,1);
    labels = cell(n_om,1);
    for k = 1:n_om
        marker_k = marker_list{mod(k-1, numel(marker_list))+1};
        h(k) = plot(nan, nan, marker_k, 'MarkerFaceColor', 'k', ...
                    'MarkerEdgeColor', 'k', 'MarkerSize', 6, 'LineStyle', 'none');
        slope_k = pcell{k}(1);
        labels{k} = sprintf('$%s=%g,\\ \\mathrm{slope}=%.2f$', outerloop_legend_name, outerloop_pvals(k), slope_k);
    end
    
        if skip_count > 0 && false
        % x-bounds of excluded region
        x_fill = [min(values_eps), max(values_eps(1+skip_count:end))];
        % y-bounds: full vertical axis
        yl = [1e-9 1e2]; 
        ylim_old = get(gca,'YLim');
        patch([x_fill(1) x_fill(2) x_fill(2) x_fill(1)], ...
              [yl(1) yl(1) yl(2) yl(2)], ...
              'k', 'FaceAlpha', 0.05, 'EdgeColor', 'none');
        set(gca, 'YLim', ylim_old);  
        end
    uistack(findobj(gca,'Type','patch'),'bottom'); % send patch behind data
    legend(h, labels, 'Interpreter', 'latex', 'FontSize', 10, 'Location', 'southeast', 'Box', 'off');
    
    hold off;
end