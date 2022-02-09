
function plot_flux(idx, flux_name, cmap, data, mut, st, av, va, min_glc, max_glc)

    figure
    Ndata = length(data);
    for i = 1:Ndata
        [maxp(i), exp_true(i), ~, ~] = plot_marg(mut(idx,i), sqrt(st(idx,i)), i, data(i).lb(idx), data(i).ub(idx), ...
            data, av(idx,i), cmap, min_glc, max_glc);
        hold on
    end
    if any(exp_true == 1)
        ylim([0,max(maxp)]);
    elseif any(exp_true == 0)
        ylim([0, max(maxp(find(~exp_true)))  + std(maxp(find(~exp_true)))]);
    else
        ;
    end
    if any(exp_true == 0)
        lb = min(av(idx,:) - 3*sqrt(va(idx,:)));
        ub = max(av(idx,:) + 3*sqrt(va(idx,:))); 
    else
        lb = mean(av(idx,:) - 5*sqrt(va(idx,:)));
        ub = mean(av(idx,:) + 5*sqrt(va(idx,:)));
    end
    xlim([lb,ub])
    colormap(cmap)
    h = colorbar('ticks',linspace(0,1,10),'yticklabel', round(linspace(min_glc, max_glc, 10), 2)); 
    ylabel(h, 'Glucose uptake')
    h0 = sprintf('\\nu_{%s}', flux_name);
    xlabel(h0);
    h0 = sprintf('P(\\nu_{%s})', flux_name);
    ylabel(h0);

end