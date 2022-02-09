% plot marginals

function [maxp, exp_true, lbf, ubf] = plot_marg(mu, sigma, exp_idx, lb, ub, data, av, allcolors, min_glc, max_glc)

    glc = data(exp_idx).glc;
    Ndata = length(data);
    glc_all = linspace(min_glc, max_glc, Ndata);
    dist = Inf;
    for i = 1:Ndata
        if ((glc - glc_all(i))^2 < dist)
            dist = (glc - glc_all(i))^2;
            i_color = i;
        end
    end
    colormax = allcolors(i_color,:);
    x0 = lb;
    x1 = ub;
    npoints = 1e5;
    exp_true = 0;
    pd = makedist('Normal','mu',mu,'sigma',sigma);
    xlinsp = linspace(lb,ub, npoints);

    try
        t = truncate(pd, lb, ub);
        aux =  pdf(t, xlinsp);
        lbf = min(xlinsp(find(aux > 1e-6)));
        ubf = max(xlinsp(find(aux > 1e-6)));
    catch
        aux = xlinsp .* 0;
        if( sum(isnan(aux)) || sum(isinf(aux)) || nnz(aux) == 0 )
            npoints = 5e4;
            xlinsp = linspace(x0,x1, npoints);
            pd = makedist('Exponential','mu',av);
            t = truncate(pd, x0, x1);
            exp_true = 1;
            aux =  pdf(t, xlinsp);
            lbf = min(xlinsp(find(aux > 1e-6)));
            ubf = max(xlinsp(find(aux > 1e-6)));
        end
    end
    plot(xlinsp, pdf(t, xlinsp),'Color',colormax);
    hold on
    area(xlinsp, pdf(t, xlinsp),'LineStyle','-','EdgeColor', colormax,'FaceColor',colormax,'FaceAlpha',0.5)
    hold on
    title('EP marginal')
    [maxp, aux] = max( pdf(t, xlinsp));
    %idx_maxp = xlinsp(aux);
    
end
