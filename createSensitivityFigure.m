function createSensitivityFigure(gammas, metricValues, datasets, markers, yLabelText, outputFile)
fig = figure('Color','w', 'Position',[100 100 800 600]);
ax = axes(fig);
ax.Color = 'w';
hold(ax,'on');

for d = 1:numel(datasets)
    plot(ax, gammas, metricValues(d,:), '-', 'Marker',markers{d}, 'LineWidth',1.8, 'MarkerSize',8, 'DisplayName',char(datasets(d)));
end

xlabel(ax, '\gamma', 'Interpreter','tex', 'FontSize',14, 'FontWeight','bold', 'Color','k');
ylabel(ax, yLabelText, 'FontSize',14, 'FontWeight','bold', 'Color','k');
xticks(ax,gammas);
xticklabels(ax,compose('%.1f',gammas));
xlim(ax,[0 1]);

allValues = metricValues(:);
allValues = allValues(isfinite(allValues));

if ~isempty(allValues)
    minY = min(allValues);
    maxY = max(allValues);
    rangeY = maxY - minY;

    if rangeY < eps
        rangeY = 0.1;
    end

    padding = 0.10 * rangeY;
    lowerLimit = minY - padding;
    upperLimit = maxY + padding;

    if ~strcmpi(yLabelText,'MCC')
        lowerLimit = max(0,lowerLimit);
        upperLimit = min(1,upperLimit);
    else
        lowerLimit = max(-1,lowerLimit);
        upperLimit = min(1,upperLimit);
    end

    if lowerLimit < upperLimit
        ylim(ax,[lowerLimit upperLimit]);
    end
end

ax.Color = 'w';
ax.XColor = 'k';
ax.YColor = 'k';
ax.FontSize = 14;
ax.LineWidth = 1;
ax.TickDir = 'out';
ax.GridColor = [0.75 0.75 0.75];
ax.GridAlpha = 0.35;
grid(ax,'on');
box(ax,'on');

lgd = legend(ax, 'Location','northeast', 'FontSize',11, 'Box','on');
lgd.Color = 'w';
lgd.TextColor = 'k';
lgd.EdgeColor = 'k';

fig.Color = 'w';
set(fig, 'InvertHardcopy','off', 'PaperPositionMode','auto');
hold(ax,'off');
exportgraphics(fig, outputFile, 'BackgroundColor','white', 'Resolution',300);
end