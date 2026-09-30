function plot_sensitivity_figures(csvFile,figureFolder,datasetOrder)

if nargin < 3 || isempty(datasetOrder)
    datasetOrder = ["NASA","RELINK","SOFTLAB","AEEEM"];
end

datasetOrder = string(datasetOrder(:).');

assert(isfile(csvFile),'Summary CSV not found: %s',csvFile);

T = readtable(csvFile,'TextType','string','VariableNamingRule','preserve');

assert(all(ismember(["Dataset","Gamma"],string(T.Properties.VariableNames))),'CSV must contain Dataset and Gamma columns.');
assert(height(T) > 0,'Summary CSV is empty: %s',csvFile);

T.Dataset = strtrim(string(T.Dataset));
T.Gamma = round(double(T.Gamma),10);

[~,~,groupIdx] = unique(T.Dataset + "|" + string(T.Gamma));
assert(max(accumarray(groupIdx,1)) == 1,'CSV contains duplicate Dataset/Gamma rows.');

csvDatasets = unique(T.Dataset,'stable').';
datasets = [datasetOrder(ismember(datasetOrder,csvDatasets)),csvDatasets(~ismember(csvDatasets,datasetOrder))];
gammas = unique(T.Gamma).';

assert(all(isfinite(gammas) & gammas >= 0 & gammas <= 1),'Gamma values must lie in [0,1].');

fprintf('Loaded %s\n',csvFile);
fprintf('Datasets: %s | gammas: %s\n',char(strjoin(datasets,', ')),mat2str(gammas));

if ~exist(figureFolder,'dir')
    mkdir(figureFolder);
end

classifiers = ["RF","LR"];
metricLabels = {'ROC-AUC','G-mean','F1-score','MCC'};
metricFiles = {'ROC_AUC','G_Mean','F1','MCC'};
metricColumns = {"ROC_AUC",["G_Mean","G_Measure"],"F1","MCC"};
markers = {'s','o','^','d'};

for cc = 1:numel(classifiers)
    for mm = 1:numel(metricLabels)
        column = findMetricColumn(T,classifiers(cc),metricColumns{mm});
        metricValues = nan(numel(datasets),numel(gammas));
        for r = 1:height(T)
            d = find(datasets == T.Dataset(r),1);
            g = find(gammas == T.Gamma(r),1);
            metricValues(d,g) = double(T.(column)(r));
        end
        outputFile = fullfile(figureFolder,sprintf('%s_Gamma_%s.eps', ...
            char(classifiers(cc)),metricFiles{mm}));
        createSensitivityFigure(gammas,metricValues,datasets, ...
            markers,metricLabels{mm},outputFile);
        fprintf('Saved figure: %s (+ 600 dpi .png)\n',outputFile);
    end
end

fprintf('\nFigures saved in:\n%s\n',figureFolder);

end

function column = findMetricColumn(T,classifier,suffixes)
names = string(T.Properties.VariableNames);
candidates = classifier + "_" + suffixes;
column = candidates(find(ismember(candidates,names),1));
assert(~isempty(column),'Column %s not found in the CSV.',char(strjoin(candidates," or ")));
column = char(column);
end

function createSensitivityFigure(gammas, metricValues, datasets, markers, yLabelText, outputFile)
colors = [0 0 0; 0.85 0 0; 0 0 0.85; 0 0.55 0];
fig = figure('Color','w','Units','inches','Position',[1 1 6 4.5],'InvertHardcopy','off');
ax = axes('Parent',fig,'Units','inches','Position',[1.0 0.8 4.8 3.5]);
hold(ax,'on');
for d = 1:numel(datasets)
    c = colors(mod(d-1,size(colors,1))+1,:);
    marker = markers{mod(d-1,numel(markers))+1};
    plot(ax,gammas,metricValues(d,:),'-','Color',c,'LineWidth',2,'Marker',marker,'MarkerSize',8,'MarkerFaceColor',c,'MarkerEdgeColor',c,'DisplayName',char(datasets(d)));
end
set(ax,'FontName','Arial','FontSize',14,'LineWidth',1.5,'Box','on','TickDir','in','TickLength',[0.015 0.015],'XColor','k','YColor','k','Color','w','Layer','top','XMinorTick','on','YMinorTick','on');
xlabel(ax,'\it\gamma','FontName','Arial','FontSize',16,'Color','k');
ylabel(ax,yLabelText,'FontName','Arial','FontSize',16,'Color','k');
xlim(ax,[-0.04 1.04]);
xticks(ax,0:0.2:1);
xticklabels(ax,compose('%.1f',0:0.2:1));
ax.XAxis.MinorTickValues = 0:0.1:1;
lgd = legend(ax,'Location','northeast','NumColumns',2,'FontName','Arial','FontSize',13,'TextColor','k');
set(lgd,'Box','on','EdgeColor','k','Color','w','LineWidth',1);
lgd.Units = 'inches';
drawnow;
values = metricValues(isfinite(metricValues));
if ~isempty(values)
    yMin = min(values);
    yMax = max(values);
    span = max(yMax - yMin,max(abs(yMax),1) * 0.01);
    yLow = yMin - 0.08 * span;
    legendBand = (ax.Position(2) + ax.Position(4) - lgd.Position(2)) / ax.Position(4) + 0.04;
    yTop = (yMax - legendBand * yLow) / (1 - legendBand);
    step = niceStep((yTop - yLow) / 6);
    yLow = floor(yLow / step) * step;
    yTop = ceil(yTop / step) * step;
    while yMax > yTop - legendBand * (yTop - yLow)
        yTop = yTop + step;
    end
    ylim(ax,[yLow yTop]);
    yticks(ax,yLow:step:yTop);
    decimals = 0;
    while abs(round(step * 10^decimals) - step * 10^decimals) > 1e-9
        decimals = decimals + 1;
    end
    yticklabels(ax,compose(sprintf('%%.%df',decimals),yLow:step:yTop));
    ax.YAxis.MinorTickValues = yLow:step/2:yTop;
end
hold(ax,'off');
drawnow;
[folder,name] = fileparts(outputFile);
exportgraphics(fig,outputFile,'ContentType','vector','BackgroundColor','white');
exportgraphics(fig,fullfile(folder,[name '.png']),'Resolution',600,'BackgroundColor','white');
close(fig);
end

function step = niceStep(rawStep)
magnitude = 10^floor(log10(rawStep));
candidates = [1 2 2.5 5 10] * magnitude;
step = candidates(find(candidates >= rawStep,1));
end
