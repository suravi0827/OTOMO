clear;
close all;
clc;

summaryFolder = fullfile(pwd,'Results','Metric_Summaries');
outputFolder = fullfile(pwd,'Results','Significance_Tests');
gammaTag = "";
datasets = ["AEEEM","RELINK","NASA","SOFTLAB"];
alpha = 0.05;
tieCorrection = false;
highlightMethod = "OTOMO";

classifierFiles = ["Random_Forest","Logistic_Regression"];
classifierTags = ["RF","LR"];
metricFiles = ["ROC_AUC","G_Measure","F1_Score","MCC_Score"];
metricLabels = ["ROC-AUC","G-mean","F1-score","MCC"];

assert(isfolder(summaryFolder),'Summary folder not found: %s',summaryFolder);

files = dir(fullfile(summaryFolder,'*_All_Datasets_*_gamma_*.csv'));
assert(~isempty(files),'No metric summary CSVs found in %s',summaryFolder);

tags = unique(string(regexp({files.name},'gamma_.+(?=\.csv$)','match','once')));

if strlength(gammaTag) == 0
    assert(isscalar(tags),'Several gamma settings found (%s); set gammaTag.',char(strjoin(tags,', ')));
    gammaTag = tags;
end

inputFiles = strings(numel(classifierFiles),numel(metricFiles));

for cc = 1:numel(classifierFiles)
    for mm = 1:numel(metricFiles)
        inputFiles(cc,mm) = fullfile(summaryFolder,sprintf('%s_All_Datasets_%s_%s.csv',classifierFiles(cc),metricFiles(mm),gammaTag));
        assert(isfile(inputFiles(cc,mm)),'Metric summary not found: %s',inputFiles(cc,mm));
    end
end

outputFolder = fullfile(outputFolder,char(gammaTag));

if ~exist(outputFolder,'dir')
    mkdir(outputFolder);
end

summaryRows = {};
summaryOrder = zeros(0,3);
allMethods = strings(1,0);
groups = [datasets,"All"];

for cc = 1:numel(classifierFiles)

    panelMethods = cell(1,numel(metricFiles));
    panelStats = cell(1,numel(metricFiles));

    for mm = 1:numel(metricFiles)

        file = inputFiles(cc,mm);

        [methods,family,values] = readMetricSummary(file);

        if isempty(allMethods)
            allMethods = methods;
            if ~any(allMethods == highlightMethod)
                warning('highlightMethod "%s" is not a column of the summaries (%s): nothing will be highlighted.', ...
                    highlightMethod,char(strjoin(allMethods,', ')));
            end
        end
        assert(isequal(methods,allMethods),'Method columns in %s differ from the other summary files.',file);

        familiesFound = unique(family,'stable').';
        missingFamilies = setdiff(datasets,familiesFound,'stable');
        assert(isempty(missingFamilies),'%s has no rows for %s (families found: %s).', ...
            file,char(strjoin(missingFamilies,', ')),char(strjoin(familiesFound,', ')));
        extraFamilies = setdiff(familiesFound,datasets,'stable');
        if ~isempty(extraFamilies)
            warning('%s: families not listed in datasets are excluded: %s',file,char(strjoin(extraFamilies,', ')));
        end

        fprintf('\n%s | %s | %s\n',classifierTags(cc),metricLabels(mm),file);

        for dd = 1:numel(groups)

            if groups(dd) == "All"
                rowMask = ismember(family,datasets);
            else
                rowMask = family == groups(dd);
            end

            [groupMethods,X] = completeBlocks(methods,values(rowMask,:),groups(dd));

            S = friedmanNemenyi(X,alpha,tieCorrection);

            fprintf('  %-8s eta=%2d L=%d | chi2F=%8.4f p=%.3g | FF=%8.4f Fcrit=%.4f reject H0=%d | CD=%.4f\n', ...
                groups(dd),S.N,S.k,S.chi2,S.pChi2,S.FF,S.Fcrit,S.rejectH0,S.cd);

            summaryRows(end+1,:) = summaryRow(metricLabels(mm),classifierTags(cc),groups(dd),S,groupMethods,methods,highlightMethod,alpha); %#ok<SAGROW>
            summaryOrder(end+1,:) = [mm cc dd]; %#ok<SAGROW>

            if groups(dd) == "All"
                panelMethods{mm} = groupMethods;
                panelStats{mm} = S;
            end

        end

    end

    nPanels = numel(metricFiles);
    kMax = max(cellfun(@(S) S.k,panelStats));
    cliqueRows = max(cellfun(@(S) size(findCliques(sort(S.avgRank),S.cd),1),panelStats));
    L = diagramLayout(kMax,cliqueRows,max(strlength(allMethods)));

    figureWidth = 12 * diff(L.xLim) / 10.2;
    figureHeight = ceil(nPanels / 2) * diff(L.yLim) * (figureWidth / 2) / diff(L.xLim) + 0.4;
    fig = figure('Visible','off','Color','w','Units','inches','Position',[1 1 figureWidth figureHeight],'InvertHardcopy','off');
    layout = tiledlayout(fig,ceil(nPanels / 2),2,'TileSpacing','tight','Padding','tight');

    for mm = 1:nPanels
        ax = nexttile(layout);
        caption = sprintf('(%c) %s',char('a' + mm - 1),metricLabels(mm));
        if ~panelStats{mm}.rejectH0
            caption = sprintf('%s (Friedman n.s., p = %.3f)',caption,panelStats{mm}.pFF);
        end
        drawCriticalDifference(ax,panelMethods{mm},panelStats{mm}.avgRank,panelStats{mm}.cd,highlightMethod,caption,L);
    end

    outputFile = fullfile(outputFolder,sprintf('%s_Overall_Friedman_Nemenyi_CD.eps',classifierTags(cc)));
    exportgraphics(fig,outputFile,'ContentType','vector','BackgroundColor','white');
    exportgraphics(fig,strrep(outputFile,'.eps','.png'),'Resolution',600,'BackgroundColor','white');
    close(fig);

    fprintf('\nSaved figure (all %d pairs, one panel per metric): %s (+ 600 dpi .png)\n',panelStats{1}.N,outputFile);

end

[~,rowOrder] = sortrows(summaryOrder);
summaryRows = summaryRows(rowOrder,:);

header = [{'Metric','Classifier','Dataset','N_Pairs','K_Methods','Friedman_Chi2','Friedman_p','Friedman_Significant', ...
    'ImanDavenport_F','ImanDavenport_Fcrit','ImanDavenport_p','Reject_H0','Nemenyi_q','CD','Best_Method'}, ...
    cellstr("AvgRank_" + allMethods), cellstr(highlightMethod + "_vs_" + allMethods(allMethods ~= highlightMethod))];

summaryFile = fullfile(outputFolder,'Friedman_Nemenyi_Summary.csv');
writecell([header;summaryRows],summaryFile);

fprintf('\nSummary table: %s\n',summaryFile);
fprintf('Figures saved in: %s\n',outputFolder);

function [methods,family,values] = readMetricSummary(file)

lines = splitlines(string(regexprep(fileread(file),'^\x{FEFF}','')));
lines = strtrim(lines(strlength(strtrim(lines)) > 0));

header = strtrim(split(lines(1),",")).';
assert(numel(header) >= 4 && header(1) == "Source" && header(2) == "Target",'Unexpected header in %s',file);

methods = header(3:end);
methods = methods(1:find(methods ~= "",1,'last'));
k = numel(methods);
assert(k >= 2 && all(methods ~= "") && numel(unique(methods)) == k,'Empty or duplicate method names in the header of %s',file);

family = strings(0,1);
pairKeys = strings(0,1);
values = zeros(0,k);
current = "";

for i = 2:numel(lines)

    cells = strtrim(split(lines(i),",")).';
    cells(end+1:k+2) = "";

    if cells(1) ~= "" && all(cells(2:end) == "")
        current = cells(1);
        continue;
    end

    if cells(1) == "" || cells(2) == ""
        continue;
    end

    assert(current ~= "",'Pair row before any dataset row in %s: %s',file,lines(i));

    row = nan(1,k);

    for j = 1:k
        token = regexp(cells(j+2),'^(NaN|[-+]?(\d+\.?\d*|\.\d+)([eE][-+]?\d+)?)','match','once');
        if strlength(token) > 0
            row(j) = str2double(token);
        end
    end

    family(end+1,1) = current; %#ok<AGROW>
    pairKeys(end+1,1) = current + ": " + cells(1) + " -> " + cells(2); %#ok<AGROW>
    values(end+1,:) = row; %#ok<AGROW>

end

assert(~isempty(values),'No source/target rows found in %s',file);

[~,firstRow] = unique(pairKeys,'stable');
duplicateRows = setdiff(1:numel(pairKeys),firstRow);
assert(isempty(duplicateRows),'Duplicate source/target rows in %s: %s',file,char(strjoin(unique(pairKeys(duplicateRows)),', ')));

end

function [methods,X] = completeBlocks(methods,X,groupName)

keepMethod = any(isfinite(X),1);

if any(~keepMethod)
    warning('%s: dropping methods with no finite values: %s',groupName,char(strjoin(methods(~keepMethod),', ')));
end

methods = methods(keepMethod);
X = X(:,keepMethod);

keepRow = all(isfinite(X),2);

if any(~keepRow)
    warning('%s: dropping %d source/target pairs with missing values.',groupName,sum(~keepRow));
end

X = X(keepRow,:);

assert(size(X,1) >= 2 && size(X,2) >= 2,'%s: not enough complete data for the Friedman test.',groupName);

end

function S = friedmanNemenyi(X,alpha,tieCorrection)

% Notation as in the paper: eta cross-project pairs (rows of X), L methods (columns),
% R_j the average rank of method j over all pairs (rank 1 = best, ties get average ranks).
[eta,L] = size(X);

ranks = zeros(eta,L);

for i = 1:eta
    ranks(i,:) = tiedrank(-X(i,:));
end

R = mean(ranks,1);

S.N = eta;
S.k = L;
S.avgRank = R;

% Friedman statistic: chi2_F = 12*eta/(L*(L+1)) * (sum_j R_j^2 - L*(L+1)^2/4), with L-1 degrees of freedom.
% tieCorrection = true instead divides by the tie-corrected rank variance (MATLAB friedman / R / SciPy form).
if ~tieCorrection
    S.chi2 = 12 * eta / (L * (L + 1)) * (sum(R .^ 2) - L * (L + 1) ^ 2 / 4);
else
    ssRank = sum(ranks(:) .^ 2) - eta * L * (L + 1) ^ 2 / 4;
    if ssRank > 0
        S.chi2 = (L - 1) * eta ^ 2 * sum((R - (L + 1) / 2) .^ 2) / ssRank;
    else
        S.chi2 = 0;
    end
end

S.pChi2 = chi2cdf(S.chi2,L - 1,'upper');

% Iman-Davenport: F_F = (eta-1)*chi2_F / (eta*(L-1) - chi2_F) against F(L-1, (L-1)*(eta-1));
% reject H0 when F_F exceeds the upper-tail critical value at alpha.
S.Fcrit = finv(1 - alpha,L - 1,(L - 1) * (eta - 1));

if eta * (L - 1) - S.chi2 > 0
    S.FF = (eta - 1) * S.chi2 / (eta * (L - 1) - S.chi2);
    S.pFF = fcdf(S.FF,L - 1,(L - 1) * (eta - 1),'upper');
else
    S.FF = Inf;
    S.pFF = 0;
end

S.rejectH0 = S.FF > S.Fcrit;

% Nemenyi: CD = q_{alpha,L} * sqrt(L*(L+1)/(6*eta)); two methods differ when |R_i - R_j| > CD.
S.q = nemenyiQ(L,alpha);
S.cd = S.q * sqrt(L * (L + 1) / (6 * eta));

end

function q = nemenyiQ(k,alpha)

rangeCdf = @(w) k * integral(@(z) normpdf(z) .* (normcdf(z) - normcdf(z - w)) .^ (k - 1),-Inf,Inf);

w = fzero(@(w) rangeCdf(w) - (1 - alpha),[0.01 20]);

q = w / sqrt(2);

end

function row = summaryRow(metric,classifier,groupName,S,groupMethods,allMethods,highlightMethod,alpha)

avgRank = nan(1,numel(allMethods));
[found,loc] = ismember(groupMethods,allMethods);
avgRank(loc(found)) = S.avgRank(found);

best = abs(S.avgRank - min(S.avgRank)) < 1e-9;

others = allMethods(allMethods ~= highlightMethod);
verdicts = repmat({'n/a'},1,numel(others));
h = find(groupMethods == highlightMethod,1);

if ~S.rejectH0
    verdicts(:) = {'not tested (Friedman n.s.)'};
elseif ~isempty(h)
    for j = 1:numel(others)
        o = find(groupMethods == others(j),1);
        if isempty(o)
            continue;
        end
        difference = S.avgRank(o) - S.avgRank(h);
        if difference > S.cd
            verdicts{j} = 'significantly better';
        elseif -difference > S.cd
            verdicts{j} = 'significantly worse';
        else
            verdicts{j} = 'no significant difference';
        end
    end
end

row = [{char(metric),char(classifier),char(groupName),S.N,S.k,S.chi2,S.pChi2,char(string(S.pChi2 < alpha)), ...
    S.FF,S.Fcrit,S.pFF,char(string(S.rejectH0)),S.q,S.cd,char(strjoin(groupMethods(best),' / '))},num2cell(avgRank),verdicts];

end

function cliques = findCliques(sortedRank,cd)

k = numel(sortedRank);
cliques = zeros(0,2);

for i = 1:k
    j = find(sortedRank - sortedRank(i) <= cd,1,'last');
    if j > i && (isempty(cliques) || j > cliques(end,2))
        cliques(end+1,:) = [i j]; %#ok<AGROW>
    end
end

end

function L = diagramLayout(k,cliqueRows,maxChars)

xPad = max(2.1,0.5 + 0.24 * maxChars);

L.rowGap = 0.34;
L.cliqueGap = 0.13;
L.yCD = 0.78;
L.yClique = -0.36;
L.yLabel = L.yClique - max(cliqueRows,1) * L.cliqueGap - 0.14;
L.yBottom = L.yLabel - (ceil(k / 2) - 1) * L.rowGap;
L.yCaption = L.yBottom - 0.42;
L.xLim = [-xPad k - 1 + xPad];
L.yLim = [L.yBottom - 0.95 L.yCD + 0.55];

end

function drawCriticalDifference(ax,methods,avgRank,cd,highlightMethod,caption,L)

fontName = 'Times New Roman';
k = numel(methods);

[sortedRank,order] = sort(reshape(avgRank,1,[]),'ascend');
sortedNames = methods(order);

xOf = @(r) k - r;

cliqueColors = [0 0 1; 1 0 0; 0.3 0.85 0.3; 0 0 0; 1 0.5 0; 0.6 0 0.8];

hold(ax,'on');
axis(ax,'off');

line(ax,[0 k-1],[0 0],'Color','k','LineWidth',1.2);

for r = 1:k
    line(ax,[xOf(r) xOf(r)],[0 -0.22],'Color','k','LineWidth',1.2);
    text(ax,xOf(r),0.07,sprintf('%d',r),'HorizontalAlignment','center','VerticalAlignment','bottom', ...
        'FontName',fontName,'FontSize',12,'FontWeight','bold','Color','k');
end

for r = 1.5:1:k-0.5
    line(ax,[xOf(r) xOf(r)],[0 -0.12],'Color','k','LineWidth',1.2);
end

yCD = L.yCD;
line(ax,[0 cd],[yCD yCD],'Color','k','LineWidth',1.6);
line(ax,[0 0],[yCD-0.14 yCD+0.14],'Color','k','LineWidth',1.6);
line(ax,[cd cd],[yCD-0.14 yCD+0.14],'Color','k','LineWidth',1.6);
text(ax,cd/2,yCD+0.16,sprintf('CD=%.4f',cd),'HorizontalAlignment','center','VerticalAlignment','bottom', ...
    'FontName',fontName,'FontSize',11,'FontWeight','bold','Color','k');

cliques = findCliques(sortedRank,cd);

for c = 1:size(cliques,1)
    y = L.yClique - (c - 1) * L.cliqueGap;
    xRight = xOf(sortedRank(cliques(c,1))) + 0.08;
    xLeft = xOf(sortedRank(cliques(c,2))) - 0.08;
    line(ax,[xLeft xRight],[y y],'Color',cliqueColors(mod(c-1,size(cliqueColors,1))+1,:),'LineWidth',2.8);
end

nRight = ceil(k / 2);
xLabelRight = k - 1 + 0.3;
xLabelLeft = -0.3;

for i = 1:k

    x = xOf(sortedRank(i));

    if i <= nRight
        y = L.yLabel - (i - 1) * L.rowGap;
        xEnd = xLabelRight;
        textX = xEnd + 0.08;
        alignment = 'left';
    else
        y = L.yLabel - (k - i) * L.rowGap;
        xEnd = xLabelLeft;
        textX = xEnd - 0.08;
        alignment = 'right';
    end

    if sortedNames(i) == highlightMethod
        labelColor = [0 0 1];
    else
        labelColor = [0 0 0];
    end

    line(ax,[x x],[0 y],'Color','k','LineWidth',0.9);
    line(ax,[x xEnd],[y y],'Color','k','LineWidth',0.9);
    text(ax,textX,y,char(sortedNames(i)),'HorizontalAlignment',alignment,'VerticalAlignment','middle', ...
        'FontName',fontName,'FontSize',12,'FontWeight','bold','Color',labelColor);

end

text(ax,(k - 1) / 2,L.yCaption,caption,'HorizontalAlignment','center','VerticalAlignment','top', ...
    'FontName',fontName,'FontSize',16,'Color','k');

set(ax,'DataAspectRatio',[1 1 1],'XLim',L.xLim,'YLim',L.yLim);
hold(ax,'off');

end
