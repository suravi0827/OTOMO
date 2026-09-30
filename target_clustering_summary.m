clear;
clc;
projectRoot = pwd;
datasets = ["NASA","RELINK","SOFTLAB","AEEEM"];
seeds = 1:30;
outputFolder = fullfile(projectRoot,'target clustering summary');
assert(isfolder(fullfile(projectRoot,'data')),'Run this script from the project folder containing data.');
assert(~isempty(seeds) && isvector(seeds) && numel(unique(seeds)) == numel(seeds),'Provide distinct seeds.');
assert(all(isfinite(seeds) & seeds >= 0 & seeds <= 2^32-1 & seeds == floor(seeds)),'Seeds must be integers between 0 and 2^32-1.');
fileLists = cell(numel(datasets),1);
nProjects = 0;
for dd = 1:numel(datasets)
    files = dir(fullfile(projectRoot,'data',char(datasets(dd)),'*.xlsx'));
    files = files(~startsWith({files.name},'~$'));
    [~,order] = sort({files.name});
    fileLists{dd} = files(order);
    nProjects = nProjects + numel(files);
end
assert(nProjects > 0,'No project spreadsheets found.');
header = {'Dataset','Project','Purity (%) mean +/- SD','ARI (%) mean +/- SD','NMI (%) mean +/- SD'};
result = cell(nProjects+1,numel(header));
result(1,:) = header;
row = 1;
for dd = 1:numel(datasets)
    files = fileLists{dd};
    if isempty(files)
        warning('No project files found for %s.',char(datasets(dd)));
    end
    for pp = 1:numel(files)
        row = row + 1;
        [~,projectName] = fileparts(files(pp).name);
        result{row,1} = char(datasets(dd));
        result{row,2} = projectName;
        scores = nan(numel(seeds),3);
        fprintf('%s | %s | %d seeds\n',char(datasets(dd)),projectName,numel(seeds));
        try
            A = readmatrix(fullfile(files(pp).folder,files(pp).name));
            assert(size(A,2) >= 2,'Input must contain features and a final label column.');
            A = A(all(isfinite(A),2),:);
            assert(size(A,1) >= 2,'At least two valid rows are required.');
            Xt = safeZscore(A(:,1:end-1));
            Yt = double(A(:,end));
            assert(all(ismember(Yt,[0 1])),'Labels must be binary 0 and 1.');
            for rr = 1:numel(seeds)
                try
                    [~,~,purity,ari,nmi] = otomo_clustering(Xt,Yt,seeds(rr));
                    scores(rr,:) = [purity,ari,nmi];
                catch ME
                    warning('%s | %s | seed=%d failed: %s',char(datasets(dd)),projectName,seeds(rr),ME.message);
                end
            end
        catch ME
            warning('%s | %s failed: %s',char(datasets(dd)),projectName,ME.message);
        end
        for metric = 1:3
            formatted = formatPercentageMeanSD(scores(:,metric));
            result{row,metric+2} = formatted;
        end
    end
end
if ~isfolder(outputFolder)
    mkdir(outputFolder);
end
summaryFile = fullfile(outputFolder,'OTOMO_Project_Clustering_Summary.csv');
writecell(result,summaryFile);
fprintf('Summary saved: %s\n',summaryFile);

function Z = safeZscore(X)
X = double(X);
mu = mean(X,1,'omitnan');
sd = std(X,0,1,'omitnan');
sd(~isfinite(sd) | sd < eps) = 1;
Z = (X-mu)./sd;
Z(~isfinite(Z)) = 0;
end

function [clusterIdx,Tmin,purity,ARI,NMI] = otomo_clustering(Xt,Yt,seed)
if nargin < 2
    Yt = [];
end
if nargin < 3 || isempty(seed)
    seed = 1;
end
Xt = real(double(Xt));
assert(ismatrix(Xt) && size(Xt,1) >= 2 && size(Xt,2) >= 1,'Xt must contain at least two rows and one feature.');
assert(all(isfinite(Xt(:))),'Xt must contain only finite values.');
assert(isempty(Yt) || numel(Yt) == size(Xt,1),'Yt must contain one label per target row.');
purity = NaN;
ARI = NaN;
NMI = NaN;
rng(seed,'twister');
clusterIdx = kmeans(Xt,2,'Distance','sqeuclidean','Replicates',20,'MaxIter',1000,'Display','off');
clusterIdx = double(clusterIdx(:));
clusterSizes = [sum(clusterIdx == 1),sum(clusterIdx == 2)];
[~,minorityCluster] = min(clusterSizes);
Tmin = Xt(clusterIdx == minorityCluster,:);
if ~isempty(Yt)
    [purity,ARI,NMI] = clusteringMetrics(clusterIdx,Yt);
end
end

function [purity,ARI,NMI] = clusteringMetrics(predicted,original)
predicted = double(predicted(:));
original = double(original(:));
valid = isfinite(predicted) & isfinite(original);
predicted = predicted(valid);
original = original(valid);
if isempty(original)
    purity = NaN;
    ARI = NaN;
    NMI = NaN;
    return;
end
predClass = unique(predicted);
trueClass = unique(original);
C = zeros(numel(predClass),numel(trueClass));
for i = 1:numel(predClass)
    for j = 1:numel(trueClass)
        C(i,j) = sum(predicted == predClass(i) & original == trueClass(j));
    end
end
n = sum(C(:));
purity = sum(max(C,[],2)) ./ n;
c2 = @(x) x .* (x - 1) ./ 2;
cellPairs = sum(c2(C(:)));
rowPairs = sum(c2(sum(C,2)));
colPairs = sum(c2(sum(C,1)));
totalPairs = c2(n);
expected = rowPairs .* colPairs ./ max(totalPairs,eps);
maxIndex = (rowPairs + colPairs) ./ 2;
den = maxIndex - expected;
if abs(den) <= eps
    ARI = double(abs(cellPairs - expected) <= eps);
else
    ARI = (cellPairs - expected) ./ den;
end
P = C ./ n;
pPred = sum(P,2);
pTrue = sum(P,1);
MI = 0;
for i = 1:size(P,1)
    for j = 1:size(P,2)
        if P(i,j) > 0
            MI = MI + P(i,j) .* log(P(i,j) ./ (pPred(i) .* pTrue(j)));
        end
    end
end
Hpred = -sum(pPred(pPred > 0) .* log(pPred(pPred > 0)));
Htrue = -sum(pTrue(pTrue > 0) .* log(pTrue(pTrue > 0)));
den = sqrt(Hpred .* Htrue);
if den <= eps
    NMI = double(Hpred <= eps && Htrue <= eps);
else
    NMI = MI ./ den;
end
end

function formatted = formatPercentageMeanSD(values)

values = double(values(:));
values = values(isfinite(values));
nValid = numel(values);

if nValid == 0
    formatted = 'NaN +/- NaN';
    return;
end

meanPercentage = 100 * mean(values);

if nValid < 2
    formatted = sprintf('%.2f +/- NaN',meanPercentage);
    return;
end

sdPercentage = 100 * std(values,0);
formatted = sprintf('%.2f +/- %.2f',meanPercentage,sdPercentage);

end
