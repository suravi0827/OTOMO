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
