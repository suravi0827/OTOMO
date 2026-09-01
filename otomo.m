function [Xs_bal,Ys_bal,purity,ARI,NMI] = otomo(Xs,Ys,Xt,Yt,gamma,ratio,conf)

if nargin < 5 || isempty(gamma)
    gamma = 0.5;
end

if nargin < 6 || isempty(ratio)
    ratio = 1.0;
end

if nargin < 7
    beta = [];
end

if nargin < 8 || isempty(conf)
    conf = 0.995;
end

gamma = max(0,min(1,gamma));
ratio = max(0,ratio);

purity = NaN;
ARI = NaN;
NMI = NaN;

Xs = real(double(Xs));
Xt = real(double(Xt));
Ys = double(Ys(:));

if ~isempty(Yt)
    Yt = double(Yt(:));
end

idxPos = find(Ys == 1);
idxNeg = find(Ys == 0);
Smin = Xs(idxPos,:);
Nmin = numel(idxPos);
Nmaj = numel(idxNeg);

if Nmin < 2 || Nmaj == 0 || size(Xt,1) < 2
    Xs_bal = Xs;
    Ys_bal = Ys;
    return;
end

Nsyn = max(0,floor(Nmaj * ratio) - Nmin);

if Nsyn == 0
    Xs_bal = Xs;
    Ys_bal = Ys;
    return;
end

clusterIdx = kmeans(Xt,2,'Distance','sqeuclidean','Replicates',20,'MaxIter',1000,'Display','off');
clusterIdx = double(clusterIdx(:));

clusterSizes = [sum(clusterIdx == 1),sum(clusterIdx == 2)];
[~,minorityCluster] = min(clusterSizes);
idxTargetMin = clusterIdx == minorityCluster;
Tmin = Xt(idxTargetMin,:);

if isempty(Tmin)
    Xs_bal = Xs;
    Ys_bal = Ys;
    return;
end

if ~isempty(Yt) && numel(Yt) == numel(clusterIdx)
    [purity,ARI,NMI] = clusteringMetrics(clusterIdx,Yt);
    fprintf('K-means agreement | Purity=%.4f | ARI=%.4f | NMI=%.4f\n',purity,ARI,NMI);
end

[isOutlier,md2Target,thrTarget] = mahalanobisOutliersRegularized(Tmin,conf);

TminClean = Tmin(~isOutlier,:);

if isempty(TminClean)
    warning('All target pseudo-minority samples exceeded the Mahalanobis reference threshold. Using the original pseudo-minority cluster.');
    TminClean = Tmin;
end

targetCentroid = mean(TminClean,1);

Ds = pdist2(Smin,Smin,'euclidean');
dC = sqrt(sum((Smin - targetCentroid).^2,2));

N = size(Smin,1);

Dt = 0.5 .* (repmat(dC,1,N) + repmat(dC',N,1));
H = gamma .* Ds + (1 - gamma) .* Dt;

H(1:N+1:end) = inf;

[~,neighbors] = sort(H,2,'ascend');

synthetic = zeros(Nsyn,size(Xs,2));

count = 0;

while count < Nsyn
    for i = 1:N
        if count >= Nsyn
            break;
        end
        nn_idx = randi([1, 7]);
                j = neighbors(i,nn_idx);
        r = rand;
        xSyn = Smin(j,:) + r .* (Smin(i,:) - Smin(j,:));
        count = count + 1;
        synthetic(count,:) = xSyn;
    end
end

Xs_bal = [Xs;synthetic];
Ys_bal = [Ys;ones(Nsyn,1)];

fprintf('OTOMO | gamma=%.2f | ratio=%.2f | minority=%d | majority=%d | synthetic=%d\n',gamma,ratio,Nmin,Nmaj,Nsyn);
fprintf('Target cluster sizes: C1=%d | C2=%d | selected=%d | pseudo-minority=%d/%d\n',clusterSizes(1),clusterSizes(2),minorityCluster,size(Tmin,1),size(Xt,1));
fprintf('Target pseudo-minority Mahalanobis outliers: %d/%d | threshold=%.6f\n',sum(isOutlier),numel(isOutlier),thrTarget);
fprintf('Cleaned target pseudo-minority: %d/%d retained\n',size(TminClean,1),size(Tmin,1));

end


function [isOutlier,md2,thr] = mahalanobisOutliersRegularized(X,conf)

X = real(double(X));
[n,d] = size(X);

md2 = zeros(n,1);
isOutlier = false(n,1);

thr = chi2inv(conf,max(d,1));

if n < 2 || d < 1
    return;
end

mu = mean(X,1);
Z = X - mu;
S = cov(X);

if d == 1
    S = reshape(S,1,1);
end

S(~isfinite(S)) = 0;
S = (S + S') ./ 2;

scale = trace(S) ./ max(d,1);

if ~isfinite(scale) || scale <= eps
    featureVar = var(X,0,1);
    validVar = featureVar(isfinite(featureVar) & featureVar > eps);

    if isempty(validVar)
        scale = 1;
    else
        scale = mean(validVar);
    end
end

rho = 1e-4;
ridgeValue = rho .* scale;
ridgeValue = max(ridgeValue,1e-12);
Sreg = S + ridgeValue .* eye(d);
Sinv = pinv(Sreg);

md2 = sum((Z * Sinv) .* Z,2);
md2(~isfinite(md2) | md2 < 0) = 0;

isOutlier = md2 > thr;

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