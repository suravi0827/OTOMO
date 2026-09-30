function [Xs_bal,Ys_bal] = tomo(Xs,Ys,Xt,lambda,ratio,seed)

if nargin < 4 || isempty(lambda)
    lambda = 0.4;
end

if nargin < 5 || isempty(ratio)
    ratio = 1.0;
end

if nargin < 6 || isempty(seed)
    seed = 1;
end

k = 7;

lambda = max(0,min(1,lambda));

Xs = real(double(Xs));
Xt = real(double(Xt));
Ys = double(Ys(:));

Xs_bal = Xs;
Ys_bal = Ys;

Xp = Xs(Ys == 1,:);
Xn = Xs(Ys == 0,:);

np = size(Xp,1);
nn = size(Xn,1);

if np < 2 || nn == 0 || size(Xt,1) < 2
    return;
end

n0 = floor(nn * ratio) - np;

if n0 <= 0
    return;
end

kEffective = min(k,np - 1);

rng(seed,'twister');

clusterIdx = kmeans(Xt,2,'Distance','sqeuclidean','Replicates',20,'MaxIter',1000);

clusterSize1 = sum(clusterIdx == 1);
clusterSize2 = sum(clusterIdx == 2);

if clusterSize1 <= clusterSize2
    minorityCluster = 1;
else
    minorityCluster = 2;
end

XtMinorityPotential = Xt(clusterIdx == minorityCluster,:);

Cmin = mean(XtMinorityPotential,1);

distanceToTarget = pdist2(Xp,Cmin,'euclidean');

[distanceToTarget,order] = sort(distanceToTarget,'ascend');

SP = Xp(order,:);

distSS = pdist2(SP,SP,'euclidean');
distTS = repmat(distanceToTarget',np,1);

distSSN = normalizeDistanceRows(distSS);
distTSN = normalizeDistanceRows(distTS);

distH = lambda .* distSSN + (1 - lambda) .* distTSN;

distH(1:np+1:end) = inf;

[~,neighborOrder] = sort(distH,2,'ascend');

NeigInd = neighborOrder(:,1:kEffective);

synthetic = zeros(n0,size(Xs,2));

for s = 1:n0

    i = mod(s - 1,np) + 1;

    neighborPosition = randi(kEffective);

    nbr = NeigInd(i,neighborPosition);

    r = rand();

    synthetic(s,:) = SP(i,:) - r .* (SP(nbr,:) - SP(i,:));

end

Xs_bal = [Xs; synthetic];
Ys_bal = [Ys; ones(n0,1)];

fprintf('TOMO | k=%d | lambda=%.2f | ratio=%.2f | minority=%d | majority=%d | synthetic=%d\n',kEffective,lambda,ratio,np,nn,n0);

end

function D = normalizeDistanceRows(D)

rowSum = sum(D,2);

rowSum(~isfinite(rowSum) | rowSum <= eps) = 1;

D = D ./ rowSum;

D(~isfinite(D)) = 0;

end