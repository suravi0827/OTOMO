function [X_bal,Y_bal] = dpp_oversample(X,Y,ratio,sigma)

if nargin < 3 || isempty(ratio)
    ratio = 1.0;
end

X = real(double(X));
Y = double(Y(:));

idxPos = find(Y == 1);
idxNeg = find(Y == 0);

Xp = X(idxPos,:);
Xn = X(idxNeg,:);

np = size(Xp,1);
nn = size(Xn,1);

if np == 0 || nn == 0

    X_bal = X;
    Y_bal = Y;
    return;

end

desiredMinority = floor(ratio * nn);

G = desiredMinority - np;

if G <= 0

    X_bal = X;
    Y_bal = Y;
    return;

end


D = pdist2(Xp,Xp,'euclidean');


if nargin < 4 || isempty(sigma)

    vals = D(triu(true(size(D)),1));

    vals = vals(isfinite(vals) & vals > 0);

    if isempty(vals)

        sigma = 0.1;

    else

        sigma = median(vals);

    end

end

sigma = max(sigma,eps);


S = exp( ...
    -(D.^2) ./ ...
    (2*sigma^2));


q = ones(np,1);


L = ...
    diag(q) * ...
    S * ...
    diag(q);

%% Numerical stabilization

L = ...
    (L + L') ./ 2;

L = ...
    L + ...
    1e-10 .* eye(np);


K = ...
    L / ...
    (eye(np) + L);

prob = real(diag(K));

prob(~isfinite(prob) | prob < 0) = 0;

if sum(prob) <= eps

    prob = ones(np,1) ./ np;

else

    prob = prob ./ sum(prob);

end

sampled = zeros( ...
    G, ...
    size(X,2));

for g = 1:G

    idx = weightedRandomIndex(prob);

    sampled(g,:) = Xp(idx,:);

end

X_bal = ...
    [X;
     sampled];

Y_bal = ...
    [Y;
     ones(G,1)];

end


function idx = weightedRandomIndex(prob)

prob = double(prob(:));

prob(~isfinite(prob) | prob < 0) = 0;

if sum(prob) <= eps

    prob = ...
        ones(numel(prob),1) ./ ...
        numel(prob);

else

    prob = ...
        prob ./ sum(prob);

end

cdf = cumsum(prob);

r = rand();

idx = find( ...
    r <= cdf, ...
    1, ...
    'first');

if isempty(idx)
    idx = numel(prob);
end

end