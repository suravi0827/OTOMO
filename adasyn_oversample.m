function [X_bal, Y_bal] = adasyn_oversample(X, Y, ratio, k)

%% ============================================================
% ADASYN oversampling
%
% X     : feature matrix
% Y     : binary labels (0 = majority, 1 = minority)
% ratio : desired minority / majority ratio after oversampling
%         ratio = 1.0 -> full balance
% k     : number of nearest neighbors
%
% Example:
% [Xb,Yb] = adasyn_oversample(Xs,Ys,1.0,5);
% ============================================================

if nargin < 3 || isempty(ratio)
    ratio = 1.0;
end

if nargin < 4 || isempty(k)
    k = 5;
end

X = real(double(X));
Y = double(Y(:));

%% ------------------------------------------------------------
% Separate classes
% -------------------------------------------------------------

idxPos = find(Y == 1);
idxNeg = find(Y == 0);

Xp = X(idxPos,:);
Xn = X(idxNeg,:);

np = size(Xp,1);
nn = size(Xn,1);

%% ------------------------------------------------------------
% Basic checks
% -------------------------------------------------------------

if np < 2 || nn == 0

    X_bal = X;
    Y_bal = Y;
    return;

end

%% ------------------------------------------------------------
% Number of synthetic samples
% -------------------------------------------------------------

desiredMinority = floor(nn * ratio);

G = desiredMinority - np;

if G <= 0

    X_bal = X;
    Y_bal = Y;
    return;

end

%% ============================================================
% 1. Determine difficulty of each minority sample
%
% r_i = number of majority neighbors / k
% ============================================================

n = size(X,1);

kUse = min(k, n-1);

Dall = pdist2(Xp, X, 'euclidean');

r = zeros(np,1);

for i = 1:np

    % Exclude the same sample from itself
    originalIdx = idxPos(i);

    Dall(i,originalIdx) = inf;

    [~,ord] = sort(Dall(i,:),'ascend');

    neigh = ord(1:kUse);

    r(i) = sum(Y(neigh) == 0) / kUse;

end

%% ------------------------------------------------------------
% If every minority point is surrounded only by minority points,
% fall back to equal generation probability
% -------------------------------------------------------------

if sum(r) <= eps

    r = ones(np,1) / np;

else

    r = r / sum(r);

end

%% ============================================================
% 2. Number of synthetic samples assigned to each minority point
% ============================================================

g = round(r * G);

%% ------------------------------------------------------------
% Correct rounding so total equals G
% -------------------------------------------------------------

difference = G - sum(g);

while difference ~= 0

    if difference > 0

        [~,ord] = sort(r,'descend');

        for j = 1:min(difference,np)
            g(ord(j)) = g(ord(j)) + 1;
        end

    else

        [~,ord] = sort(g,'descend');

        remaining = -difference;

        for j = 1:np

            if remaining <= 0
                break;
            end

            ii = ord(j);

            if g(ii) > 0

                g(ii) = g(ii) - 1;
                remaining = remaining - 1;

            end

        end

    end

    difference = G - sum(g);

end

%% ============================================================
% 3. Minority-to-minority neighbors
% ============================================================

Dpos = pdist2(Xp, Xp, 'euclidean');

Dpos(1:np+1:end) = inf;

kMinor = min(k,np-1);

[~,posOrder] = sort(Dpos,2,'ascend');

%% ============================================================
% 4. Generate synthetic samples
% ============================================================

synthetic = zeros(G,size(X,2));

cnt = 0;

for i = 1:np

    for j = 1:g(i)

        cnt = cnt + 1;

        % Random minority neighbor
        nbrList = posOrder(i,1:kMinor);

        nbr = nbrList(randi(kMinor));

        alpha = rand();

        synthetic(cnt,:) = ...
            Xp(i,:) + ...
            alpha .* ...
            (Xp(nbr,:) - Xp(i,:));

    end

end

%% ============================================================
% 5. Output
% ============================================================

X_bal = [X; synthetic];

Y_bal = [Y; ones(G,1)];

end