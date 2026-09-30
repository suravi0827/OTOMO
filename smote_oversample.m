function [X_bal, Y_bal] = smote_oversample(X, Y, ratio, k)

%% ============================================================
% SMOTE oversampling
%
% X     : feature matrix
% Y     : binary labels (0 = majority, 1 = minority)
% ratio : desired minority / majority ratio after oversampling
%         ratio = 1.0 -> full 1:1 balance
% k     : number of minority nearest neighbors
%
% Example:
% [Xb,Yb] = smote_oversample(Xs,Ys,1.0,5);
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
% Number of synthetic samples required
% -------------------------------------------------------------

desiredMinority = floor(nn * ratio);

G = desiredMinority - np;

if G <= 0

    X_bal = X;
    Y_bal = Y;
    return;

end

%% ============================================================
% Minority-to-minority distances
% ============================================================

Dpos = pdist2( ...
    Xp, ...
    Xp, ...
    'euclidean');

% Prevent self-neighbor
Dpos(1:np+1:end) = inf;

%% ------------------------------------------------------------
% Number of available neighbors
% -------------------------------------------------------------

kUse = min(k,np-1);

[~,neighborOrder] = ...
    sort( ...
        Dpos, ...
        2, ...
        'ascend');

%% ============================================================
% Generate synthetic samples
% ============================================================

synthetic = zeros( ...
    G, ...
    size(X,2));

for g = 1:G

    %% --------------------------------------------------------
    % Random minority seed
    % ---------------------------------------------------------

    i = randi(np);

    %% --------------------------------------------------------
    % Randomly select one of k nearest minority neighbors
    % ---------------------------------------------------------

    neighborList = ...
        neighborOrder(i,1:kUse);

    nbr = ...
        neighborList( ...
            randi(kUse));

    %% --------------------------------------------------------
    % Interpolation coefficient
    % ---------------------------------------------------------

    alpha = rand();

    %% --------------------------------------------------------
    % Standard SMOTE interpolation
    %
    % x_new = x_i + alpha(x_neighbor - x_i)
    % ---------------------------------------------------------

    synthetic(g,:) = ...
        Xp(i,:) + ...
        alpha .* ...
        (Xp(nbr,:) - Xp(i,:));

end

%% ============================================================
% Return augmented dataset
% ============================================================

X_bal = ...
    [X;
     synthetic];

Y_bal = ...
    [Y;
     ones(G,1)];

end