clear;
close all;
clc;

datasets = ["NASA","RELINK","SOFTLAB","AEEEM"];

gammas = 0:0.1:1.0;


seeds = 1:30;

nTrees   = 100;
lambdaLR = 1e-4;

samplingRatio = 1.0;
otomoConf     = 0.995;


outputRoot = fullfile(pwd,'OTOMO_Gamma_Sensitivity');

if ~exist(outputRoot,'dir')
    mkdir(outputRoot);
end

nDatasets = numel(datasets);
nGamma    = numel(gammas);

RF_ROC = nan(nDatasets,nGamma);
RF_G   = nan(nDatasets,nGamma);
RF_F1  = nan(nDatasets,nGamma);
RF_MCC = nan(nDatasets,nGamma);

LR_ROC = nan(nDatasets,nGamma);
LR_G   = nan(nDatasets,nGamma);
LR_F1  = nan(nDatasets,nGamma);
LR_MCC = nan(nDatasets,nGamma);

ValidPairs = zeros(nDatasets,nGamma);

for d = 1:nDatasets

    dataset = datasets(d);

    fprintf('\n');
    fprintf('############################################################\n');
    fprintf('DATASET: %s\n',char(dataset));
    fprintf('############################################################\n');

    files = dir(fullfile('data',char(dataset),'*.xlsx'));

    if isempty(files)
        warning('No files found for dataset %s.',char(dataset));
        continue;
    end

    [~,idx] = sort({files.name});
    files = files(idx);

    nProj = numel(files);

    X = cell(nProj,1);
    Y = cell(nProj,1);
    names = strings(nProj,1);

    for i = 1:nProj

        T = readmatrix( fullfile(files(i).folder,files(i).name));

        Xi = double(T(:,1:end-1));
        Yi = double(T(:,end));
        Yi = Yi(:);

        valid = all(isfinite(Xi),2) & isfinite(Yi);

        Xi = Xi(valid,:);
        Yi = Yi(valid);

        if any(~ismember(unique(Yi),[0 1]))

            error( 'Dataset %s project %s contains labels other than 0/1.', char(dataset), files(i).name);

        end

        X{i} = Xi;
        Y{i} = Yi;

        names(i) = erase(string(files(i).name),".xlsx");

        fprintf( 'Loaded %-10s | N=%4d | Defect=%4d | Non-defect=%4d | Features=%d\n', char(names(i)), size(Xi,1), sum(Yi == 1), sum(Yi == 0), size(Xi,2));

    end

    for g = 1:nGamma

        gamma = gammas(g);

        fprintf('\n');
        fprintf('============================================================\n');
        fprintf('%s | GAMMA = %.1f\n',char(dataset),gamma);
        fprintf('============================================================\n');

        nPairs = nProj * (nProj - 1);

        pairRF_ROC = nan(nPairs,1);
        pairRF_G   = nan(nPairs,1);
        pairRF_F1  = nan(nPairs,1);
        pairRF_MCC = nan(nPairs,1);

        pairLR_ROC = nan(nPairs,1);
        pairLR_G   = nan(nPairs,1);
        pairLR_F1  = nan(nPairs,1);
        pairLR_MCC = nan(nPairs,1);

        pairCount = 0;

        for s = 1:nProj

            for t = 1:nProj

                if s == t
                    continue;
                end

                pairCount = pairCount + 1;

                fprintf( '\n%s | %s -> %s | gamma=%.1f\n', char(dataset), char(names(s)), char(names(t)), gamma);

                Xs = safeZscore(X{s});
                Ys = Y{s};

                Xt = safeZscore(X{t});
                Yt = Y{t};

                rfROC = nan(numel(seeds),1);
                rfG   = nan(numel(seeds),1);
                rfF1  = nan(numel(seeds),1);
                rfMCC = nan(numel(seeds),1);

                lrROC = nan(numel(seeds),1);
                lrG   = nan(numel(seeds),1);
                lrF1  = nan(numel(seeds),1);
                lrMCC = nan(numel(seeds),1);

                for rr = 1:numel(seeds)

                    seed = seeds(rr);

                    try

                        rng(seed,'twister');

                        [Xb,Yb,~,~,~] = otomo( Xs, Ys, Xt, Yt, gamma, samplingRatio, otomoConf,seed);

                        Xtrain = coral(Xb,Xt);

                        Ytrain = double(Yb(:));
                        Xtest  = Xt;

                        Xtrain = real(double(Xtrain));
                        Xtest  = real(double(Xtest));

                        Xtrain(~isfinite(Xtrain)) = 0;
                        Xtest(~isfinite(Xtest))   = 0;

                        if isempty(Xtrain)
                            error('Empty training data.');
                        end

                        if size(Xtrain,1) ~= numel(Ytrain)

                            error( 'Training rows/labels mismatch.');

                        end

                        if size(Xtrain,2) ~= size(Xtest,2)

                            error( 'Train/test feature mismatch.');

                        end

                        if numel(unique(Ytrain)) ~= 2

                            error( 'Training data does not contain two classes.');

                        end

                        rng(seed,'twister');

                        mtry = max( 1, floor( sqrt(size(Xtrain,2))));

                        rf = TreeBagger( nTrees, Xtrain, Ytrain, 'Method','classification', 'NumPredictorsToSample',mtry, 'MinLeafSize',1, 'OOBPrediction','off');

                        [predRF,scoreRF] = predict(rf,Xtest);

                        predRF = str2double(string(predRF));

                        predRF = predRF(:);

                        posRF = find( string(rf.ClassNames) == "1", 1);

                        if isempty(posRF)
                            error('RF positive class not found.');
                        end

                        MR = metrics( Yt, predRF, scoreRF(:,posRF));

                        rfROC(rr) = MR.roc;
                        rfG(rr)   = MR.g;
                        rfF1(rr)  = MR.f1;
                        rfMCC(rr) = MR.mcc;

                        rng(seed,'twister');

                        lr = fitclinear( Xtrain, Ytrain, 'Learner','logistic', 'Regularization','ridge', 'Lambda',lambdaLR, 'Solver','lbfgs', 'ClassNames',[0 1]);

                        [predLR,scoreLR] = predict(lr,Xtest);

                        predLR = double(predLR(:));

                        posLR = find( double(lr.ClassNames) == 1, 1);

                        if isempty(posLR)
                            error('LR positive class not found.');
                        end

                        ML = metrics( Yt, predLR, scoreLR(:,posLR));

                        lrROC(rr) = ML.roc;
                        lrG(rr)   = ML.g;
                        lrF1(rr)  = ML.f1;
                        lrMCC(rr) = ML.mcc;

                    catch ME

                        warning( '%s | %s -> %s | gamma=%.1f | seed=%d failed: %s', char(dataset), char(names(s)), char(names(t)), gamma, seed, ME.message);

                    end

                end

                pairRF_ROC(pairCount) = safeMean(rfROC);
                pairRF_G(pairCount)   = safeMean(rfG);
                pairRF_F1(pairCount)  = safeMean(rfF1);
                pairRF_MCC(pairCount) = safeMean(rfMCC);

                pairLR_ROC(pairCount) = safeMean(lrROC);
                pairLR_G(pairCount)   = safeMean(lrG);
                pairLR_F1(pairCount)  = safeMean(lrF1);
                pairLR_MCC(pairCount) = safeMean(lrMCC);

                fprintf( ['PAIR AVG | RF: AUC=%.4f G=%.4f F1=%.4f MCC=%.4f' ' | LR: AUC=%.4f G=%.4f F1=%.4f MCC=%.4f\n'], pairRF_ROC(pairCount), pairRF_G(pairCount), pairRF_F1(pairCount), pairRF_MCC(pairCount), pairLR_ROC(pairCount), pairLR_G(pairCount), pairLR_F1(pairCount), pairLR_MCC(pairCount));

            end
        end

        RF_ROC(d,g) = safeMean(pairRF_ROC);
        RF_G(d,g)   = safeMean(pairRF_G);
        RF_F1(d,g)  = safeMean(pairRF_F1);
        RF_MCC(d,g) = safeMean(pairRF_MCC);

        LR_ROC(d,g) = safeMean(pairLR_ROC);
        LR_G(d,g)   = safeMean(pairLR_G);
        LR_F1(d,g)  = safeMean(pairLR_F1);
        LR_MCC(d,g) = safeMean(pairLR_MCC);

        validPairMask = isfinite(pairRF_ROC) & isfinite(pairRF_G)   & isfinite(pairRF_F1)  & isfinite(pairRF_MCC);

        ValidPairs(d,g) = sum(validPairMask);

        fprintf('\n');
        fprintf('------------------------------------------------------------\n');
        fprintf('OVERALL DATASET AVERAGE\n');
        fprintf('Dataset = %s | gamma = %.1f\n', char(dataset),gamma);

        fprintf( 'RF | ROC-AUC=%.4f | G=%.4f | F1=%.4f | MCC=%.4f\n', RF_ROC(d,g), RF_G(d,g), RF_F1(d,g), RF_MCC(d,g));

        fprintf( 'LR | ROC-AUC=%.4f | G=%.4f | F1=%.4f | MCC=%.4f\n', LR_ROC(d,g), LR_G(d,g), LR_F1(d,g), LR_MCC(d,g));

        fprintf('------------------------------------------------------------\n');

    end
end

DatasetColumn = strings(nDatasets*nGamma,1);
GammaColumn   = zeros(nDatasets*nGamma,1);

RF_ROC_Column = zeros(nDatasets*nGamma,1);
RF_G_Column   = zeros(nDatasets*nGamma,1);
RF_F1_Column  = zeros(nDatasets*nGamma,1);
RF_MCC_Column = zeros(nDatasets*nGamma,1);

LR_ROC_Column = zeros(nDatasets*nGamma,1);
LR_G_Column   = zeros(nDatasets*nGamma,1);
LR_F1_Column  = zeros(nDatasets*nGamma,1);
LR_MCC_Column = zeros(nDatasets*nGamma,1);

ValidPairColumn = zeros(nDatasets*nGamma,1);

row = 0;

for d = 1:nDatasets

    for g = 1:nGamma

        row = row + 1;

        DatasetColumn(row) = datasets(d);
        GammaColumn(row)   = gammas(g);

        RF_ROC_Column(row) = RF_ROC(d,g);
        RF_G_Column(row)   = RF_G(d,g);
        RF_F1_Column(row)  = RF_F1(d,g);
        RF_MCC_Column(row) = RF_MCC(d,g);

        LR_ROC_Column(row) = LR_ROC(d,g);
        LR_G_Column(row)   = LR_G(d,g);
        LR_F1_Column(row)  = LR_F1(d,g);
        LR_MCC_Column(row) = LR_MCC(d,g);

        ValidPairColumn(row) = ValidPairs(d,g);

    end
end

Summary = table( DatasetColumn, GammaColumn, ValidPairColumn, RF_ROC_Column, RF_G_Column, RF_F1_Column, RF_MCC_Column, LR_ROC_Column, LR_G_Column, LR_F1_Column, LR_MCC_Column, 'VariableNames', { 'Dataset', 'Gamma', 'Valid_Pairs', 'RF_ROC_AUC', 'RF_G_Mean', 'RF_F1', 'RF_MCC', 'LR_ROC_AUC', 'LR_G_Mean', 'LR_F1', 'LR_MCC'});

summaryFile = fullfile( outputRoot, 'OTOMO_Gamma_Sensitivity_Summary.csv');

writetable(Summary,summaryFile);

save( fullfile(outputRoot,'OTOMO_Gamma_Sensitivity.mat'), 'datasets', 'gammas', 'RF_ROC', 'RF_G', 'RF_F1', 'RF_MCC', 'LR_ROC', 'LR_G', 'LR_F1', 'LR_MCC', 'ValidPairs');

% Draw every figure from the saved CSV, exactly as sensitivity_figure_from_csv.m does.
figureFolder = fullfile(outputRoot,'Sensitivity_Figures');
plot_sensitivity_figures(summaryFile,figureFolder,datasets);

fprintf('\n');
fprintf('============================================================\n');
fprintf('GAMMA SENSITIVITY EXPERIMENT COMPLETE\n');
fprintf('============================================================\n');
fprintf('Summary CSV:\n%s\n',summaryFile);
fprintf('\nFigures saved in:\n%s\n',figureFolder);
fprintf('============================================================\n');

function Z = safeZscore(X)

X = double(X);

mu = mean(X,1,'omitnan');
sd = std(X,0,1,'omitnan');

sd(~isfinite(sd) | sd < eps) = 1;

Z = (X - mu) ./ sd;

Z(~isfinite(Z)) = 0;

end

function M = metrics(y,pred,score)

y     = double(y(:));
pred  = double(pred(:));
score = double(score(:));

TP = sum(pred == 1 & y == 1);
TN = sum(pred == 0 & y == 0);
FP = sum(pred == 1 & y == 0);
FN = sum(pred == 0 & y == 1);

if TP + FP == 0
    precision = 0;
else
    precision = TP / (TP + FP);
end

if TP + FN == 0
    recall = 0;
else
    recall = TP / (TP + FN);
end

if TN + FP == 0
    specificity = 0;
else
    specificity = TN / (TN + FP);
end

if precision + recall == 0
    M.f1 = 0;
else

    M.f1 = 2 * precision * recall / (precision + recall);

end

if recall + specificity == 0
    M.g = 0;
else
    M.g = sqrt(recall * specificity);
end

denominator = sqrt( double( (TP + FP) * (TP + FN) * (TN + FP) * (TN + FN)));

if denominator == 0

    M.mcc = 0;

else

    M.mcc = (TP * TN - FP * FN) / denominator;

end

M.roc = NaN;

if numel(unique(y)) == 2 && all(isfinite(score))

    if numel(unique(score)) == 1

        M.roc = 0.5;

    else

        [~,~,~,M.roc] = perfcurve(y,score,1);

    end

end

end

function value = safeMean(values)

values = double(values(:));

values = values(isfinite(values));

if isempty(values)

    value = NaN;

else

    value = mean(values);

end

end
