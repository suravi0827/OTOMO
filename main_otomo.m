clear;
close all;
clc;

dataset = "SOFTLAB";
% gamma = 0.4;
gamma = 0;

seeds = 1:30;
nTrees = 100;
lambdaLR = 1e-4;
samplingRatio = 1.0;
otomoConf = 0.995;

outputRoot = fullfile(pwd,'OTOMO_Results');

if ~exist(outputRoot,'dir')
    mkdir(outputRoot);
end

files = dir(fullfile('data',char(dataset),'*.xlsx'));

if isempty(files)
    error('No files found for dataset %s.',char(dataset));
end

[~,idx] = sort({files.name});
files = files(idx);
nProj = numel(files);

X = cell(nProj,1);
Y = cell(nProj,1);
names = strings(nProj,1);

fprintf('\n====================================================\n');
fprintf('DATASET: %s | GAMMA: %.2f\n',char(dataset),gamma);
fprintf('====================================================\n');

for i = 1:nProj

    T = readmatrix(fullfile(files(i).folder,files(i).name));
    Xi = double(T(:,1:end-1));
    Yi = double(T(:,end));
    Yi = Yi(:);

    valid = all(isfinite(Xi),2) & isfinite(Yi);
    Xi = Xi(valid,:);
    Yi = Yi(valid);

    if any(~ismember(unique(Yi),[0 1]))
        error('Dataset %s project %s contains labels other than 0 and 1.',char(dataset),files(i).name);
    end

    X{i} = Xi;
    Y{i} = Yi;
    names(i) = erase(string(files(i).name),".xlsx");

    fprintf('Loaded %-12s | N=%d | Defect=%d | Non-defect=%d | Features=%d\n',char(names(i)),size(Xi,1),sum(Yi == 1),sum(Yi == 0),size(Xi,2));

end

Source = strings(0,1);
Target = strings(0,1);
Gamma = zeros(0,1);
Valid_Runs = zeros(0,1);

Purity = strings(0,1);
ARI_Value = strings(0,1);
NMI_Value = strings(0,1);

Purity_Mean = zeros(0,1);
Purity_SD = zeros(0,1);
ARI_Mean = zeros(0,1);
ARI_SD = zeros(0,1);
NMI_Mean = zeros(0,1);
NMI_SD = zeros(0,1);

RF_ROC_AUC = strings(0,1);
RF_PR_AUC = strings(0,1);
RF_G_Measure = strings(0,1);
RF_F1_Score = strings(0,1);
RF_MCC = strings(0,1);

LR_ROC_AUC = strings(0,1);
LR_PR_AUC = strings(0,1);
LR_G_Measure = strings(0,1);
LR_F1_Score = strings(0,1);
LR_MCC = strings(0,1);

RF_ROC_Mean = zeros(0,1);
RF_ROC_SD = zeros(0,1);
RF_PR_Mean = zeros(0,1);
RF_PR_SD = zeros(0,1);
RF_G_Mean = zeros(0,1);
RF_G_SD = zeros(0,1);
RF_F1_Mean = zeros(0,1);
RF_F1_SD = zeros(0,1);
RF_MCC_Mean = zeros(0,1);
RF_MCC_SD = zeros(0,1);

LR_ROC_Mean = zeros(0,1);
LR_ROC_SD = zeros(0,1);
LR_PR_Mean = zeros(0,1);
LR_PR_SD = zeros(0,1);
LR_G_Mean = zeros(0,1);
LR_G_SD = zeros(0,1);
LR_F1_Mean = zeros(0,1);
LR_F1_SD = zeros(0,1);
LR_MCC_Mean = zeros(0,1);
LR_MCC_SD = zeros(0,1);

row = 0;

for s = 1:nProj

    for t = 1:nProj

        if s == t
            continue;
        end

        Xs = safeZscore(X{s});
        Ys = Y{s};
        Xt = safeZscore(X{t});
        Yt = Y{t};

        fprintf('\n====================================================\n');
        fprintf('OTOMO | %s -> %s | gamma=%.2f\n',char(names(s)),char(names(t)),gamma);
        fprintf('====================================================\n');

        rfROC = nan(numel(seeds),1);
        rfPR = nan(numel(seeds),1);
        rfG = nan(numel(seeds),1);
        rfF1 = nan(numel(seeds),1);
        rfMCC = nan(numel(seeds),1);

        lrROC = nan(numel(seeds),1);
        lrPR = nan(numel(seeds),1);
        lrG = nan(numel(seeds),1);
        lrF1 = nan(numel(seeds),1);
        lrMCC = nan(numel(seeds),1);

        purityValues = nan(numel(seeds),1);
        ariValues = nan(numel(seeds),1);
        nmiValues = nan(numel(seeds),1);

        for rr = 1:numel(seeds)

            seed = seeds(rr);

            try

                rng(seed,'twister');

                [Xb,Yb,purityValue,ariValue,nmiValue] = otomo(Xs,Ys,Xt,Yt,gamma,samplingRatio,otomoConf);

                purityValues(rr) = purityValue;
                ariValues(rr) = ariValue;
                nmiValues(rr) = nmiValue;

              
                Xtrain = coral(Xb,Xt);
                Ytrain = double(Yb(:));
                Xtest = Xt;

                Xtrain = real(double(Xtrain));
                Xtest = real(double(Xtest));

                Xtrain(~isfinite(Xtrain)) = 0;
                Xtest(~isfinite(Xtest)) = 0;

                if isempty(Xtrain)
                    error('Empty training data.');
                end

                if size(Xtrain,1) ~= numel(Ytrain)
                    error('Training rows/labels mismatch: %d vs %d.',size(Xtrain,1),numel(Ytrain));
                end

                if size(Xtrain,2) ~= size(Xtest,2)
                    error('Feature mismatch: train=%d, test=%d.',size(Xtrain,2),size(Xtest,2));
                end

                if numel(unique(Ytrain)) ~= 2
                    error('Training data does not contain two classes.');
                end

                rng(seed,'twister');

                mtry = max(1,floor(sqrt(size(Xtrain,2))));
                rf = TreeBagger(nTrees,Xtrain,Ytrain,'Method','classification','NumPredictorsToSample',mtry,'MinLeafSize',1,'OOBPrediction','off');

                [predRF,scoreRF] = predict(rf,Xtest);

                predRF = str2double(string(predRF));
                predRF = predRF(:);

                posRF = find(string(rf.ClassNames) == "1",1);

                if isempty(posRF)
                    error('RF positive class 1 not found.');
                end

                MR = metrics(Yt,predRF,scoreRF(:,posRF));

                rfROC(rr) = MR.roc;
                rfPR(rr) = MR.pr;
                rfG(rr) = MR.g;
                rfF1(rr) = MR.f1;
                rfMCC(rr) = MR.mcc;

                rng(seed,'twister');

                lr = fitclinear(Xtrain,Ytrain,'Learner','logistic','Regularization','ridge','Lambda',lambdaLR,'Solver','lbfgs','ClassNames',[0 1]);

                [predLR,scoreLR] = predict(lr,Xtest);

                predLR = double(predLR(:));

                posLR = find(double(lr.ClassNames) == 1,1);

                if isempty(posLR)
                    error('LR positive class 1 not found.');
                end

                ML = metrics(Yt,predLR,scoreLR(:,posLR));

                lrROC(rr) = ML.roc;
                lrPR(rr) = ML.pr;
                lrG(rr) = ML.g;
                lrF1(rr) = ML.f1;
                lrMCC(rr) = ML.mcc;

                fprintf('Seed=%2d | RF AUC=%.4f G=%.4f F1=%.4f MCC=%.4f | LR AUC=%.4f G=%.4f F1=%.4f MCC=%.4f\n',seed,MR.roc,MR.g,MR.f1,MR.mcc,ML.roc,ML.g,ML.f1,ML.mcc);

            catch ME

                warning('OTOMO | %s -> %s | seed=%d failed: %s',char(names(s)),char(names(t)),seed,ME.message);

            end

        end

        [purityMean,puritySD] = calculateMeanSD(purityValues);
        [ariMean,ariSD] = calculateMeanSD(ariValues);
        [nmiMean,nmiSD] = calculateMeanSD(nmiValues);

        [rfROCMean,rfROCSD] = calculateMeanSD(rfROC);
        [rfPRMean,rfPRSD] = calculateMeanSD(rfPR);
        [rfGMean,rfGSD] = calculateMeanSD(rfG);
        [rfF1Mean,rfF1SD] = calculateMeanSD(rfF1);
        [rfMCCMean,rfMCCSD] = calculateMeanSD(rfMCC);

        [lrROCMean,lrROCSD] = calculateMeanSD(lrROC);
        [lrPRMean,lrPRSD] = calculateMeanSD(lrPR);
        [lrGMean,lrGSD] = calculateMeanSD(lrG);
        [lrF1Mean,lrF1SD] = calculateMeanSD(lrF1);
        [lrMCCMean,lrMCCSD] = calculateMeanSD(lrMCC);

        validMask = isfinite(rfROC) & isfinite(rfG) & isfinite(rfF1) & isfinite(rfMCC) & isfinite(lrROC) & isfinite(lrG) & isfinite(lrF1) & isfinite(lrMCC);

        row = row + 1;

        Source(row,1) = names(s);
        Target(row,1) = names(t);
        Gamma(row,1) = gamma;
        Valid_Runs(row,1) = sum(validMask);

        Purity(row,1) = formatMeanSD(purityMean,puritySD);
        ARI_Value(row,1) = formatMeanSD(ariMean,ariSD);
        NMI_Value(row,1) = formatMeanSD(nmiMean,nmiSD);

        Purity_Mean(row,1) = purityMean;
        Purity_SD(row,1) = puritySD;
        ARI_Mean(row,1) = ariMean;
        ARI_SD(row,1) = ariSD;
        NMI_Mean(row,1) = nmiMean;
        NMI_SD(row,1) = nmiSD;

        RF_ROC_AUC(row,1) = formatMeanSD(rfROCMean,rfROCSD);
        RF_PR_AUC(row,1) = formatMeanSD(rfPRMean,rfPRSD);
        RF_G_Measure(row,1) = formatMeanSD(rfGMean,rfGSD);
        RF_F1_Score(row,1) = formatMeanSD(rfF1Mean,rfF1SD);
        RF_MCC(row,1) = formatMeanSD(rfMCCMean,rfMCCSD);

        LR_ROC_AUC(row,1) = formatMeanSD(lrROCMean,lrROCSD);
        LR_PR_AUC(row,1) = formatMeanSD(lrPRMean,lrPRSD);
        LR_G_Measure(row,1) = formatMeanSD(lrGMean,lrGSD);
        LR_F1_Score(row,1) = formatMeanSD(lrF1Mean,lrF1SD);
        LR_MCC(row,1) = formatMeanSD(lrMCCMean,lrMCCSD);

        RF_ROC_Mean(row,1) = rfROCMean;
        RF_ROC_SD(row,1) = rfROCSD;
        RF_PR_Mean(row,1) = rfPRMean;
        RF_PR_SD(row,1) = rfPRSD;
        RF_G_Mean(row,1) = rfGMean;
        RF_G_SD(row,1) = rfGSD;
        RF_F1_Mean(row,1) = rfF1Mean;
        RF_F1_SD(row,1) = rfF1SD;
        RF_MCC_Mean(row,1) = rfMCCMean;
        RF_MCC_SD(row,1) = rfMCCSD;

        LR_ROC_Mean(row,1) = lrROCMean;
        LR_ROC_SD(row,1) = lrROCSD;
        LR_PR_Mean(row,1) = lrPRMean;
        LR_PR_SD(row,1) = lrPRSD;
        LR_G_Mean(row,1) = lrGMean;
        LR_G_SD(row,1) = lrGSD;
        LR_F1_Mean(row,1) = lrF1Mean;
        LR_F1_SD(row,1) = lrF1SD;
        LR_MCC_Mean(row,1) = lrMCCMean;
        LR_MCC_SD(row,1) = lrMCCSD;

        fprintf('AVERAGE | %s -> %s | RF AUC=%.4f G=%.4f F1=%.4f MCC=%.4f | LR AUC=%.4f G=%.4f F1=%.4f MCC=%.4f\n',char(names(s)),char(names(t)),rfROCMean,rfGMean,rfF1Mean,rfMCCMean,lrROCMean,lrGMean,lrF1Mean,lrMCCMean);

    end

end

Result = table(Source,Target,Gamma,Valid_Runs,Purity,Purity_Mean,Purity_SD,ARI_Value,ARI_Mean,ARI_SD,NMI_Value,NMI_Mean,NMI_SD,RF_ROC_AUC,RF_PR_AUC,RF_G_Measure,RF_F1_Score,RF_MCC,LR_ROC_AUC,LR_PR_AUC,LR_G_Measure,LR_F1_Score,LR_MCC,RF_ROC_Mean,RF_ROC_SD,RF_PR_Mean,RF_PR_SD,RF_G_Mean,RF_G_SD,RF_F1_Mean,RF_F1_SD,RF_MCC_Mean,RF_MCC_SD,LR_ROC_Mean,LR_ROC_SD,LR_PR_Mean,LR_PR_SD,LR_G_Mean,LR_G_SD,LR_F1_Mean,LR_F1_SD,LR_MCC_Mean,LR_MCC_SD);

outputFile = fullfile(outputRoot,sprintf('%s_OTOMO_Gamma_%.2f_Results.csv',char(dataset),gamma));

writetable(Result,outputFile);

fprintf('\n====================================================\n');
fprintf('EXPERIMENT COMPLETE\n');
fprintf('Dataset: %s\n',char(dataset));
fprintf('Gamma: %.2f\n',gamma);
fprintf('Results: %s\n',outputFile);
fprintf('====================================================\n');


function Z = safeZscore(X)

X = double(X);
mu = mean(X,1,'omitnan');
sd = std(X,0,1,'omitnan');
sd(~isfinite(sd) | sd < eps) = 1;
Z = (X - mu) ./ sd;
Z(~isfinite(Z)) = 0;

end


function M = metrics(y,pred,score)

y = double(y(:));
pred = double(pred(:));
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

M.g = sqrt(recall * specificity);

denominator = sqrt(double((TP + FP) * (TP + FN) * (TN + FP) * (TN + FN)));

if denominator == 0
    M.mcc = 0;
else
    M.mcc = (TP * TN - FP * FN) / denominator;
end

if numel(unique(y)) == 2 && numel(unique(score)) > 1
    try
        [~,~,~,M.roc] = perfcurve(y,score,1);
    catch
        M.roc = NaN;
    end
else
    M.roc = NaN;
end

if numel(unique(y)) == 2 && numel(unique(score)) > 1
    try
        [~,~,~,M.pr] = perfcurve(y,score,1,'XCrit','reca','YCrit','prec');
    catch
        M.pr = NaN;
    end
else
    M.pr = NaN;
end

end


function [mu,sd] = calculateMeanSD(values)

values = double(values(:));
values = values(isfinite(values));

if isempty(values)
    mu = NaN;
    sd = NaN;
    return;
end

mu = mean(values);

if numel(values) == 1
    sd = 0;
else
    sd = std(values,0);
end

end


function txt = formatMeanSD(mu,sd)

if ~isfinite(mu) || ~isfinite(sd)
    txt = "NaN";
else
    txt = sprintf('%.4f +/- %.4f',mu,sd);
end

end