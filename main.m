clear;
close all;
clc;

datasets = ["NASA","RELINK","SOFTLAB","AEEEM"];
gammas = 0.1;
seeds = 1:30;
nTrees = 100;
lambdaLR = 1e-4;

methods = ["OTOMO","CORAL","TCA+","ADASYN","DPP","SMOTE","TOMO"];
% methods = ["OTOMO"]

methods = string(methods(:).');

assert(~isempty(methods),'Select at least one method.');
assert(numel(unique(methods)) == numel(methods),'Method names must be unique.');

nMethods = numel(methods);

samplingRatio = 1.0;
kSMOTE = 7;
kADASYN = 7;
tomoLambda = 0.1;
otomoBeta = 0.1;
otomoConf = 0.995;

outputRoot = fullfile(pwd,'Results');

if ~exist(outputRoot,'dir')
    mkdir(outputRoot);
end

ROC_all = cell(1,numel(datasets));
F1_all = cell(size(ROC_all));
G_all = cell(size(ROC_all));
MCC_all = cell(size(ROC_all));
names_all = cell(size(ROC_all));

completed = false(size(datasets));

for dd = 1:numel(datasets)

    dataset = datasets(dd);

    fprintf('\n');
    fprintf('====================================================\n');
    fprintf('DATASET: %s\n',char(dataset));
    fprintf('====================================================\n');

    files = dir(fullfile('data',char(dataset),'*.xlsx'));

    if isempty(files)
        warning('No files found for %s.',char(dataset));
        continue;
    end

    if numel(files) < 2
        warning('At least two projects are required for %s; skipping dataset.',char(dataset));
        continue;
    end

    [~,idx] = sort({files.name});
    files = files(idx);

    nProj = numel(files);

    X = cell(nProj,1);
    Y = cell(nProj,1);
    names = strings(nProj,1);

    for i = 1:nProj

        T = readmatrix(fullfile(files(i).folder,files(i).name));

        Xi = double(T(:,1:end-1));
        Yi = double(T(:,end));
        Yi = Yi(:);

        valid = all(isfinite(Xi),2) & isfinite(Yi);

        Xi = Xi(valid,:);
        Yi = Yi(valid);

        X{i} = Xi;
        Y{i} = Yi;

        names(i) = erase(string(files(i).name),".xlsx");

        fprintf('Loaded %-12s | N=%d | Defect=%d | Features=%d\n',char(names(i)),size(Xi,1),sum(Yi == 1),size(Xi,2));

    end

    nGamma = numel(gammas);
    nSeed = numel(seeds);

    ROC = nan(nMethods,nGamma,nProj,nProj,nSeed,2);
    F1 = nan(nMethods,nGamma,nProj,nProj,nSeed,2);
    G = nan(nMethods,nGamma,nProj,nProj,nSeed,2);
    MCC = nan(nMethods,nGamma,nProj,nProj,nSeed,2);

    for gg = 1:nGamma

        gamma = gammas(gg);

        fprintf('\n');
        fprintf('====================================================\n');
        fprintf('GAMMA = %.2f\n',gamma);
        fprintf('====================================================\n');

        for mm = 1:nMethods

            method = methods(mm);

            fprintf('\n');
            fprintf('----------------------------------------------------\n');
            fprintf('METHOD: %s\n',char(method));
            fprintf('----------------------------------------------------\n');

            for s = 1:nProj

                for t = 1:nProj

                    if s == t
                        continue;
                    end

                    Xs = X{s};
                    Ys = Y{s};

                    Xt = X{t};
                    Yt = Y{t};

                    Xs = safeZscore(Xs);
                    Xt = safeZscore(Xt);

                    fprintf('%-8s | %s -> %s\n',char(method),char(names(s)),char(names(t)));

                    for rr = 1:nSeed

                        seed = seeds(rr);

                        try

                            rng(seed,'twister');

                            [Xtrain,Ytrain,Xtest] = applyMethod(method,Xs,Ys,Xt,Yt,gamma,seed,samplingRatio,kSMOTE,kADASYN,tomoLambda,otomoConf);

                            Xtrain = real(double(Xtrain));
                            Ytrain = double(Ytrain(:));
                            Xtest = real(double(Xtest));

                            if isempty(Xtrain)
                                error('%s returned empty training data.',char(method));
                            end

                            if size(Xtrain,1) ~= numel(Ytrain)
                                error('Training data/label mismatch: %d rows vs %d labels.',size(Xtrain,1),numel(Ytrain));
                            end

                            if size(Xtrain,2) ~= size(Xtest,2)
                                error('Feature mismatch after %s: train=%d, test=%d.',char(method),size(Xtrain,2),size(Xtest,2));
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

                            M = metrics(Yt,predRF,scoreRF(:,posRF));

                            ROC(mm,gg,s,t,rr,1) = M.roc;
                            F1(mm,gg,s,t,rr,1) = M.f1;
                            G(mm,gg,s,t,rr,1) = M.g;
                            MCC(mm,gg,s,t,rr,1) = M.mcc;

                            rng(seed,'twister');

                            lr = fitclinear(Xtrain,Ytrain,'Learner','logistic','Regularization','ridge','Lambda',lambdaLR,'Solver','lbfgs','ClassNames',[0 1]);

                            [predLR,scoreLR] = predict(lr,Xtest);

                            predLR = double(predLR(:));

                            posLR = find(double(lr.ClassNames) == 1,1);

                            if isempty(posLR)
                                error('LR positive class 1 not found.');
                            end

                            M = metrics(Yt,predLR,scoreLR(:,posLR));

                            ROC(mm,gg,s,t,rr,2) = M.roc;
                            F1(mm,gg,s,t,rr,2) = M.f1;
                            G(mm,gg,s,t,rr,2) = M.g;
                            MCC(mm,gg,s,t,rr,2) = M.mcc;

                        catch ME

                            warning('%s | %s -> %s | gamma=%.2f | seed=%d failed: %s',char(method),char(names(s)),char(names(t)),gamma,seed,ME.message);

                        end

                    end

                end

            end

        end

    end

    ROC_all{dd} = ROC;
    F1_all{dd} = F1;
    G_all{dd} = G;
    MCC_all{dd} = MCC;
    names_all{dd} = names;

    completed(dd) = true;

    fprintf('\nCompleted %s; results retained for combined summaries.\n',char(dataset));

end

summaryOrder = ["NASA","RELINK","SOFTLAB","AEEEM"];

[found,summaryIdx] = ismember(summaryOrder,datasets);

summaryIdx = summaryIdx(found);

summaryIdx = [summaryIdx,setdiff(1:numel(datasets),summaryIdx,'stable')];

summaryIdx = summaryIdx(completed(summaryIdx));

if isempty(summaryIdx)

    warning('No datasets were processed; combined summaries were not created.');

else

    summaryFolder = fullfile(outputRoot,'Metric_Summaries');

    if ~any(methods == "OTOMO")
        fprintf('OTOMO is not selected; summaries include averages without W/T/L rows.\n');
    end

    writeMethodSummary(ROC_all(summaryIdx),F1_all(summaryIdx),G_all(summaryIdx),MCC_all(summaryIdx),gammas,names_all(summaryIdx),methods,1,fullfile(summaryFolder,'Random_Forest_All_Datasets.csv'),datasets(summaryIdx));

    writeMethodSummary(ROC_all(summaryIdx),F1_all(summaryIdx),G_all(summaryIdx),MCC_all(summaryIdx),gammas,names_all(summaryIdx),methods,2,fullfile(summaryFolder,'Logistic_Regression_All_Datasets.csv'),datasets(summaryIdx));

    if any(~completed)
        warning('Combined summaries exclude skipped datasets: %s',char(strjoin(datasets(~completed),', ')));
    end

end

fprintf('\n');
fprintf('====================================================\n');
fprintf('ALL EXPERIMENTS COMPLETE\n');
fprintf('====================================================\n');


function [Xtrain,Ytrain,Xtest] = applyMethod(method,Xs,Ys,Xt,Yt,gamma,seed,samplingRatio,kSMOTE,kADASYN,tomoLambda,otomoConf)

rng(seed,'twister');

method = upper(string(method));

switch method

    case "CORAL"

        Xtrain = coral(Xs,Xt);
        Ytrain = Ys;
        Xtest = Xt;

    case "TCA+"

        [Xtrain,Xtest] = tca_plus(Xs,Xt);
        Ytrain = Ys;

    case "SMOTE"

        rng(seed,'twister');

        [Xb,Yb] = smote_oversample(Xs,Ys,samplingRatio,kSMOTE);

        Xtrain = coral(Xb,Xt);
        Ytrain = Yb;
        Xtest = Xt;

    case "ADASYN"

        rng(seed,'twister');

        [Xb,Yb] = adasyn_oversample(Xs,Ys,samplingRatio,kADASYN);

        Xtrain = coral(Xb,Xt);
        Ytrain = Yb;
        Xtest = Xt;

    case "DPP"

        rng(seed,'twister');

        [Xb,Yb] = dpp_oversample(Xs,Ys,samplingRatio,0.1);

        Xtrain = coral(Xb,Xt);
        Ytrain = Yb;
        Xtest = Xt;

    case "TOMO"

        rng(seed,'twister');

        [Xb,Yb] = tomo(Xs,Ys,Xt,tomoLambda,samplingRatio,seed);

        Xtrain = coral(Xb,Xt);
        Ytrain = Yb;
        Xtest = Xt;

    case "OTOMO"

        rng(seed,'twister');

        [Xb,Yb,purity,ari,nmi] = otomo(Xs,Ys,Xt,Yt,gamma,samplingRatio,otomoConf,seed);

        Xtrain = coral(Xb,Xt);
        Ytrain = Yb;
        Xtest = Xt;

    otherwise

        error('Unknown method: %s',char(method));

end

Xtrain = real(double(Xtrain));
Ytrain = double(Ytrain(:));
Xtest = real(double(Xtest));

Xtrain(~isfinite(Xtrain)) = 0;
Xtest(~isfinite(Xtest)) = 0;

end


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

if recall + specificity == 0
    M.g = 0;
else
    M.g = sqrt(recall * specificity);
end

denominator = sqrt(double((TP + FP) * (TP + FN) * (TN + FP) * (TN + FN)));

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


function writeMethodSummary(ROC,F1,G,MCC,gammas,names,methods,classifierIdx,file,datasetNames)

if ~iscell(ROC)

    ROC = {ROC};
    F1 = {F1};
    G = {G};
    MCC = {MCC};

    names = {names};

end

nDatasets = numel(ROC);

if nargin < 10
    datasetNames = "Dataset " + string(1:nDatasets);
end

datasetNames = string(datasetNames);
methods = string(methods(:).');

preferred = ["OTOMO","CORAL","TCA+","DPP","ADASYN","SMOTE","TOMO"];

[present,order] = ismember(preferred,methods);

order = [order(present),setdiff(1:numel(methods),order(present),'stable')];

methodNames = methods(order);
nMethods = numel(methodNames);

reference = find(methodNames == "OTOMO",1);

assert(nMethods > 0,'Select at least one method.');
assert(numel(unique(methods)) == numel(methods),'Method names must be unique.');
assert(numel(datasetNames) == nDatasets,'Provide one name per dataset.');
assert(iscell(names) && numel(names) == nDatasets,'Provide one project-name vector per dataset.');

metricData = {ROC,F1,G,MCC};

metricNames = ["ROC_AUC","F1_Score","G_Measure","MCC_Score"];

[folder,base,~] = fileparts(file);

if isempty(folder)
    folder = pwd;
end

if ~exist(folder,'dir')
    mkdir(folder);
end

for kk = 1:numel(metricNames)

    data = metricData{kk};

    assert(numel(data) == nDatasets,'Metric dataset counts must match.');

    for gg = 1:numel(gammas)

        result = [{'Source','Target'},cellstr(methodNames)];

        allMeans = zeros(0,nMethods);

        for dd = 1:nDatasets

            A = data{dd};

            projectNames = string(names{dd});
            nProj = numel(projectNames);

            assert(size(A,1) == numel(methods) && size(A,2) == numel(gammas),'Method/gamma dimensions do not match.');

            assert(size(A,3) == nProj && size(A,4) == nProj,'Project dimensions do not match.');

            assert(classifierIdx >= 1 && classifierIdx <= size(A,6) && classifierIdx == floor(classifierIdx),'Invalid classifier index.');

            result(end+1,:) = [{char(datasetNames(dd))},repmat({''},1,nMethods + 1)];

            means = nan(nProj * (nProj - 1),nMethods);
            counts = zeros(size(means));

            rr = 0;

            for s = 1:nProj

                for t = 1:nProj

                    if s == t
                        continue;
                    end

                    rr = rr + 1;

                    entry = cell(1,nMethods + 2);

                    entry{1} = char(projectNames(s));
                    entry{2} = char(projectNames(t));

                    for jj = 1:nMethods

                        values = reshape(A(order(jj),gg,s,t,:,classifierIdx),[],1);

                        values = values(isfinite(values));

                        counts(rr,jj) = numel(values);

                        if isempty(values)

                            entry{jj+2} = 'NaN';

                        else

                            means(rr,jj) = mean(values);

                            entry{jj+2} = sprintf('%.4f +/- %.4f',mean(values),std(values,0));

                        end

                    end

                    result(end+1,:) = entry;

                end

            end

            if any(counts(:) < size(A,5))

                warning('writeMethodSummary:MissingRuns','%s, %s, gamma=%g: some tasks have missing runs; means use available finite runs.',char(datasetNames(dd)),char(metricNames(kk)),gammas(gg));

            end

            result = [result; summaryRows(means,reference,false)];

            allMeans = [allMeans; means];

        end

        result = [result; summaryRows(allMeans,reference,true)];

        outputFile = fullfile(folder,sprintf('%s_%s_gamma_%03d_%g.csv',base,char(metricNames(kk)),gg,gammas(gg)));

        writecell(result,outputFile);

        fprintf('Wrote: %s\n',outputFile);

    end

end

end


function rows = summaryRows(means,reference,isOverall)

nMethods = size(means,2);

hasReference = ~isempty(reference);

averageRow = 1 + double(hasReference);

rows = repmat({''},averageRow,nMethods + 2);

if isOverall

    rows{averageRow,1} = 'Overall Average (%)';

    if hasReference
        rows{1,1} = 'Overall W/T/L';
    end

else

    rows{averageRow,1} = 'Average (%)';

    if hasReference
        rows{1,1} = 'W/T/L';
    end

end

roundedMeans = round(means,4);

for jj = 1:nMethods

    if hasReference

        if jj == reference

            rows{1,jj+2} = '';

        else

            valid = isfinite(roundedMeans(:,reference)) & isfinite(roundedMeans(:,jj));

            difference = roundedMeans(valid,reference) - roundedMeans(valid,jj);

            if isempty(difference)

                rows{1,jj+2} = 'N/A';

            else

                rows{1,jj+2} = sprintf('%d/%d/%d',sum(difference > 0),sum(difference == 0),sum(difference < 0));

            end

        end

    end

    values = means(:,jj);

    values = values(isfinite(values));

    if isempty(values)
        rows{averageRow,jj+2} = 'NaN';
    else
        rows{averageRow,jj+2} = sprintf('%.2f',100 * mean(values));
    end

end

if any(~isfinite(means(:)))

    coverage = [{'Valid tasks (average)',''},repmat({''},1,nMethods)];

    for jj = 1:nMethods

        coverage{jj+2} = sprintf('%d/%d',sum(isfinite(means(:,jj))),size(means,1));

    end

    rows = [rows; coverage];

    if hasReference

        comparisons = [{'Valid comparisons vs OTOMO',''},repmat({''},1,nMethods)];

        for jj = 1:nMethods

            comparisons{jj+2} = sprintf('%d/%d',sum(isfinite(means(:,reference)) & isfinite(means(:,jj))),size(means,1));

        end

        comparisons{reference+2} = '';

        rows = [rows; comparisons];

    end

end

end