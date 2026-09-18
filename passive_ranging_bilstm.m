clc;
clear;
close all;

%% =====================================================
% BiLSTM Residual Learning
% Uses 20 trajectories generated from RunID
%% =====================================================

T = readtable('passive_ranging_dataset_big.csv');

%% FEATURES

Xall = [ ...
    T.KFX ...
    T.KFY ...
    T.KFZ ...
    T.ConditionNumber ...
    T.Az1 T.El1 ...
    T.Az2 T.El2 ...
    T.Az3 T.El3 ...
    T.Az4 T.El4 ...
    T.Az5 T.El5 ...
    T.Az6 T.El6 ];

%% TARGETS (Residual Learning)

Yall = [ ...
    T.TrueX - T.KFX ...
    T.TrueY - T.KFY ...
    T.TrueZ - T.KFZ ];

%% NORMALIZATION

muX = mean(Xall);
sigX = std(Xall);

Xall = (Xall - muX) ./ sigX;

%% CREATE SEQUENCES

numRuns = max(T.RunID);

XSeq = cell(numRuns,1);
YSeq = cell(numRuns,1);

for run = 1:numRuns

    idx = (T.RunID == run);

    Xrun = Xall(idx,:);
    Yrun = Yall(idx,:);

    XSeq{run} = Xrun';
    YSeq{run} = Yrun';

end

%% TRAIN / TEST SPLIT

trainRuns = 1:16;
testRuns = 17:20;

XTrain = XSeq(trainRuns);
YTrain = YSeq(trainRuns);

XTest = XSeq(testRuns);
YTest = YSeq(testRuns);

%% NETWORK

numFeatures = size(XTrain{1},1);

layers = [

    sequenceInputLayer(numFeatures)

    bilstmLayer(128,'OutputMode','sequence')

    fullyConnectedLayer(64)

    reluLayer

    fullyConnectedLayer(3)

    regressionLayer

];

%% TRAINING OPTIONS

options = trainingOptions('adam', ...
    'MaxEpochs',200, ...
    'MiniBatchSize',4, ...
    'InitialLearnRate',1e-3, ...
    'GradientThreshold',1, ...
    'Shuffle','every-epoch', ...
    'Verbose',true, ...
    'Plots','training-progress');

%% TRAIN NETWORK

net = trainNetwork( ...
    XTrain, ...
    YTrain, ...
    layers, ...
    options);

%% PREDICTION

YPred = predict(net,XTest);

%% EVALUATION

allKF = [];
allTruth = [];
allBi = [];

for r = 1:length(testRuns)

    runID = testRuns(r);

    idx = (T.RunID == runID);

    KF = [ ...
        T.KFX(idx) ...
        T.KFY(idx) ...
        T.KFZ(idx)];

    Truth = [ ...
        T.TrueX(idx) ...
        T.TrueY(idx) ...
        T.TrueZ(idx)];

    ResidualPred = YPred{r}';

    BiPos = KF + ResidualPred;

    allKF = [allKF ; KF];
    allTruth = [allTruth ; Truth];
    allBi = [allBi ; BiPos];

end

%% RMSE

kfRMSE = sqrt(mean(sum((allKF-allTruth).^2,2)));

biRMSE = sqrt(mean(sum((allBi-allTruth).^2,2)));

fprintf('\n====================================\n');
fprintf('KALMAN RMSE = %.3f m\n',kfRMSE);
fprintf('BiLSTM RMSE = %.3f m\n',biRMSE);
fprintf('Improvement = %.2f %%\n', ...
    100*(kfRMSE-biRMSE)/kfRMSE);
fprintf('====================================\n');

%% POSITION ERROR

kfErr = vecnorm(allKF-allTruth,2,2);
biErr = vecnorm(allBi-allTruth,2,2);

%% PLOT 1

figure;

plot(kfErr,'b');
hold on;

plot(biErr,'r');

grid on;

legend('Kalman','BiLSTM');

xlabel('Sample');
ylabel('Error (m)');

title('Error Comparison');

%% PLOT 2

figure;

histogram(kfErr,30);

hold on;

histogram(biErr,30);

legend('Kalman','BiLSTM');

title('Error Distribution');

%% SAVE

save('bilstm_model_v2.mat','net');

disp('Training Complete');