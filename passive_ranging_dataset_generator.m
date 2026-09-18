clc;
clear;
close all;

numRuns = 20;
allData = table();

for run = 1:numRuns

    rng(run);

    fprintf('\nRunning Scenario %d/%d\n',run,numRuns);

%% =====================================================
% TERRESTRIAL PASSIVE RANGING PROJECT
% MOVING SENSORS + TRIANGULATION + KALMAN FILTER
% FULL RESEARCH BASELINE
%% =====================================================

N = 1000;
dt = 1;
numSensors = 6;

%% TARGET MOTION

targetPos = zeros(N,3);
targetVel = zeros(N,3);

targetPos(1,:) = [1000 1000 20];
targetVel(1,:) = [2 1 0];

accStd = 0.03 + 0.05*rand;

for k = 2:N

    if mod(k,200)==0
        targetVel(k-1,:) = targetVel(k-1,:) + [0.5 -0.3 0];
    end

    a = accStd*randn(1,3);

    targetVel(k,:) = targetVel(k-1,:) + a*dt;
    targetPos(k,:) = targetPos(k-1,:) + targetVel(k,:)*dt;

end

%% MOVING SENSOR TRAJECTORIES

sensorTrue = zeros(N,3,numSensors);

Rmove = 1000 + 500*rand;

for k = 1:N

    ang = 0.004*k;

    sensorTrue(k,:,1) = [2000 + Rmove*cos(ang), ...
                         2000 + Rmove*sin(ang), ...
                         20];

    sensorTrue(k,:,2) = [2000 - Rmove*cos(ang), ...
                         2000 - Rmove*sin(ang), ...
                         20];

    sensorTrue(k,:,3) = [1000 + 300*sin(0.01*k), ...
                         3500, ...
                         30];

    sensorTrue(k,:,4) = [3500, ...
                         1000 + 300*cos(0.01*k), ...
                         25];
    sensorTrue(k,:,5) = [1500*cos(0.008*k) ...
        1500*sin(0.008*k) ...
        100];

    sensorTrue(k,:,6) = [-1500*cos(0.008*k) ...
        -1500*sin(0.008*k) ...
        150];
end

%% SENSOR POSITION ERROR

gpsSigma = 0.5 + rand;

sensorPos = sensorTrue + gpsSigma*randn(size(sensorTrue));

%% ANGLE NOISE + BIAS

azNoise = deg2rad(0.2 + 0.3*rand);
elNoise = deg2rad(0.1 + 0.2*rand);

azBias = deg2rad(0.15)*sin((1:N)'/250);
elBias = deg2rad(0.10)*cos((1:N)'/300);

%% STORAGE

triPos = zeros(N,3);
azStore = zeros(N,numSensors);
elStore = zeros(N,numSensors);
conditionNum = zeros(N,1);

%% TRIANGULATION

for k = 1:N

    A = [];
    b = [];

    tx = targetPos(k,1);
    ty = targetPos(k,2);
    tz = targetPos(k,3);

    for s = 1:numSensors

        sx = sensorPos(k,1,s);
        sy = sensorPos(k,2,s);
        sz = sensorPos(k,3,s);

        dx = tx-sx;
        dy = ty-sy;
        dz = tz-sz;

        az = atan2(dy,dx);
        el = atan2(dz,sqrt(dx^2+dy^2));

        azMeas = az + azNoise*randn + azBias(k);
        elMeas = el + elNoise*randn + elBias(k);

        if rand < 0.01
            azMeas = azMeas + deg2rad(2);
            elMeas = elMeas + deg2rad(1);
        end

        azStore(k,s) = azMeas;
        elStore(k,s) = elMeas;

        ux = cos(elMeas)*cos(azMeas);
        uy = cos(elMeas)*sin(azMeas);
        uz = sin(elMeas);

        u = [ux uy uz];

        Ai = eye(3) - (u'*u);

        w = 1/(azNoise^2 + elNoise^2);

        A = [A ; sqrt(w)*Ai];
        b = [b ; sqrt(w)*Ai*[sx;sy;sz]];

    end

    conditionNum(k) = cond(A);

    if conditionNum(k) > 1e6
        triPos(k,:) = [NaN NaN NaN];
    else
        triPos(k,:) = (A\b)';
    end

end

%% REMOVE BAD GEOMETRY POINTS

bad = isnan(triPos(:,1));

triPos = triPos(~bad,:);
targetValid = targetPos(~bad,:);

%% TRIANGULATION ERROR

errorVec = vecnorm(triPos-targetValid,2,2);
RMSE = sqrt(mean(errorVec.^2));

%% KALMAN FILTER

M = size(triPos,1);

xKF = zeros(6,M);

xKF(:,1) = [triPos(1,1);
            triPos(1,2);
            triPos(1,3);
            0;0;0];

F = [1 0 0 dt 0 0;
     0 1 0 0 dt 0;
     0 0 1 0 0 dt;
     0 0 0 1 0 0;
     0 0 0 0 1 0;
     0 0 0 0 0 1];

H = [1 0 0 0 0 0;
     0 1 0 0 0 0;
     0 0 1 0 0 0];

P = eye(6);
Q = 0.05*eye(6);
Rk = 25*eye(3);

for k = 2:M

    xPred = F*xKF(:,k-1);
    PPred = F*P*F' + Q;

    z = triPos(k,:)';

    K = PPred*H'/(H*PPred*H' + Rk);

    xKF(:,k) = xPred + K*(z - H*xPred);

    P = (eye(6)-K*H)*PPred;

end

kfPos = xKF(1:3,:)';

kfError = vecnorm(kfPos-targetValid,2,2);
kfRMSE = sqrt(mean(kfError.^2));

%% PERFORMANCE METRICS

meanTri = mean(errorVec);
meanKF = mean(kfError);

maxTri = max(errorVec);
maxKF = max(kfError);

fprintf('\n=====================================\n');
fprintf('PERFORMANCE COMPARISON\n');
fprintf('=====================================\n');
fprintf('Triangulation RMSE     : %.3f m\n',RMSE);
fprintf('Kalman RMSE            : %.3f m\n',kfRMSE);
fprintf('Triangulation Mean Err : %.3f m\n',meanTri);
fprintf('Kalman Mean Err        : %.3f m\n',meanKF);
fprintf('Triangulation Max Err  : %.3f m\n',maxTri);
fprintf('Kalman Max Err         : %.3f m\n',maxKF);
fprintf('Improvement            : %.2f %%\n',...
100*(RMSE-kfRMSE)/RMSE);
fprintf('=====================================\n');

%% DATASET

T = table();

T.TrueX = targetValid(:,1);
T.TrueY = targetValid(:,2);
T.TrueZ = targetValid(:,3);

T.TriX = triPos(:,1);
T.TriY = triPos(:,2);
T.TriZ = triPos(:,3);

T.KFX = kfPos(:,1);
T.KFY = kfPos(:,2);
T.KFZ = kfPos(:,3);

T.TriError = errorVec;
T.KFError = kfError;

T.ConditionNumber = conditionNum(~bad);

for s = 1:numSensors
    T.(sprintf('Az%d',s)) = azStore(~bad,s);
    T.(sprintf('El%d',s)) = elStore(~bad,s);
end

T.RunID = run*ones(height(T),1);

allData = [allData ; T];
fprintf('Rows so far = %d\n',height(allData));
%{
%% PLOT 1 - COMBINED

figure;
plot3(targetValid(:,1),targetValid(:,2),targetValid(:,3),'g','LineWidth',2);
hold on;
plot3(triPos(:,1),triPos(:,2),triPos(:,3),'r.');
plot3(kfPos(:,1),kfPos(:,2),kfPos(:,3),'b','LineWidth',2);
grid on;
legend('True Target','Triangulation','Kalman Filter');
title('Combined Comparison');
xlabel('X'); ylabel('Y'); zlabel('Z');

%% PLOT 2 - TRIANGULATION ONLY

figure;
plot3(targetValid(:,1),targetValid(:,2),targetValid(:,3),'g','LineWidth',2);
hold on;
plot3(triPos(:,1),triPos(:,2),triPos(:,3),'r.');
grid on;
legend('True Target','Triangulation');
title('Triangulation Only');
xlabel('X'); ylabel('Y'); zlabel('Z');

%% PLOT 3 - KALMAN ONLY

figure;
plot3(targetValid(:,1),targetValid(:,2),targetValid(:,3),'g','LineWidth',2);
hold on;
plot3(kfPos(:,1),kfPos(:,2),kfPos(:,3),'b','LineWidth',2);
grid on;
legend('True Target','Kalman Filter');
title('Kalman Filter Only');
xlabel('X'); ylabel('Y'); zlabel('Z');

%% PLOT 4 - ERROR COMPARISON

figure;
plot(errorVec,'r');
hold on;
plot(kfError,'b','LineWidth',1.5);
grid on;
legend('Triangulation','Kalman');
title('Error Comparison');
xlabel('Time Step');
ylabel('Position Error (m)');

%% PLOT 5 - RMSE BAR CHART

figure;
bar([RMSE kfRMSE]);
xticklabels({'Triangulation','Kalman'});
ylabel('RMSE (m)');
title('RMSE Comparison');
grid on;

figure;
scatter(conditionNum,errorVec);
xlabel('Condition Number');
ylabel('Triangulation Error');
title('Geometry vs Error');
grid on;
end
%}
close all
end

writetable(allData,'passive_ranging_dataset_big.csv');

save('passive_ranging_dataset_big.mat','allData');

fprintf('\n====================================\n');
fprintf('DATASET GENERATED SUCCESSFULLY\n');
fprintf('ROWS = %d\n',height(allData));
fprintf('====================================\n');