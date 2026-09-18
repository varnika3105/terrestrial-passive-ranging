clc;
clear;
close all;

%% =====================================================
% REALISTIC TERRESTRIAL PASSIVE RANGING DATASET
% MULTI SENSOR 3D TRIANGULATION
% NO FILTERS
%% =====================================================

rng(1);

%% =====================================================
% SIMULATION PARAMETERS
%% =====================================================

N = 1000;
dt = 1;

numSensors = 4;

%% =====================================================
% TARGET TRAJECTORY
%% =====================================================

targetPos = zeros(N,3);
targetVel = zeros(N,3);

targetPos(1,:) = [500 500 20];
targetVel(1,:) = [8 3 0];

accelStd = 0.2;

for k = 2:N

    acc = accelStd*randn(1,3);

    targetVel(k,:) = ...
        targetVel(k-1,:) + acc*dt;

    targetPos(k,:) = ...
        targetPos(k-1,:) + targetVel(k,:)*dt;

end

%% =====================================================
% SENSOR LOCATIONS
%% =====================================================

sensorTrue = zeros(N,3,numSensors);

for k = 1:N

    sensorTrue(k,:,1) = [0 0 5];

    sensorTrue(k,:,2) = [2000 0 5];

    sensorTrue(k,:,3) = [2000 2000 5];

    sensorTrue(k,:,4) = [0 2000 5];

end

%% =====================================================
% SENSOR POSITION ERROR
%% =====================================================

gpsSigma = 2;

sensorPos = sensorTrue + ...
    gpsSigma*randn(size(sensorTrue));

%% =====================================================
% ANGLE ERRORS
%% =====================================================

azNoise = deg2rad(0.5);
elNoise = deg2rad(0.3);

azBias = deg2rad(0.2)*sin((1:N)'/200);
elBias = deg2rad(0.15)*cos((1:N)'/300);

%% =====================================================
% STORAGE
%% =====================================================

triPos = zeros(N,3);

azStore = zeros(N,numSensors);
elStore = zeros(N,numSensors);

%% =====================================================
% TRIANGULATION
%% =====================================================

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

        %% TRUE ANGLES

        az = atan2(dy,dx);

        el = atan2( ...
            dz,...
            sqrt(dx^2+dy^2));

        %% ADD REALISTIC ERRORS

        azMeas = ...
            az + ...
            azNoise*randn + ...
            azBias(k);

        elMeas = ...
            el + ...
            elNoise*randn + ...
            elBias(k);

        azStore(k,s) = azMeas;
        elStore(k,s) = elMeas;

        %% LOS VECTOR

        ux = cos(elMeas)*cos(azMeas);
        uy = cos(elMeas)*sin(azMeas);
        uz = sin(elMeas);

        u = [ux uy uz];

        %% TRIANGULATION CONSTRAINT

        Ai = eye(3) - (u'*u);

        A = [A; Ai];

        b = [b; Ai*[sx;sy;sz]];

    end

    triPos(k,:) = (A\b)';

end

%% =====================================================
% ERROR ANALYSIS
%% =====================================================

errorVec = vecnorm( ...
    triPos-targetPos,...
    2,...
    2);

RMSE = sqrt(mean(errorVec.^2));

fprintf('\n');
fprintf('=====================================\n');
fprintf('TRIANGULATION RESULTS\n');
fprintf('=====================================\n');
fprintf('Samples : %d\n',N);
fprintf('RMSE    : %.3f m\n',RMSE);
fprintf('=====================================\n');

%% =====================================================
% SAVE COMPLETE DATASET
%% =====================================================

data = struct();

data.targetPos = targetPos;
data.targetVel = targetVel;

data.sensorTrue = sensorTrue;
data.sensorPos = sensorPos;

data.azimuth = azStore;
data.elevation = elStore;

data.triangulatedPos = triPos;
data.error = errorVec;
data.rmse = RMSE;

save('terrestrial_passive_ranging_dataset.mat','data');

%% =====================================================
% CSV EXPORT
%% =====================================================

T = table();

T.TrueX = targetPos(:,1);
T.TrueY = targetPos(:,2);
T.TrueZ = targetPos(:,3);

T.TriX = triPos(:,1);
T.TriY = triPos(:,2);
T.TriZ = triPos(:,3);

T.Error = errorVec;

for s = 1:numSensors

    T.(['Az' num2str(s)]) = azStore(:,s);
    T.(['El' num2str(s)]) = elStore(:,s);

end

writetable(T,...
'terrestrial_passive_ranging_dataset.csv');

%% =====================================================
% PLOTS
%% =====================================================

figure;

plot3( ...
targetPos(:,1), ...
targetPos(:,2), ...
targetPos(:,3), ...
'g','LineWidth',2);

hold on;

plot3( ...
triPos(:,1), ...
triPos(:,2), ...
triPos(:,3), ...
'r.');

grid on;

legend('True Target',...
       'Triangulation');

xlabel('X (m)');
ylabel('Y (m)');
zlabel('Z (m)');

title('Terrestrial Passive Ranging');

view(3);

%% =====================================================
% ERROR PLOT
%% =====================================================

figure;

plot(errorVec);

grid on;

xlabel('Time Step');
ylabel('Position Error (m)');

title('Triangulation Error');
