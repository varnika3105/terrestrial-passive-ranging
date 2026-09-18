clc;
clear;
close all;

rng(1);

%% =====================================================
% RESEARCH GRADE MOVING SENSOR PASSIVE RANGING
% ANGLE ONLY TRIANGULATION
%% =====================================================

N = 20000;
dt = 1;

numSensors = 4;

%% =====================================================
% TARGET MOTION
%% =====================================================

targetPos = zeros(N,3);
targetVel = zeros(N,3);

targetPos(1,:) = [1000 1000 20];
targetVel(1,:) = [2 1 0];

accStd = 0.05;

for k = 2:N

    if mod(k,200)==0
        targetVel(k-1,:) = ...
            targetVel(k-1,:) + [0.5 -0.3 0];
    end

    a = accStd*randn(1,3);

    targetVel(k,:) = ...
        targetVel(k-1,:) + a*dt;

    targetPos(k,:) = ...
        targetPos(k-1,:) + targetVel(k,:)*dt;

end

%% =====================================================
% MOVING SENSOR TRAJECTORIES
%% =====================================================

sensorTrue = zeros(N,3,numSensors);

R = 1200;

for k = 1:N

    ang = 0.004*k;

    % Sensor 1

    sensorTrue(k,:,1) = ...
        [2000 + R*cos(ang), ...
         2000 + R*sin(ang), ...
         20];

    % Sensor 2

    sensorTrue(k,:,2) = ...
        [2000 - R*cos(ang), ...
         2000 - R*sin(ang), ...
         20];

    % Sensor 3

    sensorTrue(k,:,3) = ...
        [1000 + 300*sin(0.01*k), ...
         3500, ...
         30];

    % Sensor 4

    sensorTrue(k,:,4) = ...
        [3500,...
         1000 + 300*cos(0.01*k),...
         25];

end

%% =====================================================
% SENSOR POSITION ERROR
%% =====================================================

gpsSigma = 1;

sensorPos = sensorTrue + ...
            gpsSigma*randn(size(sensorTrue));

%% =====================================================
% MEASUREMENT ERROR MODEL
%% =====================================================

azNoise = deg2rad(0.3);
elNoise = deg2rad(0.2);

azBias = deg2rad(0.15)* ...
         sin((1:N)'/250);

elBias = deg2rad(0.10)* ...
         cos((1:N)'/300);

%% =====================================================
% STORAGE
%% =====================================================

triPos = zeros(N,3);

azStore = zeros(N,numSensors);
elStore = zeros(N,numSensors);

conditionNum = zeros(N,1);

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

        %% RANDOM NOISE

        azMeas = ...
            az + ...
            azNoise*randn + ...
            azBias(k);

        elMeas = ...
            el + ...
            elNoise*randn + ...
            elBias(k);

        %% OCCASIONAL OUTLIERS

        if rand < 0.01

            azMeas = ...
                azMeas + deg2rad(2);

            elMeas = ...
                elMeas + deg2rad(1);

        end

        azStore(k,s) = azMeas;
        elStore(k,s) = elMeas;

        %% LOS VECTOR

        ux = cos(elMeas)*cos(azMeas);
        uy = cos(elMeas)*sin(azMeas);
        uz = sin(elMeas);

        u = [ux uy uz];

        %% WEIGHTED CONSTRAINT

        Ai = eye(3) - (u'*u);

        w = 1/(azNoise^2 + elNoise^2);

        A = [A ; sqrt(w)*Ai];

        b = [b ; sqrt(w)*Ai*[sx;sy;sz]];

    end

    %% GEOMETRY CHECK

    conditionNum(k) = cond(A);

    if conditionNum(k) > 1e6

        triPos(k,:) = [NaN NaN NaN];

    else

        triPos(k,:) = (A\b)';

    end

end

%% =====================================================
% REMOVE BAD GEOMETRY POINTS
%% =====================================================

bad = isnan(triPos(:,1));

triPos(bad,:) = [];

targetValid = targetPos(~bad,:);

%% =====================================================
% ERROR ANALYSIS
%% =====================================================

errorVec = vecnorm( ...
    triPos-targetValid,...
    2,...
    2);

RMSE = sqrt(mean(errorVec.^2));

fprintf('\n');
fprintf('===========================\n');
fprintf('RMSE = %.3f m\n',RMSE);
fprintf('===========================\n');

%% =====================================================
% SAVE EVERYTHING
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

data.conditionNumber = conditionNum;

save('moving_sensor_dataset.mat','data');

%% =====================================================
% ENHANCED CSV EXPORT
%% =====================================================

T = table();

%% -----------------------------------------------------
% TRUE TARGET POSITION
%% -----------------------------------------------------

T.TrueX = targetValid(:,1);
T.TrueY = targetValid(:,2);
T.TrueZ = targetValid(:,3);

%% -----------------------------------------------------
% TARGET VELOCITY
%% -----------------------------------------------------

T.Vx = targetVel(~bad,1);
T.Vy = targetVel(~bad,2);
T.Vz = targetVel(~bad,3);

%% -----------------------------------------------------
% TRIANGULATED POSITION
%% -----------------------------------------------------

T.TriX = triPos(:,1);
T.TriY = triPos(:,2);
T.TriZ = triPos(:,3);

%% -----------------------------------------------------
% POSITION ERROR
%% -----------------------------------------------------

T.Error = errorVec;

%% -----------------------------------------------------
% GEOMETRY QUALITY
%% -----------------------------------------------------

T.ConditionNumber = conditionNum(~bad);

%% -----------------------------------------------------
% SENSOR POSITIONS
%% -----------------------------------------------------

for s = 1:numSensors

    T.(['S' num2str(s) '_X']) = ...
        squeeze(sensorPos(~bad,1,s));

    T.(['S' num2str(s) '_Y']) = ...
        squeeze(sensorPos(~bad,2,s));

    T.(['S' num2str(s) '_Z']) = ...
        squeeze(sensorPos(~bad,3,s));

end

%% -----------------------------------------------------
% AZIMUTH / ELEVATION
%% -----------------------------------------------------

for s = 1:numSensors

    T.(['Az' num2str(s)]) = ...
        azStore(~bad,s);

    T.(['El' num2str(s)]) = ...
        elStore(~bad,s);

end

%% -----------------------------------------------------
% LOS UNIT VECTORS
%% -----------------------------------------------------

for s = 1:numSensors

    az = azStore(~bad,s);
    el = elStore(~bad,s);

    ux = cos(el).*cos(az);
    uy = cos(el).*sin(az);
    uz = sin(el);

    T.(['Ux' num2str(s)]) = ux;
    T.(['Uy' num2str(s)]) = uy;
    T.(['Uz' num2str(s)]) = uz;

end

%% -----------------------------------------------------
% SAVE CSV
%% -----------------------------------------------------

writetable(T,...
    'moving_sensor_dataset_full.csv');

fprintf('\n===================================\n');
fprintf('Enhanced Dataset Saved\n');
fprintf('Rows : %d\n',height(T));
fprintf('Cols : %d\n',width(T));
fprintf('===================================\n');

%% =====================================================
% VISUALIZATION
%% =====================================================

figure;

plot3( ...
    targetValid(:,1), ...
    targetValid(:,2), ...
    targetValid(:,3), ...
    'g','LineWidth',2);

hold on;

plot3( ...
    triPos(:,1), ...
    triPos(:,2), ...
    triPos(:,3), ...
    'r.');

grid on;

xlabel('X (m)');
ylabel('Y (m)');
zlabel('Z (m)');

legend('True Target',...
       'Triangulated');

title('Moving Sensor Passive Ranging');

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

%% =====================================================
% GEOMETRY HEALTH
%% =====================================================

figure;

semilogy(conditionNum);

grid on;

xlabel('Time Step');
ylabel('Condition Number');

title('Geometry Quality');