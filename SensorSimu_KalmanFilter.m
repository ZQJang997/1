%% Aerospace Autonomy 497/597
% ZHANG,ZQ
% Homework 7, Question 2
clear all; close all hidden; clc;
%% SETTING UP OUR SIMULATION OF REALITY USING NON-LINEAR MODEL
dt = 0.01;
time = 0;
endtime = 150;
i = 1;

% restate our control gains, for convenience, as also used in derivative function
Kpsibyx = 0.011; % deg/(ft ⋅sec)
Kpsi = 0.0667;   % deg/(deg⋅sec)
Kv = 1;          % kt/kt/sec

% initialize our time history of state, with state vector containing:
% [Positionx(East)_feet; Positiony(North)_feet; V_knots; heading_degrees]
x(:, i) = [0 0 60 0]';

% initialize our controls
% x=750^', speed of 60 knots, and heading of 0
u = [750 60 0]';

%% SETTING UP OUR SIMULATION OF THE MEASUREMENT
%parameterize the statistics of our 4 measures
sigmaPositionSensor = 10; % feet; 2 of these, one to measure x and one to measure 6
sigmaVelocitySensor = 5;   % knots
sigmaHeadingSensor = 5;    % degrees

%set the time between measurements
mInterval = 1;
mCntr = 1;

%% INITIALIZING OUR KALMAN FILTER
% modelling the disturbances on V and Psi,
% note no direct disturbance on x & y
SigmaVDisturbance = 5;    % std dev of disturbance on v
SigmaPsiDisturbance = 12; % std dev of disturbance on psi

G = [0 0; 0 0; 1 0; 0 1];            

Q = [SigmaVDisturbance^2           0
              0         SigmaPsiDisturbance^2]; 

R = [100^2    0            0                         0
       0    100^2          0                         0
       0      0     sigmaVelocitySensor^2            0
       0      0            0              sigmaHeadingSensor^2];

e = [100;100;sigmaVelocitySensor;sigmaHeadingSensor];

% modeling the plant in LINEAR FORM to find the Kalman filter gain
% NOTE: conversion of knots to fps in elements multiplied by V
A = [0       0     0     u(2) * 6076/3600
     0       0  6076/3600      0
     0       0    -Kv          0
    -Kpsibyx 0     0         -Kpsi];

B = [0      0   0 
     0      0   0
     0      Kv  0
    Kpsibyx 0 Kpsi];

C = [1 0 0 0
     0 1 0 0
     0 0 1 0
     0 0 0 1]; % Since we are only measuring x and y, we need diagonal matrix including only for x and y

D = [0 0 0
     0 0 0
     0 0 0
     0 0 0];
Plant = ss(A, [B G], C, [D zeros(4, 2)], ...
    'inputname', {'Vc' 'Psic' 'Xc' 'WonV' 'WonPsi'}, ...
    'outputname', {'Measurementofx' 'Measurementofy' 'MeasurementofV' 'Measurementofpsi'});

[kalmf, K, Pinitial] = kalman(Plant, Q, R);
xest(:, i) = [0 0 60 0]';
Pactual(:, :, i) = [(x(1, i) - xest(1, i)) * (x(1, i) - xest(1, i)) 0 0 0; ...
                    0 (x(2, i) - xest(2, i)) * (x(2, i) - xest(2, i)) 0 0; ...
                    0 0 (x(3, i) - xest(3, i)) * (x(3, i) - xest(3, i)) 0; ...
                    0 0 0 (x(4, i) - xest(4, i)) * (x(4, i) - xest(4, i))];
P(:, :, i) = Pinitial;

while time < endtime
    %% SIMULATING REALITY USING NON-LINEAR MODEL
    % calculate the 4 estimates of derivative over the time interval from t to t+dt
    % in this case, x varies throughout the interval, but the elements of the control vector 'u' are constant
    xdot1 = derivativefunc(x(:, i), u, time);
    xdot2 = derivativefunc(x(:, i)+xdot1*dt/2, u, time+dt/2);  
    xdot3 = derivativefunc(x(:, i)+xdot2*dt/2, u, time+dt/2); 
    xdot4 = derivativefunc(x(:, i)+xdot3*dt, u, time+dt);  

    % calculate our 4th order Runge Kutta estimate of derivative
    totalxdot(:, i) = (xdot1 + 2 * xdot2 + 2 * xdot3 + xdot4) / 6;

    % adding disturbances to the dynamics of the aircraft
    w(1) = randn * SigmaVDisturbance;
    w(2) = randn * SigmaPsiDisturbance;
    totalxdot(:, i) = totalxdot(:, i) + G * w';

    % update state
    x(:, i + 1) = x(:, i) + totalxdot(:, i) * dt;

    %% HERE'S WHAT WOULD NEED TO BE COMPUTED IN REAL-TIME ONBOARD: THE KALMAN FILTER
    % Calculate the derivative of xest and of P, integrating using RK4
    % To calculate xest, we'll use our non-linear model
    xestdot1 = derivativefunc(xest(:, i), u, time);
    xestdot2 = derivativefunc(x(:, i)+xestdot1*dt/2, u, time+dt/2); 
    xestdot3 = derivativefunc(x(:, i)+xestdot2*dt/2, u, time+dt/2); 
    xestdot4 = derivativefunc(x(:, i)+xestdot3*dt, u, time+dt);
    totalxestdot = (xestdot1 + 2 * xestdot2 + 2 * xestdot3 + xestdot4) / 6;
    % DON'T INTEGRATE xestdot yet, we may have a measurement to include!

    % P, on the other hand, needs to know our linear model
    Pdot1 = A*P(:, :, i)*A' + G*Q*G'; 
    Pdot2 = A*(P(:, :, i)+Pdot1*(dt/2))*A' + G*Q*G'; 
    Pdot3 = A*(P(:, :, i)+Pdot2*(dt/2))*A' + G*Q*G'; 
    Pdot4 = A*(P(:, :, i)+Pdot3*dt)*A' + G*Q*G'; 
    totalPdot = (Pdot1 + 2 * Pdot2 + 2 * Pdot3 + Pdot4) / 6;
    P(:, :, i + 1) = P(:, :, i) + totalPdot * dt;

    %% IF A MEASUREMENT COMES IN, SIMULATING THE MEASUREMENT PROCESS AND USING MEASUREMENT TO CORRECT KALMAN FILTER
    if mod(time, mInterval) < dt
        time_arr_s(mCntr) = time;
        % Here's the sensor making the measurement
        y(1, mCntr) = x(1, i) + randn * sigmaPositionSensor;
        y(2, mCntr) = x(2, i) + randn * sigmaPositionSensor;
        y(3, mCntr) = x(3, i) + randn * sigmaVelocitySensor;
        y(4, mCntr) = x(4, i) + randn * sigmaHeadingSensor;

        % Here's where we correct the Kalman Filter, would be running onboard
        % Make our estimate of what the measures would be, from xest
        yest = C*xest(:, i); % !!HERE!!
        
        % calculate K
        K = P(:, :, i)* C' * (C * P(:, :, i) * C' + R)^-1; 
        
        % Add in the measurement's contribution to xestdot
        totalxestdot = totalxestdot + K*(y(:,mCntr)-yest(:)); 

        % correct our estimate of P
        P(:, :, i + 1) = (eye(4)-K*C)*P(:, :, i) ; 
        mCntr = mCntr + 1;
    end

    % Now we can finally integrate xest, now that we have both dynamics and
    % correction from measurements
    xest(:, i + 1) = xest(:, i) + totalxestdot * dt;

    % Pactual represents the actual error squared -- in real-life you wouldn't know this...
    % But here, in simulation, is a nice check to see how accurate the KF is.
    Pactual(:, :, i) = [(x(1, i) - xest(1, i)) * (x(1, i) - xest(1, i)) 0 0 0; ...
                        0 (x(2, i) - xest(2, i)) * (x(2, i) - xest(2, i)) 0 0; ...
                        0 0 (x(3, i) - xest(3, i)) * (x(3, i) - xest(3, i)) 0; ...
                        0 0 0 (x(4, i) - xest(4, i)) * (x(4, i) - xest(4, i))];

    i = i + 1;        % increment our iteration counter
    time = time + dt; % increment time
    time_arr(i) = time;

end

%% Find the statistics on the error in the estimate, compared to reality
fprintf('Error in X: \t\x03bc = % 8.3f\t\x03C3 = % 8.3f\n', mean(x(1, :) - xest(1, :)), std(x(1, :) - xest(1, :)))
fprintf('Error in Y: \t\x03bc = % 8.3f\t\x03C3 = % 8.3f\n', mean(x(2, :) - xest(2, :)), std(x(2, :) - xest(2, :)))
fprintf('Error in V: \t\x03bc = % 8.3f\t\x03C3 = % 8.3f\n', mean(x(3, :) - xest(3, :)), std(x(3, :) - xest(3, :)))
fprintf('Error in \x03C8: \t\x03bc = % 8.3f\t\x03C3 = % 8.3f\n', mean(x(4, :) - xest(4, :)), std(x(4, :) - xest(4, :)))

%% Finally, plot!
subplot(4, 1, [1, 2])
plot(x(1, :), x(2, :), y(1, :), y(2, :), xest(1, :), xest(2, :))
legend('Real Position', 'Measured position', 'Defualt KF estimate of position')
title('y by x position, feet')

subplot(4, 1, 3)
plot(time_arr, x(3, :), time_arr_s, y(3, :), time_arr, xest(3, :))
legend('Real Velocity', 'Measured Velocity', 'Defualt KF estimate of velocity')
title('Timeline of V, knots')

subplot(4, 1, 4)
plot(time_arr, x(4, :), time_arr_s, y(4, :), time_arr, xest(4, :), '-.')
legend('Real heading', 'Measured Heading', 'Defualt KF estimate of heading')
title('Timeline of \Psi, degrees')