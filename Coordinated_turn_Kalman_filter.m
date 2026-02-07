clc; clear; close all;
%IVV TH
%% Load datas

if isfile('wind_est_sim_data.mat')
    load('wind_est_sim_data.mat');
else
    error('Lütfen önce veri seti oluşturucu scripti çalıştırın!');
end

time_vec = sim_data.time;
dt = time_vec(2) - time_vec(1);
N = length(time_vec);

acc_meas  = sim_data.imu.acc;   
gyro_meas = sim_data.imu.gyro;  
gps_meas  = sim_data.gps.vel;  
pitot_meas= sim_data.pitot.tas; 

% True vals
W_N_true = sim_data.truth.wind(1);
W_E_true = sim_data.truth.wind(2);
SF_true  = sim_data.truth.sf;
Attitude_true = sim_data.att; 

%% Initals
% States -> [ua; va; wa; sf; VwN; VwE; VwD]
x = zeros(7, 1);

V_start = pitot_meas(1);     
alpha_guess_deg = 0.25;       
alpha_guess_rad = alpha_guess_deg * pi/180;

x = zeros(7, 1);
x(1) = V_start * cos(alpha_guess_rad); % u 
x(2) = 0;                              % v
x(3) = V_start * sin(alpha_guess_rad); % w 
x(4) = 1.0;                            % sf


P = eye(7);
P(1:3, 1:3) = 100 * eye(3);
P(4,4)      = 0.01;       
P(5:7, 5:7) = 10 * eye(3);

% Process Noise (Q)
Q = diag([...
    0.01, 0.01, 0.01, ... % ua, va, wa
    1e-5, ...             % Scale Factor 
    1e-5, 1e-5, 1e-5 ...  % wind
    ]);


%Measurement Noise
R = diag([...
    1.0^2, ...            % TAS 
    0.2^2, 0.2^2, 0.2^2, ... % GPS 
    0.5 ...              % v_a
    ]);
hist_x = zeros(7, N);
hist_sigma = zeros(7, N); 
Va_estimated = zeros(1,N);

%% 3. KALMAN FILTER

for k = 1:N

    ax = acc_meas(1, k); ay = acc_meas(2, k); az = acc_meas(3, k);
    p  = gyro_meas(1, k); q  = gyro_meas(2, k); r  = gyro_meas(3, k);
    
    V_tas_meas = pitot_meas(k);
    V_gps_meas = gps_meas(:, k);
    
    % degisecek
    roll  = Attitude_true(1, k);
    pitch = Attitude_true(2, k);
    yaw   = Attitude_true(3, k);
    

    ua = x(1); va = x(2); wa = x(3);
    
    F = zeros(7,7);
    F(1:3, 1:3) = [0, r, -q; -r, 0, p; q, -p, 0];
    
    u_input = zeros(7,1);
    u_input(1) = ax - 9.81 * sin(pitch);
    u_input(2) = ay + 9.81 * sin(roll) * cos(pitch);
    u_input(3) = az + 9.81 * cos(roll) * cos(pitch);
    
    x_dot = F * x + u_input;
    x = x + x_dot * dt;
    
    F_discrete = eye(7) + F * dt;
    P = F_discrete * P * F_discrete' + Q;
    

    C_BN = angle2dcm(yaw, pitch, roll)'; 
    
    ua = x(1); va = x(2); wa = x(3); sf = x(4);
    V_wind = x(5:7);
    
    V_total = sqrt(ua^2 + va^2 + wa^2);
    R(1,1) = (V_total* 0.01 + 0.05) / 3 ;
    
    %predictions
    z_pred = zeros(5,1);
    z_pred(1) = sf * V_total;                   % TAS 
    z_pred(2:4) = C_BN * [ua; va; wa] + V_wind; % GPS 
    z_pred(5) = sf * va;                        % sideslip
    
    % Measurement
    z_meas = [V_tas_meas; V_gps_meas; 0]; 
    
    
    H = zeros(5, 7);
    
    % H - TAS (V_meas = sf * V_tot)
    if V_total > 0.1
        H(1, 1:3) = (sf / V_total) * [ua, va, wa];
        H(1, 4)   = V_total; 
    end
    
    % H - GPS
    H(2:4, 1:3) = C_BN;
    H(2:4, 5:7) = eye(3);
    
    % H - sideslip
    H(5, 2) = sf;
    H(5, 4) = va;
    
    % Kalman Gain
    y = z_meas - z_pred;
    S = H * P * H' + R;
    K = P * H' / S;
    
    x = x + K * y;
    P = (eye(7) - K * H) * P;

    Va_estimated(k) = x(4) * sqrt(x(1)^2 + x(2)^2 + x(3)^2);
   
    hist_x(:, k) = x;
    hist_sigma(:, k) = sqrt(diag(P));
end

%% Plot
figure(1);

subplot(2,2,1);
plot(time_vec, hist_x(5,:), 'LineWidth', 2); hold on;
yline(W_N_true, 'r--', 'LineWidth', 2);
title('North Wind (V_{wN})'); ylabel('m/s'); grid on;
legend('Estimate', 'True');

subplot(2,2,2);
plot(time_vec, hist_x(6,:), 'LineWidth', 2); hold on;
yline(W_E_true, 'r--', 'LineWidth', 2);
title('East Wind (V_{wE})'); ylabel('m/s'); grid on;

subplot(2,2,3);
plot(time_vec, hist_x(4,:), 'b', 'LineWidth', 2); hold on;
yline(SF_true, 'r--', 'LineWidth', 2);
yline(1.0, 'k:', 'Nominal (1.0)');
title('Pitot Scale Factor (SF)'); ylabel('Ratio'); grid on;
ylim([0.9 1.1]);

% % --- 4. YÖRÜNGE ---
% subplot(2,2,4);
% plot(sim_data.gps.vel(2,:), sim_data.gps.vel(1,:), 'color', [0.7 0.7 0.7]); hold on;
% % Rüzgar vektörünü çiz
% quiver(0, 0, hist_x(6,end)*10, hist_x(5,end)*10, 'r', 'LineWidth', 2, 'MaxHeadSize', 0.5);
% text(0, 0, '  Rüzgar', 'Color', 'r');
% title('Hız Uzayı ve Rüzgar Vektörü'); 
% xlabel('Doğu Hızı'); ylabel('Kuzey Hızı'); grid on; axis equal;

subplot(2,2,4)

plot(time_vec, atan(hist_x(3,:)./hist_x(1,:)),'b',LineWidth=1.2);hold on;
plot(time_vec, sim_data.truth.AoA,"r--",LineWidth=1.2);
title("AoA")
grid on;

figure(2)
subplot(1,3,1)
plot(time_vec,sim_data.truth.Va_b(1,:) - hist_x(1,:));hold on;
plot(time_vec,3 * hist_sigma(1,:),'r--',LineWidth=2)
plot(time_vec,-3*hist_sigma(1,:),'r--',LineWidth=2)
grid on;
title("ua error")

subplot(1,3,2)
plot(time_vec,sim_data.truth.Va_b(2,:) - hist_x(2,:));hold on;
plot(time_vec,3 * hist_sigma(2,:),'r--',LineWidth=2)
plot(time_vec,-3*hist_sigma(2,:),'r--',LineWidth=2)
grid on;
title("va error")

subplot(1,3,3)
plot(time_vec,sim_data.truth.Va_b(3,:) - hist_x(3,:));hold on;
plot(time_vec,3 * hist_sigma(3,:),'r--',LineWidth=2)
plot(time_vec,-3*hist_sigma(3,:),'r--',LineWidth=2)
grid on;
title("wa error")

grid on; 