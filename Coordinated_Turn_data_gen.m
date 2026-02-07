rng(0, 'twister');

clear; clc; close all;
dt = 0.01;              
T_total = 120;           
time = 0:dt:T_total;   
N = length(time);
g = 9.81;              
Va_ref_start = 250; %air speed until step decrease
Va_ref_end = 170;   %air speed at the end of step_decrease
%True vals
Va_true_mag = zeros(1,N);     % TAS
Wind_N =-5.0;           %m/s
Wind_E = -10.0;         
Wind_D = 0.0;           
Vw_N_vec = [Wind_N; Wind_E; Wind_D]; 
%Sensor noises
SF_true =  1;     
noise_acc = 0.1;       
noise_gyro = 0.002;    
noise_gps = 0.2;       
noise_pitot = 1.0;      
%% Trajectory generation
%psi = yaw , theta = pitch, phi = roll;
phi = zeros(1, N);      
theta = zeros(1, N);    
psi = zeros(1, N);      
% coordinated turn --> 10-30 , step decrease 45 - 50;
idx_turn_start = round(10 / dt) + 1;
idx_turn_end   = round(70 / dt) + 1;
idx_step_dec_start = (100/dt) + 1;
idx_step_dec_end   = round(110/dt) + 1;
Va_true_mag(1:idx_step_dec_start) = Va_ref_start; %Step decrease başına kadar sabit 50;
Va_true_mag(idx_step_dec_start : idx_step_dec_end) = linspace(Va_ref_start,Va_ref_end,(idx_step_dec_end - idx_step_dec_start) + 1 ); % 50den 30 düşüş
Va_true_mag(idx_step_dec_end:end) = Va_ref_end;
% Yaw during trajectory
psi(1:idx_turn_start) = 0;
psi(idx_turn_start:idx_turn_end) = linspace(0, pi, idx_turn_end - idx_turn_start + 1);
psi(idx_turn_end:end) = pi;
% Yaw rate
turn_duration = time(idx_turn_end) - time(idx_turn_start);
psi_dot_val = pi / turn_duration; % rad/s
% Required Roll for coordinated turn
% tan(phi) = (V * psi_dot) / g
V_turn_avg = mean(Va_true_mag(idx_turn_start:idx_turn_end)); % Ortalama hız
phi_target = atan((V_turn_avg * psi_dot_val) / g);
%Transition time for smooth roll 
t_transition = 7.0; 
n_trans = round(t_transition / dt);
idx_flat_start = idx_turn_start + n_trans; 
idx_flat_end   = idx_turn_end - n_trans;   
phi(idx_turn_start : idx_flat_start) = linspace(0, phi_target, n_trans + 1);
phi(idx_flat_start : idx_flat_end) = phi_target;
phi(idx_flat_end : idx_turn_end) = linspace(phi_target, 0, n_trans + 1);
alpha_start =  0.3 * pi/180; % AoA before step decrease
n_load = 1./cos(phi);
alpha_rad = alpha_start * (Va_ref_start ./ Va_true_mag).^2 .* n_load;
theta = atan( tan(alpha_rad) .* cos(phi) );
%Derivatives
d_phi = gradient(phi, dt);
d_theta = gradient(theta, dt);
d_psi = gradient(psi, dt);
%% Velocity and acceleration calculations
p = zeros(1, N); q = zeros(1, N); r = zeros(1, N);
Va_b = zeros(3,N);
Va_N = zeros(3, N); % Airspeed in NED
Vg_N = zeros(3, N); % Groundspeedin NED
acc_body = zeros(3, N); % Accelerometer vals
% Velocity calc.
for k = 1:N
    
    ph = phi(k); th = theta(k); ps = psi(k);
    
    %(Euler Rates -> Body Rates)
    
    R_rates = [1,  0,       -sin(th);
               0,  cos(ph),  sin(ph)*cos(th);
               0, -sin(ph),  cos(ph)*cos(th)];
           
    euler_rates = [d_phi(k); d_theta(k); d_psi(k)];
    body_rates = R_rates * euler_rates;
    
    p(k) = body_rates(1);
    q(k) = body_rates(2);
    r(k) = body_rates(3);
    
    % Body -> NED)
    cph = cos(ph); sph = sin(ph);
    cth = cos(th); sth = sin(th);
    cps = cos(ps); sps = sin(ps);
    
    C_b_n = [cth*cps, sph*sth*cps - cph*sps, cph*sth*cps + sph*sps;
             cth*sps, sph*sth*sps + cph*cps, cph*sth*sps - sph*cps;
            -sth,     sph*cth,               cph*cth];
        
    % Body frame air speed
    Va_b(:,k) = [Va_true_mag(k) * cos(alpha_rad(k)); 0; Va_true_mag(k) * sin(alpha_rad(k))]; 
    
    % NED frame air speed
    Va_N(:, k) = C_b_n * Va_b(:,k);
    
    % Ground speed from velocity vector
    Vg_N(:, k) = Va_N(:, k) + Vw_N_vec;
end
% --- İvme Hesaplama (Reverse INS) ---
% İvme = d(Vg)/dt (Nav Frame)
% İvmeölçer yerçekimini de ölçer: f_b = C_n^b * (a_n - g_n)
% Yer Hızının Türevi (Nav Frame İvmesi)
acc_nav = gradient(Vg_N, dt); % 3xN matrix türevi
g_n = [0; 0; 9.81]; 
for k = 1:N
    %(NED -> Body) 
    ph = phi(k); th = theta(k); ps = psi(k);
    cph = cos(ph); sph = sin(ph);
    cth = cos(th); sth = sin(th);
    cps = cos(ps); sps = sin(ps);
    
    C_b_n = [cth*cps, sph*sth*cps - cph*sps, cph*sth*cps + sph*sps;
             cth*sps, sph*sth*sps + cph*cps, cph*sth*sps - sph*cps;
            -sth,     sph*cth,               cph*cth];
    
    C_n_b = C_b_n'; 
    
    
    % Real accel. extract gravity
    f_b = C_n_b * (acc_nav(:, k) - g_n);
    
    acc_body(:, k) = f_b;
end
%% Sensor datas
%INS
meas_acc = acc_body + noise_acc * randn(3, N);
meas_gyro = [p; q; r] + noise_gyro * randn(3, N);
%GPS
meas_Vg = Vg_N + noise_gps * randn(3, N);
%Pitot
meas_TAS = (Va_true_mag * SF_true) + ((Va_true_mag .* 0.01 + 0.05)./3) .* randn(1, N);
% meas_TAS = (Va_true_mag * SF_true) + noise_pitot.* randn(1, N);
%% Plot
figure(1);
subplot(2, 2, 1);
plot(time, psi * 180 / pi, 'LineWidth', 2);
grid on;
ylabel('Heading (deg)');
title('Heading Change');
subplot(2, 2, 2);
pos_E = cumsum(Va_N(2, :)) * dt;
pos_N = cumsum(Va_N(1, :)) * dt;
plot(pos_E, pos_N, 'LineWidth', 2);
grid on;
axis equal;
hold on;
quiver(pos_E(1:500:end), pos_N(1:500:end), ...
       ones(size(pos_E(1:500:end))) * Wind_E, ...
       ones(size(pos_N(1:500:end))) * Wind_N, ...
       0.5, 'r', 'LineWidth', 1);
xlabel('East (m)');
ylabel('North (m)');
title('Trajectory and Wind Vector');
legend('Trajectory', 'Wind', 'Location', 'best');
subplot(2, 2, 3);
plot(time, meas_Vg(1, :), 'r', time, meas_Vg(2, :), 'b');
grid on;
ylabel('GPS Speed (m/s)');
title('Ground Speed (NED)');
legend('V_{North}', 'V_{East}');
subplot(2, 2, 4);
plot(time, meas_TAS, 'k');
hold on;
plot(time, Va_true_mag, 'r', 'LineWidth', 1.5);
grid on;
ylabel('TAS (m/s)');
title('Measured Air Speed (with Scale Factor error)');
xlabel('Time (s)');
legend('Measured', 'True');
figure(2);
plot(time, alpha_rad * 180 / pi, 'LineWidth', 1.5);
grid on;
title('True Angle of Attack');
xlabel('Time (s)');
ylabel('Alpha (deg)');
%%  Extract data
sim_data.time = time;
sim_data.imu.acc = meas_acc;    
sim_data.imu.gyro = meas_gyro;  
sim_data.gps.vel = meas_Vg;     
sim_data.gps.sigma = noise_gps;
sim_data.pitot.tas = meas_TAS; 
sim_data.pitot.sigma = noise_pitot;
sim_data.truth.wind = Vw_N_vec;
sim_data.truth.sf = SF_true;
sim_data.truth.AoA = alpha_rad;
sim_data.att = [phi;theta;psi];
sim_data.truth.Va_b = Va_b;
save('wind_est_sim_data.mat', 'sim_data');
