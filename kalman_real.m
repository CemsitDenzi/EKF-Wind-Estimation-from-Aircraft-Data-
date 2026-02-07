%% 1. SETUP AND INITIALIZATION
clc; clear; close all;
filename = "661c4485-6ec7-4471-9f97-3ef7bb6f51f1.ulg";
if ~isfile(filename), error('File not found Check the name'); end
disp('1 Reading ULog Data');
ulog = ulogreader(filename);
gps_raw  = readTopicMsgs(ulog, 'TopicNames', {'vehicle_local_position'}, 'InstanceID', {0});
air_raw  = readTopicMsgs(ulog, 'TopicNames', {'airspeed'},               'InstanceID', {0});
att_raw  = readTopicMsgs(ulog, 'TopicNames', {'vehicle_attitude'},       'InstanceID', {0});
imu_raw  = readTopicMsgs(ulog, 'TopicNames', {'sensor_combined'},        'InstanceID', {0});
wind_raw = readTopicMsgs(ulog, 'TopicNames', {'wind'},                   'InstanceID', {0});

%% 2. SAFE DATA EXTRACTION
att_tbl = get_table_safe(att_raw, 'Attitude');
imu_tbl = get_table_safe(imu_raw, 'IMU');
gps_tbl = get_table_safe(gps_raw, 'GPS');
air_tbl = get_table_safe(air_raw, 'Airspeed');
if ~isempty(wind_raw)
    wind_tbl = get_table_safe(wind_raw, 'Wind');
else
    wind_tbl = [];
end
t_att = get_time_sec(att_tbl);
t_imu = get_time_sec(imu_tbl);
t_gps = get_time_sec(gps_tbl);
t_air = get_time_sec(air_tbl);
if ~isempty(wind_tbl), t_wind = get_time_sec(wind_tbl); else, t_wind = []; end
disp('Data and Times Extracted');

%% 3. RAW DATA PREPARATION
if ismember('q', att_tbl.Properties.VariableNames), q_val = att_tbl.q; 
else, q_val = [att_tbl.q_0, att_tbl.q_1, att_tbl.q_2, att_tbl.q_3]; end
eul_val = quat2eul(q_val, 'ZYX');
if ismember('gyro_rad', imu_tbl.Properties.VariableNames)
    gyro_val = imu_tbl.gyro_rad; acc_val = imu_tbl.accelerometer_m_s2;
else
    gyro_val = [imu_tbl.gyro_rad_0, imu_tbl.gyro_rad_1, imu_tbl.gyro_rad_2];
    acc_val  = [imu_tbl.accelerometer_m_s2_0, imu_tbl.accelerometer_m_s2_1, imu_tbl.accelerometer_m_s2_2];
end
if ismember('true_airspeed_m_s', air_tbl.Properties.VariableNames)
    tas_val = air_tbl.true_airspeed_m_s;
else
    tas_val = air_tbl.indicated_airspeed_m_s;
end

%% 4. SYNCHRONIZATION AND INTERPOLATION
disp('2 Synchronizing Data Loiter Interval');
t_start = 4260; 
t_end = 4980;
idx = t_att >= t_start & t_att <= t_end;
if sum(idx) == 0, error('No data in selected time interval'); end
t_ref = t_att(idx);
dt_target = 0.02;
time_common = (t_ref(1):dt_target:t_ref(end))';
interp = @(t,y) interp1(t, double(y), time_common, 'linear', 'extrap');
flight.time = time_common;
flight.phi   = interp(t_att, unwrap(eul_val(:,3)));
flight.theta = interp(t_att, unwrap(eul_val(:,2)));
flight.psi   = interp(t_att, unwrap(eul_val(:,1)));
flight.p = interp(t_imu, gyro_val(:,1));
flight.q = interp(t_imu, gyro_val(:,2));
flight.r = interp(t_imu, gyro_val(:,3));
flight.ax = interp(t_imu, acc_val(:,1));
flight.ay = interp(t_imu, acc_val(:,2));
flight.az = interp(t_imu, acc_val(:,3));
flight.Vn = interp(t_gps, gps_tbl.vx);
flight.Ve = interp(t_gps, gps_tbl.vy);
flight.Vd = interp(t_gps, gps_tbl.vz);
flight.TAS = interp(t_air, tas_val);
if ~isempty(t_wind)
    flight.Wn_ref = interp(t_wind, wind_tbl.windspeed_north);
    flight.We_ref = interp(t_wind, wind_tbl.windspeed_east);
else
    flight.Wn_ref = zeros(size(time_common));
    flight.We_ref = zeros(size(time_common));
end

%% 5. DATA MAPPING
disp('3 Converting Data to Kalman Format');
sim_data.time = flight.time'; 
sim_data.dt   = dt_target;
sim_data.imu.acc  = [flight.ax'; flight.ay'; flight.az'];
sim_data.imu.gyro = [flight.p'; flight.q'; flight.r'];
sim_data.gps.vel  = [flight.Vn'; flight.Ve'; flight.Vd'];
sim_data.pitot.tas = flight.TAS';
sim_data.att = [flight.phi'; flight.theta'; flight.psi'];
sim_data.truth.wind = [flight.Wn_ref'; flight.We_ref'];

%% 6. KALMAN FILTER
disp('4 Running Kalman Filter');
time_vec = sim_data.time;
dt = sim_data.dt;
N = length(time_vec);
acc_meas  = sim_data.imu.acc;   
gyro_meas = sim_data.imu.gyro;  
gps_meas  = sim_data.gps.vel;  
pitot_meas= sim_data.pitot.tas; 
W_N_ref = sim_data.truth.wind(1,:);
W_E_ref = sim_data.truth.wind(2,:);
Attitude_meas = sim_data.att; 
x = zeros(7, 1);
V_start = pitot_meas(1);     
alpha_guess_deg = 2.0;       
alpha_guess_rad = alpha_guess_deg * pi/180;
x(1) = V_start * cos(alpha_guess_rad); 
x(2) = 0;                              
x(3) = V_start * sin(alpha_guess_rad); 
x(4) = 1.0;                            
P = eye(7);
P(1:3, 1:3) = 10 * eye(3);
P(4,4)      = 0.0;       
P(5:7, 5:7) = 5 * eye(3);
Q = diag([...
    0.5, 0.5, 0.5, ...      
    0, ...                  
    1e-8, 1e-8, 1e-8 ...    
    ]);
R = diag([...
    2.0^2, ...               
    0.5^2, 0.5^2, 0.5^2, ... 
    10.0^2 ...               
    ]);
hist_x = zeros(7, N);
Va_estimated = zeros(1,N);
for k = 1:N
    ax = acc_meas(1, k); ay = acc_meas(2, k); az = acc_meas(3, k);
    p  = gyro_meas(1, k); q  = gyro_meas(2, k); r  = gyro_meas(3, k);
    V_tas_meas = pitot_meas(k);
    V_gps_meas = gps_meas(:, k);
    roll  = Attitude_meas(1, k);
    pitch = Attitude_meas(2, k);
    yaw   = Attitude_meas(3, k);
    F = zeros(7,7);
    F(1:3, 1:3) = [0, r, -q; -r, 0, p; q, -p, 0];
    u_input = zeros(7,1);
    u_input(1) = ax - 9.81 * sin(pitch);
    u_input(2) = ay + 9.81 * sin(roll) * cos(pitch);
    u_input(3) = az + 9.81 * cos(roll) * cos(pitch);
    x = x + (F * x + u_input) * dt;
    F_discrete = eye(7) + F * dt;
    P = F_discrete * P * F_discrete' + Q;
    C_BN = angle2dcm(yaw, pitch, roll)'; 
    ua = x(1); va = x(2); wa = x(3); sf = x(4); V_wind = x(5:7);
    V_total = sqrt(ua^2 + va^2 + wa^2);
    z_pred = zeros(5,1);
    z_pred(1) = sf * V_total;                   
    z_pred(2:4) = C_BN * [ua; va; wa] + V_wind; 
    z_pred(5) = sf * va;                        
    z_meas = [V_tas_meas; V_gps_meas; 0]; 
    H = zeros(5, 7);
    if V_total > 0.1
        H(1, 1:3) = (sf / V_total) * [ua, va, wa];
        H(1, 4)   = V_total; 
    end
    H(2:4, 1:3) = C_BN;
    H(2:4, 5:7) = eye(3);
    H(5, 2) = sf; H(5, 4) = va;
    y = z_meas - z_pred;
    S = H * P * H' + R;
    K = P * H' / S;
    x = x + K * y;
    P = (eye(7) - K * H) * P;
    Va_estimated(k) = x(4) * sqrt(x(1)^2 + x(2)^2 + x(3)^2);
    hist_x(:, k) = x;
end
disp('Kalman Filter Completed');

%% 7. PLOTTING RESULTS
figure('Name', 'Wind Estimation Results', 'Color', 'w');
subplot(2,2,1); plot(time_vec, hist_x(5,:), 'b', 'LineWidth', 2); hold on;
plot(time_vec, W_N_ref, 'r--', 'LineWidth', 1.5); title('North Wind (m/s)'); legend('Kalman', 'PX4'); grid on;
subplot(2,2,2); plot(time_vec, hist_x(6,:), 'b', 'LineWidth', 2); hold on;
plot(time_vec, W_E_ref, 'r--', 'LineWidth', 1.5); title('East Wind (m/s)'); grid on;
subplot(2,2,3); plot(time_vec, hist_x(4,:), 'k'); title('Scale Factor'); grid on;
subplot(2,2,4); plot(gps_meas(2,:), gps_meas(1,:), 'Color', [0.8 0.8 0.8]); hold on;
quiver(0, 0, hist_x(6,end)*20, hist_x(5,end)*20, 'r', 'LineWidth', 2);
title('Path & Wind Vector'); axis equal; grid on;

% additional functions 
function tbl = get_table_safe(raw_topic, name)
    if isempty(raw_topic), error('ERROR %s Topic is empty', name); end
    content = raw_topic.TopicMessages;
    if iscell(content), tbl = content{1}; else, tbl = content; end
end
function t_sec = get_time_sec(tbl)
    if istimetable(tbl)
        raw_time = tbl.Properties.RowTimes;
        if isduration(raw_time), t_sec = seconds(raw_time);
        else, t_sec = double(raw_time); end
    elseif ismember('timestamp', tbl.Properties.VariableNames)
        t_sec = double(tbl.timestamp) / 1e6;
    elseif ismember('Timestamp', tbl.Properties.VariableNames)
        t_sec = double(tbl.Timestamp) / 1e6;
    else
        error('Time column not found');
    end
end