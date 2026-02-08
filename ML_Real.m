%% 1. READING ULOG FILE AND PREPARATION
clc; clear; close all;
filename = "661c4485-6ec7-4471-9f97-3ef7bb6f51f1.ulg";
if ~isfile(filename), error('File not found!'); end
ulog = ulogreader(filename);
gps_raw  = readTopicMsgs(ulog, 'TopicNames', {'vehicle_local_position'}, 'InstanceID', {0});
air_raw  = readTopicMsgs(ulog, 'TopicNames', {'airspeed'},               'InstanceID', {0});
gps_tbl = get_table_safe(gps_raw, 'GPS');
air_tbl = get_table_safe(air_raw, 'Airspeed');
t_gps = get_time_sec(gps_tbl);
t_air = get_time_sec(air_tbl);
if ismember('true_airspeed_m_s', air_tbl.Properties.VariableNames)
    tas_raw = air_tbl.true_airspeed_m_s;
else
    tas_raw = air_tbl.indicated_airspeed_m_s;
end
%% 2. Take Loiter Part
t_start = 4260; 
t_end = 4980;
dt = 0.1; 
sim_time = (t_start:dt:t_end);
vx_double = double(gps_tbl.vx);
vy_double = double(gps_tbl.vy);
tas_double = double(tas_raw);
meas_vg_n = interp1(t_gps, vx_double, sim_time, 'linear', 'extrap'); 
meas_vg_e = interp1(t_gps, vy_double, sim_time, 'linear', 'extrap'); 
meas_tas  = interp1(t_air, tas_double, sim_time, 'linear', 'extrap'); 
total_meas_vg = [meas_vg_n; meas_vg_e]; 
total_measured_TAS = meas_tas;          
sig_g = 0.2;  
sig_a = 1.0;  
sig_w = 0.1;  
%% 3. OPTIMIZATION LOOP (LSQNONLIN)
N_size = 300; 
wx_hist = zeros(size(sim_time));
wy_hist = zeros(size(sim_time));
w_x0 = 0;
w_y0 = 0; 
step_stride = 5; 
loop_indices = (N_size + 1) : step_stride : size(sim_time,2) - N_size;
total_loops = length(loop_indices);
counter = 0;
for i = loop_indices
    counter = counter + 1;
    idx = i - N_size : i + N_size;
    
    meas_vg_x = total_meas_vg(1,idx).'; 
    meas_vg_y = total_meas_vg(2,idx).'; 
    meas_TAS  = total_measured_TAS(idx).';
    
    func = @(x) ML_func(x,meas_vg_x,meas_vg_y,meas_TAS,sig_a,sig_g,sig_w,w_x0,w_y0);
    
    options = optimoptions('lsqnonlin', 'Display', 'none', ...
        'FunctionTolerance', 1e-4, 'StepTolerance', 1e-4);
    
    x0 = double([meas_vg_x; meas_vg_y; w_x0; w_y0]);
    
    x_est = lsqnonlin(func, x0, [], [], options);
    
    w_x0 = x_est(end-1);
    w_y0 = x_est(end);
    
    fill_range = i : min(i + step_stride - 1, length(sim_time));
    wx_hist(fill_range) = w_x0;
    wy_hist(fill_range) = w_y0;
end

%% 4. PLOTTING RESULTS
if ~isempty(ulog.readTopicMsgs('TopicNames', {'wind'}, 'InstanceID', {0}))
    wind_raw = readTopicMsgs(ulog, 'TopicNames', {'wind'}, 'InstanceID', {0});
    wind_tbl = get_table_safe(wind_raw, 'Wind');
    t_wind_px4 = get_time_sec(wind_tbl);
    
    px4_wn_raw = double(wind_tbl.windspeed_north);
    px4_we_raw = double(wind_tbl.windspeed_east);
    
    px4_wn = interp1(t_wind_px4, px4_wn_raw, sim_time, 'linear', 'extrap');
    px4_we = interp1(t_wind_px4, px4_we_raw, sim_time, 'linear', 'extrap');
else
    px4_wn = zeros(size(sim_time)); px4_we = zeros(size(sim_time));
end
idx_calc = (wx_hist ~= 0); 
t_plot = sim_time(idx_calc) - sim_time(1);
figure('Name', 'Optimization vs PX4', 'Color', 'w');
subplot(2, 1, 1);
plot(t_plot, px4_wn(idx_calc), 'r--', 'LineWidth', 1.5); hold on;
plot(t_plot, wx_hist(idx_calc), 'b', 'LineWidth', 2);
grid on; ylabel('Speed (m/s)'); legend('PX4', 'Optimization');
title('North Wind Estimation');
subplot(2, 1, 2);
plot(t_plot, px4_we(idx_calc), 'r--', 'LineWidth', 1.5); hold on;
plot(t_plot, wy_hist(idx_calc), 'b', 'LineWidth', 2);
grid on; ylabel('Speed (m/s)'); xlabel('Time (s)');
title('East Wind Estimation');
% --- HELPER FUNCTIONS ---
function out = ML_func(x,V_gx,V_gy,V_a,sig_a,sig_g,sig_w,wx_prev,wy_prev)
    P = length(V_gx);
    v_gx_calc = x(1:P);
    v_gy_calc = x(P+1:2*P);
    v_wx = x(2*P + 1);
    v_wy = x(2*P + 2);
    
    V_a_calc = sqrt((v_gx_calc - v_wx).^2 + (v_gy_calc - v_wy).^2);
    
    err_1 = (V_gx - v_gx_calc) / (sqrt(2)*sig_g);
    err_2 = (V_gy - v_gy_calc) / (sqrt(2)*sig_g);
    err_3 = (V_a - V_a_calc)   / (sqrt(2)*sig_a);
    err_4 = (v_wx - wx_prev)   / (sqrt(2)*sig_w);
    err_5 = (v_wy - wy_prev)   / (sqrt(2)*sig_w);
    
    out = [err_1; err_2; err_3; err_4; err_5];
end
function tbl = get_table_safe(raw_topic, name)
    if isempty(raw_topic), error('ERROR: %s Topic is empty.', name); end
    content = raw_topic.TopicMessages;
    if iscell(content), tbl = content{1}; else, tbl = content; end
end
function t_sec = get_time_sec(tbl)
    if istimetable(tbl)
        raw_time = tbl.Properties.RowTimes;
        if isduration(raw_time)
            t_sec = seconds(raw_time);
        else
            t_sec = double(raw_time);
        end
    elseif ismember('timestamp', tbl.Properties.VariableNames)
        t_sec = double(tbl.timestamp) / 1e6;
    elseif ismember('Timestamp', tbl.Properties.VariableNames)
        t_sec = double(tbl.Timestamp) / 1e6;
    else
        error('Time column not found!');
    end
end