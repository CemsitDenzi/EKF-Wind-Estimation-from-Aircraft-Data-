clc;clear;

if isfile("wind_est_sim_data.mat")
    load("wind_est_sim_data.mat")
else
    error("No file")
end 

total_meas_vg = sim_data.gps.vel;
total_measured_TAS = sim_data.pitot.tas;
sig_g = sim_data.gps.sigma;
sig_a = sim_data.pitot.sigma;
sig_w = 0.5;
N_size = 150; % Window size 
sim_time = sim_data.time;
wx_hist = zeros(size(sim_time));
wy_hist = zeros(size(sim_time));
sf_hist = zeros(size(sim_time));
%initial estimates 
w_x0 = 0;
w_y0 = 0; 

for i = (N_size + 1) :200: size(sim_time,2) - N_size 
    idx = i - N_size : i + N_size;
    meas_vg_x = total_meas_vg(1,idx).';
    meas_vg_y = total_meas_vg(2,idx).';
    meas_TAS = total_measured_TAS(idx).';

    func = @(x) ML_func(x,meas_vg_x,meas_vg_y,meas_TAS,sig_a,sig_g,sig_w,w_x0,w_y0);

    options = optimoptions('lsqnonlin', 'Display', 'none', ...
    'FunctionTolerance', 1e-8, 'StepTolerance', 1e-8);

    x0 = [meas_vg_x;meas_vg_y;w_x0;w_y0];

    x_est = lsqnonlin(func,x0,[],[], options);

    w_x0 = x_est(end-1);
    w_y0 = x_est(end);
    
    wx_hist(i) = w_x0;
    wy_hist(i) = w_y0;

end

idx_calculated = (wx_hist ~= 0); 

figure(1);

subplot(2, 1, 1);
plot(sim_time, repmat(sim_data.truth.wind(1), size(sim_time)), 'r--', 'LineWidth', 2); 
hold on;
plot(sim_time(idx_calculated), wx_hist(idx_calculated), 'b', 'MarkerSize', 8);
grid on;
ylabel('Speed (m/s)');
legend('True Val.', 'Estimated Val.');
title('North Wind');

subplot(2, 1, 2);
plot(sim_time, repmat(sim_data.truth.wind(2), size(sim_time)), 'r--', 'LineWidth', 2); 
hold on;
plot(sim_time(idx_calculated), wy_hist(idx_calculated), 'b', 'MarkerSize', 8);
grid on;
ylabel('Speed (m/s)');
xlabel('Time (s)');
legend('True Val.', 'Estimated Val.');
title('East Wind');

function out = ML_func(x,V_gx,V_gy,V_a,sig_a,sig_g,sig_w,wx_prev,wy_prev)

P = length(V_gx);

v_gx_calc = x(1:P) ;
v_gy_calc = x(P+1:2*P) ;
v_wx = x(2*P + 1);
v_wy = x(2*P + 2);

V_a_calc =    sqrt((v_gx_calc - v_wx).^2 + (v_gy_calc - v_wy).^2);


err_1  = (V_gx - v_gx_calc) / (sqrt(2)*sig_g);
err_2  = (V_gy - v_gy_calc) / (sqrt(2)*sig_g);
err_3 =  (V_a - V_a_calc)/(sqrt(2)*sig_a);
err_4 = (v_wx - wx_prev) / (sqrt(2)*sig_w);
err_5 = (v_wy - wy_prev) / (sqrt(2)*sig_w);


out = [err_1;err_2;err_3;err_4;err_5];

end