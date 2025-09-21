%% Full Workflow: Extend Cathode OCP + Infer x0%, x100%, y0%, y100% for Graphite OCP
clear; clc; close all;

%% Load Data
Data_Location=fullfile(fileparts(pwd),"Data");
fileName =fullfile(Data_Location, 'Electrode_Balancing_Data.xlsx');

% Full-cell OCV vs SOC
fullcell = readtable(fileName, 'Sheet', 'SOC_Fullcell', 'VariableNamingRule', 'preserve');
SOC = fullcell.SOC;
OCV_exp = fullcell.OCV;

% Cathode Relative SoL vs OCP
cathode = readtable(fileName, 'Sheet', 'Cathode_Relative', 'VariableNamingRule', 'preserve');
rel_sol = cathode.Relative_SoL;
ocp_rel = cathode.OCP;
sol_abs_start = 0.13; % assumed y_min
sol_abs_end   = 1.00; % assumed y_max
sol_abs = sol_abs_start + (sol_abs_end - sol_abs_start) * rel_sol;

% Graphite Reference OCP
graphite = readtable(fileName, 'Sheet', 'Graphite_Literature', 'VariableNamingRule', 'preserve');
sol_graphite = graphite.SoL;
ocp_graphite_ref = graphite.OCP;
fOCP_graphite_ref = @(x) interp1(sol_graphite, ocp_graphite_ref, x, 'linear', 'extrap');

%% Extend Cathode OCP (MSMR + Linear)
sol_full = linspace(0, 1, 500);
ocp_full = zeros(size(sol_full));

% MSMR fit (low SoL)
x_low = sol_abs(sol_abs < 0.32);
y_low = ocp_rel(sol_abs < 0.32);
msmr_fun = @(p, x) p(1) + p(2)*log((1 - x)./x);
params_low = lsqcurvefit(msmr_fun, [4.2 -0.05], x_low, y_low, [], [], ...
    optimoptions('lsqcurvefit', 'Display', 'off'));

% Linear fit (high SoL)
x_high = sol_abs(sol_abs > 0.91);
y_high = ocp_rel(sol_abs > 0.91);
plin = polyfit(x_high, y_high, 1);

% Build extended cathode OCP
for i = 1:length(sol_full)
    if sol_full(i) < sol_abs_start
        ocp_full(i) = msmr_fun(params_low, sol_full(i));
    elseif sol_full(i) > sol_abs_end
        ocp_full(i) = polyval(plin, sol_full(i));
    else
        ocp_full(i) = interp1(sol_abs, ocp_rel, sol_full(i), 'linear');
    end
end
fOCP_cathode = @(s) interp1(sol_full, ocp_full, s, 'linear', 'extrap');

%% Optimization for Stoichiometric Limits
model_graphite_OCP = @(params, soc) ...
    fOCP_cathode(params(3) + (params(4)-params(3))*(1 - soc)) - OCV_exp;

cost_fun = @(params) model_graphite_OCP(params, SOC) - ...
    fOCP_graphite_ref(params(1) + (params(2)-params(1)) * SOC);

x0 = [0.01, 0.95, 0.1, 0.99];
lb = [0, 0.5, 0, 0.5];
ub = [0.1, 1.0, 0.5, 1.0];

options = optimoptions('lsqnonlin', 'Display', 'iter', ...
    'MaxFunctionEvaluations', 5000, 'MaxIterations', 500);

params_opt = lsqnonlin(cost_fun, x0, lb, ub, options);

x_min = params_opt(1); x_max = params_opt(2);
y_min = params_opt(3); y_max = params_opt(4);

fprintf('\nInferred Stoichiometry Limits:\n');
fprintf('x_min = %.4f\nx_max = %.4f\ny_min = %.4f\ny_max = %.4f\n', ...
    x_min, x_max, y_min, y_max);

%% Reconstruct Graphite OCP
x_SOC = x_min + (x_max - x_min) * SOC;
y_SOC = y_min + (y_max - y_min) * (1 - SOC);
OCP_graphite_inferred = fOCP_cathode(y_SOC) - OCV_exp;

%% Plot Experimental vs Model OCV
figure;
plot(SOC, OCV_exp, 'ko', 'MarkerFaceColor','b'); hold on;
plot(SOC, fOCP_cathode(y_SOC) - OCP_graphite_inferred, 'r-', 'LineWidth', 2);
xlabel('SOC'); ylabel('Voltage (V)');
legend('Experimental OCV','Model OCV'); grid on;

%% Fit MSMR at edges for graphite
% left edge
idx_left = (x_SOC >= x_min) & (x_SOC <= x_min + 0.03);
p_left = lsqcurvefit(msmr_fun, [0.5, -0.5], x_SOC(idx_left), ...
    OCP_graphite_inferred(idx_left), [], [], optimoptions('lsqcurvefit','Display','off'));
% right edge
idx_right = (x_SOC <= x_max) & (x_SOC >= x_max - 0.03);
p_right = lsqcurvefit(msmr_fun, [0.1, 0.1], x_SOC(idx_right), ...
    OCP_graphite_inferred(idx_right), [], [], optimoptions('lsqcurvefit','Display','off'));

% generate full OCP
SoL_full = linspace(0, 1, 500);
ocp_graphite_full = zeros(size(SoL_full));
for i = 1:length(SoL_full)
    x = SoL_full(i);
    if x < x_min
        ocp_graphite_full(i) = msmr_fun(p_left, x);
    elseif x > x_max
        ocp_graphite_full(i) = msmr_fun(p_right, x);
    else
        ocp_graphite_full(i) = interp1(x_SOC, OCP_graphite_inferred, x, 'linear', 'extrap');
    end
end

%% Enforce Monotonicity
ocp_graphite_monotonic = ocp_graphite_full;
for i = 2:length(ocp_graphite_monotonic)
    if ocp_graphite_monotonic(i) >= ocp_graphite_monotonic(i-1)
        ocp_graphite_monotonic(i) = ocp_graphite_monotonic(i-1) - 1e-6;
    end
end

%% Final OCV Validation
SOC_val = linspace(0, 1, 500);
x_SOC_val = x_min + (x_max - x_min) * SOC_val;
y_SOC_val = y_min + (y_max - y_min) * (1 - SOC_val);
OCP_graphite_val = interp1(SoL_full, ocp_graphite_monotonic, x_SOC_val, 'linear', 'extrap');
OCP_cathode_val = fOCP_cathode(y_SOC_val);
OCV_reconstructed = OCP_cathode_val - OCP_graphite_val;

figure;
plot(SOC_val, OCV_reconstructed, 'r-', 'LineWidth', 2); hold on;
plot(SOC, OCV_exp, 'b--');
xlabel('SOC'); ylabel('OCV (V)');
legend('Reconstructed OCV','Experimental OCV'); grid on;

%% Goodness-of-fit metric
rmse = sqrt(mean((OCV_reconstructed - interp1(SOC, OCV_exp, SOC_val, 'linear', 'extrap')).^2));
fprintf('RMSE between reconstructed and experimental OCV = %.4f V\n', rmse);

%% Export Results
T_anode = table(SoL_full(:), ocp_graphite_monotonic(:), ...
    'VariableNames', {'SoL', 'OCP_anode'});
T_cathode = table(sol_full(:), ocp_full(:), ...
    'VariableNames', {'SoL', 'OCP_cathode'});

export_location=fullfile(fileparts(pwd),"Results");

% Export with full path
writetable(T_anode,   fullfile(export_location, "OCP_anode_halfcell.xlsx"));
writetable(T_cathode, fullfile(export_location, "OCP_cathode_halfcell.xlsx"));

fprintf('✅ Exported anode & cathode OCPs to Excel.\n');
