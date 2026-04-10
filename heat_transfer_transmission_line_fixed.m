%% =========================================================
%  HEAT TRANSFER IN HIGH-VOLTAGE TRANSMISSION LINES
%  Finite Difference Method (Explicit FTCS Scheme)
%
%  Solves: dT/dt = alpha*d2T/dx2 + Q/rho_cp - h_eff*(T - T_amb)
%
%% =========================================================

clear; clc; close all;

%% =========================================================
%  SECTION 1: PHYSICAL PARAMETERS
%% =========================================================

% --- Conductor properties (ACSR - Aluminium Conductor Steel Reinforced) ---
alpha        = 8e-5;         % Thermal diffusivity            [m^2/s]
rho_cp       = 2.4e6;        % Volumetric heat capacity       [J/m^3/K]
D            = 0.02;         % Conductor diameter             [m]
rho_e        = 2.82e-8;      % Electrical resistivity         [Ohm*m]

% --- Geometry ---
L            = 500;          % Conductor segment length       [m]
T_amb        = 25;           % Ambient temperature            [deg C]
A_cross      = pi*(D/2)^2;   % Cross-sectional area           [m^2]
P_perim      = pi * D;       % Perimeter                      [m]

% --- Convective cooling ---
h_conv       = 20;           % Convective heat transfer coeff [W/m^2/K]
h_eff        = (h_conv * P_perim) / (rho_cp * A_cross);  % Effective cooling [1/s]

% --- Electrical load levels ---
I_levels     = [200, 400, 700, 1000];  % Currents to compare [A]
I_main       = 400;                     % Primary simulation current [A]

% --- Time ---
t_end        = 1200;         % Total simulation time          [s]

%% =========================================================
%  SECTION 2: GRID AND STABILITY
%% =========================================================

N    = 50;
dx   = L / (N - 1);
x    = linspace(0, L, N)';   % Column vector

% Stability-limited time step for FTCS
r_target  = 0.4;
dt_stable = r_target * dx^2 / alpha;


min_steps_desired = 400;
dt_resolution     = t_end / min_steps_desired;
dt                = min(dt_stable, dt_resolution);

r        = alpha * dt / dx^2;
n_steps  = ceil(t_end / dt);

fprintf('=== FDM STABILITY CHECK ===\n');
fprintf('  dx = %.4f m\n', dx);
fprintf('  dt_stable limit = %.6f s\n', dt_stable);
fprintf('  dt used         = %.6f s\n', dt);
fprintf('  r               = %.6f  (must be <= 0.5)\n', r);
fprintf('  Total steps     = %d\n', n_steps);
if r <= 0.5
    fprintf('  Status: STABLE\n\n');
else
    error('Unstable: r = %.4f > 0.5', r);
end

% --- Heat source for primary current ---
Q_vol_main = (I_main^2 * rho_e / A_cross) / rho_cp;
Q_wpm      = I_main^2 * rho_e / A_cross;   % [W/m] for display

fprintf('=== PRIMARY SIMULATION  I = %d A ===\n', I_main);
fprintf('  Heat source  Q     = %.4f W/m\n', Q_wpm);
fprintf('  Heat source  Q_vol = %.6e K/s\n\n', Q_vol_main);

%% =========================================================
%  SECTION 3: PRIMARY FDM SIMULATION
%% =========================================================

snap_times  = linspace(0, t_end, 7);
snap_times  = snap_times(2:end);
T_snaps     = zeros(N, length(snap_times));
snap_count  = 1;

pos_idx     = [floor(N/4), floor(N/2), floor(3*N/4)];

T_temporal_data = zeros(3, n_steps+1);
t_temporal_data = zeros(1, n_steps+1);
rec = 0;

T = T_amb * ones(N, 1);

for step = 0:n_steps
    t_now = min(step * dt, t_end);

    rec = rec + 1;
    t_temporal_data(rec)      = t_now;
    T_temporal_data(:, rec)   = T(pos_idx);

    while snap_count <= length(snap_times) && t_now >= snap_times(snap_count) - dt/2
        T_snaps(:, snap_count) = T;
        snap_count = snap_count + 1;
    end

    if t_now >= t_end
        break;
    end

    T_new          = T;
    T_new(2:N-1)   = T(2:N-1) ...
                   + r * (T(3:N) - 2*T(2:N-1) + T(1:N-2)) ...
                   + dt * (Q_vol_main - h_eff * (T(2:N-1) - T_amb));
    T_new(1) = T_amb;
    T_new(N) = T_amb;
    T        = T_new;
end

t_temporal_data  = t_temporal_data(1:rec);
T_temporal_data  = T_temporal_data(:, 1:rec);
T_final          = T;

while snap_count <= length(snap_times)
    T_snaps(:, snap_count) = T_final;
    snap_count = snap_count + 1;
end

fprintf('Steady-state peak DeltaT: %.4f deg C above ambient\n\n', max(T_final) - T_amb);

%% =========================================================
%  SECTION 4: PLOT 1 — SPATIAL DISTRIBUTION
%% =========================================================

figure('Name','Spatial Temperature Distribution','NumberTitle','off',...
       'Position',[100 100 860 520]);

cmap = lines(length(snap_times));
hold on;
for k = 1:size(T_snaps,2)
    plot(x, T_snaps(:,k) - T_amb, ...
        'Color', cmap(k,:), 'LineWidth', 2, ...
        'DisplayName', sprintf('t = %.0f s', snap_times(k)));
end
hold off;

xlabel('Position along conductor x (m)', 'FontSize', 12);
ylabel('\DeltaT above ambient (°C)',      'FontSize', 12);
title(sprintf('Temperature Distribution — I = %d A', I_main), 'FontSize', 13);
legend('Location','northeastoutside','FontSize',10);
grid on; box on;
set(gca,'FontSize',11);

%% =========================================================
%  SECTION 5: PLOT 2 — TIME EVOLUTION AT THREE POSITIONS
%% =========================================================

figure('Name','Temperature Time Evolution','NumberTitle','off',...
       'Position',[150 150 860 520]);

c3 = [0.10 0.37 0.65; 0.11 0.62 0.46; 0.85 0.35 0.19];
hold on;
for k = 1:3
    lbl = sprintf('x = %.0f m', x(pos_idx(k)));
    plot(t_temporal_data, T_temporal_data(k,:) - T_amb, ...
        'Color', c3(k,:), 'LineWidth', 2, 'DisplayName', lbl);
end
hold off;

xlabel('Time (s)',                    'FontSize', 12);
ylabel('\DeltaT above ambient (°C)', 'FontSize', 12);
title(sprintf('Temperature Evolution at Key Positions — I = %d A', I_main), 'FontSize', 13);
legend('Location','southeast','FontSize',10);
grid on; box on;
set(gca,'FontSize',11);

%% =========================================================
%  SECTION 6: PLOT 3 — CURRENT LEVEL COMPARISON
%% =========================================================

fprintf('=== STEADY-STATE FOR ALL CURRENT LEVELS ===\n');
T_steady   = zeros(N, length(I_levels));
peak_temps = zeros(1, length(I_levels));
Q_sources  = zeros(1, length(I_levels));

for ci = 1:length(I_levels)
    I_c        = I_levels(ci);
    Q_v        = (I_c^2 * rho_e / A_cross) / rho_cp;
    Q_sources(ci) = I_c^2 * rho_e / A_cross;
    T_c        = T_amb * ones(N,1);

    for step = 1:1000000
        T_new_c          = T_c;
        T_new_c(2:N-1)   = T_c(2:N-1) ...
                         + r*(T_c(3:N) - 2*T_c(2:N-1) + T_c(1:N-2)) ...
                         + dt*(Q_v - h_eff*(T_c(2:N-1) - T_amb));
        T_new_c(1) = T_amb;
        T_new_c(N) = T_amb;

        if max(abs(T_new_c - T_c)) < 1e-8
            T_c = T_new_c;
            break;
        end
        T_c = T_new_c;
    end

    T_steady(:,ci)  = T_c;
    peak_temps(ci)  = max(T_c) - T_amb;
    fprintf('  I = %4d A  | Peak DeltaT = %7.3f C  | Q = %.4f W/m\n', ...
            I_c, peak_temps(ci), Q_sources(ci));
end

figure('Name','Current Level Comparison','NumberTitle','off',...
       'Position',[200 200 860 520]);

c4     = [0.10 0.37 0.65; 0.11 0.62 0.46; 0.85 0.35 0.19; 0.64 0.08 0.18];
lstyle = {'-','--','-.',':'};
hold on;
for ci = 1:length(I_levels)
    plot(x, T_steady(:,ci) - T_amb, ...
        'Color', c4(ci,:), 'LineWidth', 2.2, 'LineStyle', lstyle{ci}, ...
        'DisplayName', sprintf('I = %d A  (peak %.2f°C)', I_levels(ci), peak_temps(ci)));
end
hold off;

xlabel('Position along conductor x (m)',        'FontSize', 12);
ylabel('Steady-state \DeltaT above ambient (°C)','FontSize', 12);
title('Steady-State Profile — Current Level Comparison', 'FontSize', 13);
legend('Location','north','FontSize',10);
grid on; box on;
set(gca,'FontSize',11);

%% =========================================================
%  SECTION 7: PLOT 4 — PEAK TEMPERATURE vs CURRENT
%% =========================================================

figure('Name','Peak Temp vs Current','NumberTitle','off',...
       'Position',[250 250 720 460]);

I_dense = linspace(100, 1200, 300);
T_analytical = (I_dense.^2 * rho_e / A_cross) ./ (h_conv * P_perim / A_cross);

yyaxis left
plot(I_dense, T_analytical, 'LineWidth',2,...
     'DisplayName','Analytical estimate (I^2 law)');
hold on;
scatter(I_levels, peak_temps, 90, 'filled',...
        'DisplayName','FDM simulation peaks');
hold off;
ylabel('\DeltaT above ambient (°C)','FontSize',12);

yyaxis right
plot(I_dense, I_dense.^2 * rho_e / A_cross, '--',...
     'LineWidth',1.5,...
     'DisplayName','Heat source Q (W/m)');
ylabel('Heat generation Q (W/m)','FontSize',12);

xlabel('Current I (A)','FontSize',12);
title('Peak Temperature Rise and Heat Source vs Current','FontSize',13);
legend('Location','northwest','FontSize',10);
grid on; box on;
set(gca,'FontSize',11);

%% =========================================================
%  SECTION 8: SUMMARY TABLE
%% =========================================================

fprintf('\n=== SUMMARY TABLE ===\n');
fprintf('%-14s %-14s %-18s %-16s %-12s\n',...
        'Current (A)','Q (W/m)','Peak DeltaT (C)','Ratio to 200A','I^2 ratio');
fprintf('%s\n', repmat('-',1,76));
for ci = 1:length(I_levels)
    fprintf('%-14d %-14.4f %-18.4f %-16.4f %-12.4f\n',...
            I_levels(ci), Q_sources(ci), peak_temps(ci),...
            peak_temps(ci)/peak_temps(1), (I_levels(ci)/I_levels(1))^2);
end
fprintf('\nDeltaT ratio should closely match I^2 ratio — confirms Ohmic heating law.\n');