%% empiricalFlutterBaseline.m   (B3)
%  Bohon (1966) empirical flutter formula for low-AR swept fins.
%  Compares empirical stability margin against FEM+piston theory result.
%
%  Formula (Bohon 1966, NACA TN 4021, valid for low-AR swept fins):
%    q_flutter_emp = (G_eff * t^2) / (AR_eff^3 * lambda_taper * F_sweep)
%  where:
%    G_eff = D66 * 12 / t^3    [Pa]  effective shear modulus from CLT
%    AR_eff = span^2 / A_plan  effective aspect ratio
%    F_sweep = cos^3(Lambda)
%
%  Note: ±40% accuracy. Use as order-of-magnitude check only.
%
%  Run from the studies/ directory.

clear; clc;
BASE = fullfile(fileparts(mfilename('fullpath')), '..');
addpath(BASE);
fprintf('=== Empirical Flutter Baseline — Bohon (1966) (B3) ===\n\n');

%% --- Laminate and geometry constants ---
lam   = jsondecode(fileread(fullfile(BASE, 'configs', 'lam.json')));
D     = lam.tailored_beta;
t     = lam.flutter_input.t_mm * 1e-3;

cr        = 0.300;
ct        = 0.150;
span      = 0.160;
sweep_deg = 57.4;

A_planform  = (cr + ct) / 2 * span;           % 0.036 m²
AR_eff      = span^2 / A_planform;             % = 0.7111
lambda_taper = ct / cr;                        % = 0.5
F_sweep     = cos(deg2rad(sweep_deg))^3;       % cos^3(57.4°)
G_eff       = D.D66_Nm * 12 / t^3;            % effective shear modulus [Pa]

fprintf('Geometry constants:\n');
fprintf('  A_planform = %.4f m²\n', A_planform);
fprintf('  AR_eff     = %.4f\n',    AR_eff);
fprintf('  lambda     = %.3f\n',    lambda_taper);
fprintf('  F_sweep    = cos³(%.1f°) = %.4f\n', sweep_deg, F_sweep);
fprintf('  G_eff      = D66*12/t³ = %.2f MPa\n', G_eff/1e6);
fprintf('  Formula accuracy: ±40%% (order-of-magnitude only)\n\n');

%% --- Load flight data ---
opts = detectImportOptions(fullfile(BASE, 'data', 'flight_data.csv'), 'CommentStyle', '#');
opts.VariableNames = {'time_s', 'altitude_ft', 'Vz_ms'};
tbl  = readtable(fullfile(BASE, 'data', 'flight_data.csv'), opts);
h_m  = tbl.altitude_ft * 0.3048;
V    = abs(tbl.Vz_ms);

flightPts = struct([]);
for i = 1:height(tbl)
    [rho_i, a_i, ~, ~] = aero.isaAtmosphere(h_m(i));
    M_i = V(i) / a_i;
    % M >= 1.05 ensures beta = sqrt(M^2-1) >= 0.312; piston theory validity per Lighthill (1953) and NACA TN 4021
    if M_i >= 1.05
        fp.Mach=M_i; fp.a=a_i; fp.rho=rho_i;
        fp.U=V(i); fp.q_inf=0.5*rho_i*V(i)^2;
        fp.h_m=h_m(i); fp.time=tbl.time_s(i);
        if isempty(flightPts), flightPts=fp;
        else,                  flightPts(end+1)=fp; %#ok<AGROW>
        end
    end
end
nF = numel(flightPts);
fprintf('Flight points (M >= 1.05): %d\n\n', nF);

%% --- B3.1: Empirical SM for all flight points ---
q_flutter_emp = G_eff * t^2 / (AR_eff^3 * lambda_taper * F_sweep);
SM_emp = zeros(nF, 1);
for fi = 1:nF
    SM_emp(fi) = q_flutter_emp / flightPts(fi).q_inf;
end
mach_vec = [flightPts.Mach];
q_vec    = [flightPts.q_inf];

fprintf('Empirical flutter q = %.2f MPa  (= %.0f × max q_flight)\n', ...
        q_flutter_emp/1e6, q_flutter_emp/max(q_vec));
fprintf('Empirical SM range: %.1f – %.1f\n\n', min(SM_emp), max(SM_emp));

%% --- Load FEM results (if available) ---
femFile = fullfile(BASE, 'results', 'flutter.mat');
fem_loaded = false;
SM_fem = [];
mach_fem = [];
if isfile(femFile)
    R = load(femFile, 'fl_margin', 'flightPts');
    SM_fem   = R.fl_margin;
    mach_fem = [R.flightPts.Mach];
    fem_loaded = true;
    fprintf('FEM results loaded from results/flutter.mat\n\n');
else
    fprintf('Note: results/flutter.mat not found — run mainFlutterSolver.m first for FEM overlay.\n\n');
end

%% --- B3.2: Overlay plot ---
figDir = fullfile(BASE, 'results');
[mach_s, si] = sort(mach_vec);
SM_emp_s = SM_emp(si);

figure('Color','w','Visible','off','Position',[50 50 900 500]);
hold on;

% Empirical curve
semilogy(mach_s, SM_emp_s, 'r-o', 'LineWidth', 1.8, 'MarkerSize', 4, ...
         'DisplayName', 'Bohon (1966) empirical ±40%');

% FEM curve
if fem_loaded
    [mf_s, sfi] = sort(mach_fem);
    SM_fem_s = SM_fem(sfi);
    SM_fem_plot = min(SM_fem_s, 1e4);   % cap Inf for log plot
    semilogy(mf_s, SM_fem_plot, 'b-s', 'LineWidth', 1.8, 'MarkerSize', 4, ...
             'DisplayName', 'FEM (high-fidelity, piston theory)');
    if any(isinf(SM_fem_s))
        text(0.5, 0.97, 'FEM: SM = \infty for all points', 'Units','normalized', ...
             'HorizontalAlignment','center','FontSize',9,'Color','b');
    end
end

% Instability threshold
yline(1, 'k--', 'LineWidth', 1.5, 'DisplayName', 'Instability onset  (SM = 1)');

hold off; grid on; box on;
xlabel('Mach number');
ylabel('Stability margin  q_{flutter} / q_{flight}  (log scale)');
title({'Empirical vs FEM Flutter Stability Margin', ...
       sprintf('Bohon (1966): q_{flutter} = %.1f MPa  (±40%%)', q_flutter_emp/1e6)});
legend('Location','best','FontSize',9);
ylim([0.1, min(max(SM_emp)*2, 1e5)]);

saveas(gcf, fullfile(figDir, 'empirical_vs_fem.png'));
fprintf('Plot saved → results/empirical_vs_fem.png\n');

%% --- B3.3: Interpretation text ---
[min_SM_emp, i_min] = min(SM_emp);
mach_at_min = mach_vec(i_min);

lines = {
    sprintf('Bohon empirical: min SM = %.2f at Mach %.2f', min_SM_emp, mach_at_min), ...
    ''
};

if fem_loaded && all(isinf(SM_fem))
    if min_SM_emp < 1
        interp_str = 'Empirical predicts SM < 1 at some points: FEM analysis was NECESSARY to resolve the design';
    else
        interp_str = sprintf('Empirical predicts SM = %.1f minimum: FEM confirms stability with higher fidelity', min_SM_emp);
    end
    lines{end+1} = 'FEM piston theory: SM = Inf (all points)';
    lines{end+1} = ['Interpretation: ' interp_str];
elseif fem_loaded
    [min_SM_fem, ~] = min(SM_fem(~isinf(SM_fem)));
    if isempty(min_SM_fem), min_SM_fem = Inf; end
    lines{end+1} = sprintf('FEM piston theory: min SM = %.2f', min_SM_fem);
    if min_SM_emp < 1 && min_SM_fem > 1
        lines{end+1} = 'Interpretation: Empirical predicts instability; FEM resolves with higher fidelity — FEM analysis was NECESSARY';
    else
        lines{end+1} = sprintf('Interpretation: Empirical predicts SM = %.1f minimum: FEM confirms and refines with higher confidence', min_SM_emp);
    end
else
    lines{end+1} = 'FEM results not available (run mainFlutterSolver.m first)';
    if min_SM_emp < 1
        lines{end+1} = 'Interpretation: Empirical predicts SM < 1 at some points — FEM analysis is NECESSARY';
    else
        lines{end+1} = sprintf('Interpretation: Empirical predicts SM = %.1f minimum', min_SM_emp);
    end
end

outFile = fullfile(figDir, 'empirical_baseline_notes.txt');
fid = fopen(outFile, 'w');
for k = 1:numel(lines)
    fprintf(fid, '%s\n', lines{k});
end
fclose(fid);

fprintf('\n--- Interpretation ---\n');
for k = 1:numel(lines)
    fprintf('%s\n', lines{k});
end
fprintf('\nNotes saved → results/empirical_baseline_notes.txt\n');

%% --- Save ---
save(fullfile(figDir, 'empirical_baseline.mat'), ...
     'SM_emp', 'mach_vec', 'q_flutter_emp', 'G_eff', 'AR_eff');
fprintf('Data saved → results/empirical_baseline.mat\n');
fprintf('\n=== B3 Done ===\n');
