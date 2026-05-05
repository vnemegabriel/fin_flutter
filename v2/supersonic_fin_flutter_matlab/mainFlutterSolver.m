%% mainFlutterSolver.m
%  Supersonic composite fin flutter analysis
%  FEM (Q4 Mindlin shell) + 2nd-order piston theory + p-k flutter solver
%
%  Pipeline:
%    1. Load laminate properties from data/lam.json
%    2. Load flight data from data/flight_data.csv → filter supersonic points
%    3. Build Q4 mesh for swept trapezoidal fin
%    4. Assemble K, M → modal analysis (6 modes)
%    5. Compute GAF at critical Mach via piston theory
%    6. p-k flutter/divergence solve over all supersonic flight points
%    7. Report V_flutter and V_div, plot stability envelope
%
%  Run from the supersonic_fin_flutter_matlab/ directory.

clear; clc; close all;
BASE = fileparts(mfilename('fullpath'));
addpath(fullfile(BASE, 'functions'));
fprintf('=== Supersonic Fin Flutter Solver ===\n\n');

%% -----------------------------------------------------------------------
%% 1. Load laminate (data/lam.json)
%% -----------------------------------------------------------------------
lamFile = fullfile(BASE, 'data', 'lam.json');
if ~isfile(lamFile)
    error('lam.json not found at %s', lamFile);
end
lam  = jsondecode(fileread(lamFile));
D    = lam.tailored_beta;

D26    = 0;
if isfield(D, 'D26_Nm'), D26 = D.D26_Nm; end
D_flex = [D.D11_Nm, D.D12_Nm, D.D16_Nm;
          D.D12_Nm, D.D22_Nm, D26;
          D.D16_Nm, D26,      D.D66_Nm];

t     = lam.flutter_input.t_mm  * 1e-3;     % shell thickness [m]
rho_m = lam.flutter_input.rho_mat_kgm3;     % material density [kg/m³]

% Isotropic-equivalent for membrane/shear DOFs (derived from D66)
E_eff  = 12 * D.D66_Nm / t^3;
nu_eff = 0.3;
geometry.t  = t;
material.E  = E_eff;
material.nu = nu_eff;

fprintf('Laminate (beta=20 tailored, T700/Epoxy AR1):\n');
fprintf('  D11=%.2f  D22=%.2f  D66=%.2f  D16=%.3f  [N·m]\n', ...
        D.D11_Nm, D.D22_Nm, D.D66_Nm, D.D16_Nm);
fprintf('  t=%.2f mm   rho=%.0f kg/m3   E_eff=%.2f GPa\n\n', ...
        t*1e3, rho_m, E_eff/1e9);

%% -----------------------------------------------------------------------
%% 2. Load flight data, filter supersonic points (M >= 1.05)
%% -----------------------------------------------------------------------
csvFile = fullfile(BASE, 'data', 'flight_data.csv');
if ~isfile(csvFile)
    error('flight_data.csv not found at %s', csvFile);
end
opts = detectImportOptions(csvFile, 'CommentStyle', '#');
opts.VariableNames = {'time_s', 'altitude_ft', 'Vz_ms'};
tbl = readtable(csvFile, opts);

h_m = tbl.altitude_ft * 0.3048;
V   = abs(tbl.Vz_ms);

flightPts = struct([]);
for i = 1:height(tbl)
    [rho_i, a_i, ~, ~] = isaAtmosphere(h_m(i));
    M_i = V(i) / a_i;
    if M_i >= 1.05
        fp.Mach  = M_i;   fp.a    = a_i;
        fp.rho   = rho_i; fp.U    = V(i);
        fp.q_inf = 0.5 * rho_i * V(i)^2;
        fp.h_m   = h_m(i); fp.time = tbl.time_s(i);
        if isempty(flightPts), flightPts = fp;
        else, flightPts(end+1) = fp; %#ok<AGROW>
        end
    end
end
if isempty(flightPts)
    error('No supersonic points found in flight_data.csv.');
end

[~, i_crit] = max([flightPts.q_inf]);
fp_crit = flightPts(i_crit);

fprintf('Flight data: %d supersonic points  Mach %.2f-%.2f  h %.0f-%.0f m\n', ...
        numel(flightPts), min([flightPts.Mach]), max([flightPts.Mach]), ...
        min([flightPts.h_m]), max([flightPts.h_m]));
fprintf('Critical point: Mach=%.3f  q=%.0f Pa  h=%.0f m  t=%.1f s\n\n', ...
        fp_crit.Mach, fp_crit.q_inf, fp_crit.h_m, fp_crit.time);

%% -----------------------------------------------------------------------
%% 3. Generate Q4 mesh
%% -----------------------------------------------------------------------
cr        = 0.300;   ct   = 0.150;    % root/tip chord [m]
span      = 0.160;   sweep_deg = 57.4;  % span [m], LE sweep [deg]
nx        = 24;      ny   = 12;

mesh      = GenerarMallaAleta(cr, ct, span, deg2rad(sweep_deg), nx, ny);
rootNodes = find(mesh.nodes(:, 2) < 1e-9);
fixedDOFs = reshape((rootNodes - 1) * 6 + (1:6), 1, []);

fprintf('Mesh: %d nodes, %d elements (%dx%d) | %d root DOFs clamped\n\n', ...
        size(mesh.nodes,1), size(mesh.connect,1), nx, ny, length(fixedDOFs));

%% -----------------------------------------------------------------------
%% 4. FEM assembly and modal analysis
%% -----------------------------------------------------------------------
K = assembleGlobalStiffness(mesh, geometry, material, D_flex);
M = assembleGlobalMass(mesh, rho_m, t);
[K_red, M_red, freeDOFs] = applyDirichletBCs(K, M, fixedDOFs);

nModes = 6;
[Phi_red, omega_n] = modalAnalysis(K_red, M_red, nModes);
f_n = omega_n / (2*pi);

Phi_full = zeros(size(K,1), nModes);
Phi_full(freeDOFs, :) = Phi_red;

fprintf('Natural frequencies [Hz]: ');
fprintf('%.1f  ', f_n); fprintf('\n\n');

%% -----------------------------------------------------------------------
%% 5. GAF at critical flight condition
%% -----------------------------------------------------------------------
% Reference semi-chord (MAC/2) for reduced-frequency k = omega*b_ref/U
b_ref  = (cr^2 + cr*ct + ct^2) / (3*(cr + ct));
k_vals = [0, 0.01, 0.05, 0.1, 0.2, 0.5, 1.0];

fprintf('Computing GAF (Mach=%.3f)... ', fp_crit.Mach);
Q_k      = pistonTheoryGAF(mesh, Phi_full, fp_crit.Mach, fp_crit.q_inf, ...
                            fp_crit.a, k_vals, sweep_deg);
Q_k_norm = Q_k / fp_crit.q_inf;   % per unit dynamic pressure [1/Pa]
fprintf('done\n\n');

%% -----------------------------------------------------------------------
%% 6. p-k flutter/divergence solve
%% -----------------------------------------------------------------------
fprintf('Running p-k solver over %d flight points...\n', numel(flightPts));
[V_fl, V_div, pkRes] = pkSolveFlutter(omega_n, Q_k_norm, k_vals, b_ref, flightPts);
fprintf('done\n\n');

% ── Summary ──────────────────────────────────────────────────────────────
V_flight = [flightPts.U]';
mach_vec = [flightPts.Mach]';

fl_margin  = V_fl  ./ V_flight;
div_margin = V_div ./ V_flight;

[min_fl_margin,  i_fl]  = min(fl_margin);
[min_div_margin, i_div] = min(div_margin);

fprintf('Flutter speed (p-k, g=%.2f):\n', pkRes.g_struct);
if all(isinf(V_fl))
    fprintf('  V_flutter = STABLE for all %d supersonic points\n', numel(flightPts));
    fprintf('  (wash-out design: sweep %.1f deg + D16=%.2f N.m suppresses flutter)\n', ...
            sweep_deg, D_flex(1,3));
else
    fprintf('  Min V_flutter = %.0f m/s  at Mach %.3f  (margin = %.1fx)\n', ...
            V_fl(i_fl), flightPts(i_fl).Mach, min_fl_margin);
end

fprintf('\nDivergence speed (quasi-steady, k=0):\n');
if all(isinf(V_div))
    fprintf('  V_div = STABLE for all %d supersonic points\n\n', numel(flightPts));
else
    fprintf('  Min V_div = %.0f m/s  at Mach %.3f  (margin = %.1fx)\n\n', ...
            V_div(i_div), flightPts(i_div).Mach, min_div_margin);
end

%% -----------------------------------------------------------------------
%% 7. Figures
%% -----------------------------------------------------------------------
figDir = fullfile(BASE, 'results');
[mach_s, si] = sort(mach_vec);

% ── Fig 1: Fin mesh ───────────────────────────────────────────────────────
figure('Color','w','Visible','off');
patch('Faces', mesh.connect, 'Vertices', mesh.nodes(:,1:2), ...
      'FaceColor',[0.88 0.93 1],'EdgeColor',[0.4 0.4 0.4],'LineWidth',0.5);
axis equal tight; grid on; box on;
xlabel('Chordwise x [m]'); ylabel('Spanwise y [m]');
title(sprintf('Q4 Mindlin shell mesh — %d×%d elements\ncr=%.0f mm, ct=%.0f mm, span=%.0f mm, \\Lambda=%.1f°', ...
              nx, ny, cr*1e3, ct*1e3, span*1e3, sweep_deg));
saveas(gcf, fullfile(figDir,'mesh.png'));

% ── Fig 2: Mode shapes (patch, interpolated colour) ───────────────────────
nN    = size(mesh.nodes, 1);
wDOFs = (0:nN-1)' * 6 + 3;       % out-of-plane DOF per node

figure('Color','w','Visible','off','Position',[100 100 900 700]);
for mi = 1:min(4, nModes)
    subplot(2, 2, mi);
    w_vals = Phi_full(wDOFs, mi);
    % Normalise for display
    w_vals = w_vals / max(abs(w_vals) + eps);
    patch('Faces', mesh.connect, 'Vertices', mesh.nodes(:,1:2), ...
          'FaceVertexCData', w_vals, 'FaceColor','interp', ...
          'EdgeColor',[0.3 0.3 0.3], 'EdgeAlpha', 0.25, 'LineWidth', 0.3);
    colormap(gca, 'coolwarm'); clim([-1 1]); colorbar;
    axis equal tight; grid off; box on;
    title(sprintf('Mode %d — f_n = %.1f Hz', mi, f_n(mi)));
    xlabel('x [m]'); ylabel('y [m]');
end
sgtitle('Mass-normalised mode shapes (out-of-plane w, normalised to max=1)');
saveas(gcf, fullfile(figDir,'mode_shapes.png'));

% ── Fig 3: Velocity envelope — V_flight vs V_flutter / V_div ──────────────
time_s  = [flightPts.time]';
V_flt_s = V_fl(si);
V_div_s = V_div(si);
V_flt_s = min(V_flt_s, 5000);   % cap Inf for plotting
V_div_s = min(V_div_s, 5000);

figure('Color','w','Visible','off','Position',[100 100 900 500]);
hold on;
% Safe corridor fill (between flight speed and minimum critical speed)
V_crit = min(V_flt_s, V_div_s);
fill([mach_s; flipud(mach_s)], [V_flight(si); flipud(V_crit)], ...
     [0.85 1 0.85], 'EdgeColor','none', 'FaceAlpha', 0.4, ...
     'DisplayName','Safe corridor');

plot(mach_s, V_flight(si),    'k-',  'LineWidth', 2.0, 'DisplayName','Flight speed');
plot(mach_s, V_flt_s,         'b--', 'LineWidth', 1.8, 'DisplayName','Flutter speed (p-k)');
plot(mach_s, V_div_s,         'r--', 'LineWidth', 1.8, 'DisplayName','Divergence speed');

% Annotate if capped
if any(isinf(V_fl))
    text(0.98, 0.95, 'V_{flutter} = \infty  (wash-out design)', ...
         'Units','normalized','HorizontalAlignment','right', ...
         'Color','b','FontSize',9,'FontWeight','bold');
end
if any(isinf(V_div))
    text(0.98, 0.88, 'V_{div} = \infty  (all div. eigenvalues \leq 0)', ...
         'Units','normalized','HorizontalAlignment','right', ...
         'Color','r','FontSize',9,'FontWeight','bold');
end

hold off; grid on; box on;
xlabel('Mach number  (M \geq 1.05, piston theory valid)');
ylabel('Velocity [m/s]');
title(sprintf(['Aeroelastic velocity envelope  —  \\Lambda=%.1f° swept fin + D_{16}=%.2f N·m\n' ...
               '%d/%d supersonic points flutter-safe   %d/%d divergence-safe'], ...
              sweep_deg, D_flex(1,3), ...
              sum(~isinf(V_fl) & V_fl > V_flight | isinf(V_fl)), numel(flightPts), ...
              sum(~isinf(V_div) & V_div > V_flight | isinf(V_div)), numel(flightPts)));
legend('Location','best','FontSize',9);
saveas(gcf, fullfile(figDir,'velocity_envelope.png'));

% ── Fig 4: p-k damping history at critical flight point ───────────────────
gam_crit = pkRes.gam_hist{i_crit};
q_crit   = pkRes.q_hist{i_crit};

figure('Color','w','Visible','off','Position',[100 100 900 450]);
cmap = lines(nModes);
hold on;
for mi = 1:nModes
    plot(q_crit/1e3, gam_crit(mi,:), '-', 'LineWidth', 1.5, ...
         'Color', cmap(mi,:), 'DisplayName', sprintf('Mode %d (%.0f Hz)', mi, f_n(mi)));
end
yline(0, 'k--', 'LineWidth', 1.5, 'DisplayName', '\gamma = 0  (flutter onset)');
xline(fp_crit.q_inf/1e3, 'k:', 'LineWidth', 1.2, 'DisplayName', 'q_{flight}');
hold off; grid on; box on;
xlabel('Dynamic pressure q  [kPa]');
ylabel('Growth rate \gamma = Re(p)  [rad/s]');
title(sprintf('p-k damping history — Mach=%.3f, h=%.0f m, g_{struct}=%.2f', ...
              fp_crit.Mach, fp_crit.h_m, pkRes.g_struct));
legend('Location','best','FontSize',8);
saveas(gcf, fullfile(figDir,'pk_damping.png'));

%% -----------------------------------------------------------------------
%% 8. Save results
%% -----------------------------------------------------------------------
save(fullfile(figDir, 'flutter.mat'), ...
    'V_fl', 'V_div', 'fl_margin', 'div_margin', ...
    'omega_n', 'f_n', 'Phi_full', 'mesh', ...
    'flightPts', 'fp_crit', 'Q_k', 'Q_k_norm', 'k_vals', 'b_ref', ...
    'pkRes', 'D_flex', 't', 'rho_m');

fprintf('Results saved → results/flutter.mat\n');
fprintf('Plots saved  → results/{mesh, mode_shapes, velocity_envelope, pk_damping}.png\n');
fprintf('\n=== Done ===\n');
