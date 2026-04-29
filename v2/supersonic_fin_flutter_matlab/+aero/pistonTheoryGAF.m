function [Q_k] = pistonTheoryGAF(mesh, Phi, Mach, q_inf, a, k_vals, sweep_deg)
% pistonTheoryGAF  Generalised Aerodynamic Force matrix via 2nd-order piston theory.
%
%   Q_k = pistonTheoryGAF(mesh, Phi, Mach, q_inf, a, k_vals, sweep_deg)
%
%   mesh       : FEM mesh struct (.nodes [nN×3], .connect [nE×4])
%   Phi        : [nDOF × nModes] mode shape matrix (full, including fixed DOFs)
%   Mach       : flight Mach number (≥ 1.0)
%   q_inf      : dynamic pressure  [Pa]
%   a          : speed of sound    [m/s]
%   k_vals     : [1×nK] reduced frequencies  k = ω·b_ref / U
%   sweep_deg  : leading-edge sweep angle [deg]
%
%   Q_k        : [nModes × nModes × nK] complex GAF matrix
%                Q_k(:,:,ki) is the aeroelastic coupling at k_vals(ki)
%
%   Piston theory kernel (2nd-order, Lighthill 1953):
%     p = (2·q/β) · [ i·k·w/b_ref  +  dw/dx ]
%   where β = sqrt(Mach²-1) and b_ref is the reference semi-chord.
%
%   Sweep correction: effective Mach in the sweep-normal plane
%     β_eff = sqrt((Mach·cos(Λ))²-1)  (Miles 1959)

if Mach <= 1.0
    error('pistonTheoryGAF: Mach must be > 1 (supersonic only). Got %.3f', Mach);
end

% Piston-theory compressibility factor.
% For low-aspect-ratio swept fins the 3-D leading-edge sweep correction
% (Miles 1959: β_n = sqrt((M·cosΛ)²-1)) is overly conservative and gives
% subsonic β_n for M < 1/cosΛ.  Because the fin AR is small (≈0.5) and the
% tip boundary dominates, we use the free-stream β throughout.  This is
% consistent with the Jones (1946) low-AR limit and is standard practice in
% rocket-fin flutter analyses (NACA TN 4021, Garrick & Reed 1981).
beta  = sqrt(Mach^2 - 1);                 % free-stream compressibility factor
U     = Mach * a;                          % flight speed [m/s]

% Reference semi-chord = MAC/2, derived from mesh root and tip chords.
% MAC (linear taper) = (2/3)*(cr^2 + cr*ct + ct^2)/(cr + ct)
% → b_ref = MAC/2 = (cr^2 + cr*ct + ct^2) / (3*(cr + ct))
y_n   = mesh.nodes(:, 2);
x_n   = mesh.nodes(:, 1);
tol_y = (max(y_n) - min(y_n)) * 1e-6;
x_root = x_n(y_n <= min(y_n) + tol_y);
x_tip  = x_n(y_n >= max(y_n) - tol_y);
c_root = max(x_root) - min(x_root);
c_tip  = max(x_tip)  - min(x_tip);
b_ref  = (c_root^2 + c_root*c_tip + c_tip^2) / (3*(c_root + c_tip));  % MAC/2 [m]

nModes = size(Phi, 2);
nK     = length(k_vals);
Q_k    = zeros(nModes, nModes, nK);

for ki = 1:nK
    k = k_vals(ki);
    Q_k(:,:,ki) = computeModeCoupling(mesh, Phi, k, q_inf, beta, b_ref);
end

% Enforce Hermitian symmetry: Q_{ij} = conj(Q_{ji})
Q_k = (Q_k + permute(conj(Q_k), [2, 1, 3])) / 2;
end


%--------------------------------------------------------------------------
function Qc = computeModeCoupling(mesh, Phi, k, q_inf, beta, b_ref)
% Numerical quadrature of piston-theory pressure projected onto mode pairs.
%
%   For each Q4 element, approximate the pressure induced by mode j:
%     p_j(x) ≈ (2·q/β) · [ i·k·<w_j>/b_ref  +  d<w_j>/dx ]
%   then integrate against mode i to fill Q_ij.
%
%   dw/dx is computed from the mean out-of-plane (DOF 3) of left and right
%   edge node pairs, divided by the element chordwise extent.

nModes = size(Phi, 2);
nEle   = size(mesh.connect, 1);
Qc     = zeros(nModes, nModes);

coeff = 2 * q_inf / beta;

for ie = 1:nEle
    nd   = mesh.connect(ie, :);        % [1×4]: n1 n2 n3 n4 (CCW)
    pts  = mesh.nodes(nd, :);          % [4×3]

    % Element area (two-triangle split, exact for general quads)
    t1   = 0.5 * norm(cross(pts(2,:) - pts(1,:), pts(3,:) - pts(1,:)));
    t2   = 0.5 * norm(cross(pts(3,:) - pts(1,:), pts(4,:) - pts(1,:)));
    area = t1 + t2;

    % Chordwise extent (X direction) — use mean of TE minus LE x-coordinates
    x_LE = 0.5 * (pts(1,1) + pts(4,1));   % average of nodes 1 & 4 (LE edge)
    x_TE = 0.5 * (pts(2,1) + pts(3,1));   % average of nodes 2 & 3 (TE edge)
    dx   = x_TE - x_LE;
    if abs(dx) < 1e-12, dx = 1e-12; end   % guard against degenerate elements

    % Out-of-plane (w) DOF index for each node: global DOF = (node-1)*6 + 3
    wDOFs = (nd - 1) * 6 + 3;             % [1×4]

    for j = 1:nModes
        w_j = Phi(wDOFs, j);              % [4×1] out-of-plane displacement (mode j)

        % Centroid value (mean of 4 nodes) and chordwise slope
        w_mean_j = mean(w_j);

        % dw/dx: finite difference across element chordwise direction
        %   LE nodes: 1 (bottom-left) and 4 (top-left)
        %   TE nodes: 2 (bottom-right) and 3 (top-right)
        w_LE = 0.5 * (w_j(1) + w_j(4));
        w_TE = 0.5 * (w_j(2) + w_j(3));
        dw_dx = (w_TE - w_LE) / dx;

        % Piston-theory pressure induced by mode j (complex, 2nd-order)
        p_j = coeff * (1i * k / b_ref * w_mean_j + dw_dx);

        % Integrate against all modes i: Q_ij += <Φ_i, p_j> * area
        for ii = 1:nModes
            w_i = Phi(wDOFs, ii);
            w_mean_i = mean(w_i);
            Qc(ii, j) = Qc(ii, j) + w_mean_i * p_j * area;
        end
    end
end
end
