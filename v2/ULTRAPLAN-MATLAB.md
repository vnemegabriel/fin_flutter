
## 📁 1. PROJECT DIRECTORY STRUCTURE
```
supersonic_fin_flutter_matlab/
├── +core/
│   ├── getTransformationMatrix.m      # PROVIDED: Unchanged
│   ├── CalcularRigidezQLLL.m          # PROVIDED: Modified internally for anisotropy
│   └── CalcularEsfuerzosCompletos.m   # PROVIDED: Unchanged
├── +fem/
│   ├── assembleGlobalStiffness.m
│   ├── assembleGlobalMass.m
│   ├── applyDirichletBCs.m
│   └── modalAnalysis.m
├── +aero/
│   ├── pistonTheoryGAF.m
│   └── isaAtmosphere.m
├── +stability/
│   ├── loewnerInterpolation.m
│   └── solveFlutterPL.m
├── configs/
│   ├── lam.json                       # PROVIDED INPUT
│   └── flight_data.csv                # PROVIDED INPUT
├── tests/
│   ├── testFEMAssembly.m
│   └── testPLSolver.m
└── mainFlutterSolver.m                # ORCHESTRATOR
```

---

## 🔹 2. INPUT DATA SPECIFICATIONS

### `configs/lam.json`
- **Path:** `configs/lam.json`
- **Load command:** `lam = jsondecode(fileread('configs/lam.json'));`
- **Key fields to extract:**
  ```matlab
  D = lam.tailored_beta;
  D_flex = [D.D11_Nm, D.D12_Nm, 0, D.D16_Nm;
            D.D12_Nm, D.D22_Nm, 0, 0;
            0,        0,        D.D66_Nm, 0;
            D.D16_Nm, 0,        0, D.D66_Nm];
  t = lam.flutter_input.t_mm * 1e-3; % [m]
  ```

### `configs/flight_data.csv`
- **Columns (after skipping `#` comments):** `time_s`, `altitude_ft`, `Vz_ms`
- **Processing:** 
  1. Convert `altitude_ft` → `h_m = altitude_ft * 0.3048`
  2. Velocity: `V = abs(Vz_ms)`
  3. ISA model → `rho, a, T, P`
  4. `Mach = V ./ a; q_inf = 0.5 * rho .* V.^2;`
  5. Filter: `Mach >= 1.0` (supersonic regime only)
  6. Return struct array with fields: `time, alt_m, V, Mach, q_inf, rho, a`

---

## 🔹 3. CORE FEM IMPLEMENTATION (`+core/`)

### `+core/getTransformationMatrix.m`
- **Use EXACTLY as provided in `getTransformationMatrix.txt`**
- Signature: `function [T, nodes2D] = getTransformationMatrix(nodes3D)`
- `nodes3D`: `[4x3]` double. Returns `T` `[24x24]`, `nodes2D` `[4x2]`.

### `+core/CalcularRigidezQLLL.m` (MODIFIED FOR ANISOTROPY)
- **Base:** Use `CalcularRigidezQLLL.txt` exactly.
- **Modification:** Replace the isotropic flexion matrix construction with anisotropic input.
- **New Signature:** `function [Ke_global] = CalcularRigidezQLLL(nodes3D, geometry, D_flex_4x4, integrationType)`
- **Internal change (line ~15-20):**
  ```matlab
  % REPLACE:
  % Db = (E*t^3/(12*(1-nu^2))) * [1, nu, 0; nu, 1, 0; 0, 0, (1-nu)/2];
  % WITH:
  Db = D_flex_4x4(1:3, 1:3); % Use [11,22,12; 22,22,12; 12,12,66] block
  Dm = (geometry.t) * [lam_eff(1), lam_eff(3), 0; 
                       lam_eff(3), lam_eff(2), 0; 
                       0, 0, lam_eff(4)]; % Membrane (isotropic equivalent for now)
  ```
- **Keep:** Selective integration (`nGaussS=1`), drilling stiffness (`k_drill=1e-3*E*t*Area_elem`), transformation `Ke_global = T' * Ke_local * T`.
- **Note:** For membrane, use equivalent `E_eff = 12*D.D66_Nm/t^3`, `nu_eff=0.3` to maintain compatibility. Pass `geometry` with `.t` and `.E_eff`.

### `+core/CalcularEsfuerzosCompletos.m`
- **Use EXACTLY as provided in `CalcularEsfuerzosCompletos .txt`**
- No modifications needed for flutter pipeline.

---

## 🔹 4. ASSEMBLY & MODAL ANALYSIS (`+fem/`)

### `+fem/assembleGlobalStiffness.m`
```matlab
function [K] = assembleGlobalStiffness(mesh, geometry, D_flex_4x4)
nNodes = size(mesh.nodes,1); nDOF = 6*nNodes;
K = sparse(nDOF, nDOF);
nEle = size(mesh.connect,1);
for i=1:nEle
    nodes3D = mesh.nodes(mesh.connect(i,:), :);
    Ke = CalcularRigidezQLLL(nodes3D, geometry, D_flex_4x4, 'selective');
    idx = reshape((mesh.connect(i,:)'-1)*6 + (1:6)', [], 1);
    K(idx,idx) = K(idx,idx) + Ke;
end
end
```

### `+fem/assembleGlobalMass.m`
```matlab
function [M] = assembleGlobalMass(mesh, rho, t)
nNodes = size(mesh.nodes,1); nDOF = 6*nNodes;
M = sparse(nDOF, nDOF);
nEle = size(mesh.connect,1);
for i=1:nEle
    nodes3D = mesh.nodes(mesh.connect(i,:), :);
    area = 0.5 * norm(cross(nodes3D(2,:)-nodes3D(1,:), nodes3D(4,:)-nodes3D(1,:)));
    me = (rho*t*area/180) * ... % Consistent mass matrix scaling (simplified for 4-node Q4)
         [2 0 0 1 0 0; 0 2 0 0 1 0; 0 0 2 0 0 1; 1 0 0 2 0 0; 0 1 0 0 2 0; 0 0 1 0 0 2]; % 3x3 block per node
    % Actually, use lumped mass for stability: M_ii = rho*t*area/4 for each DOF
    M_local = diag(ones(1,24) * (rho*t*area/4));
    idx = reshape((mesh.connect(i,:)'-1)*6 + (1:6)', [], 1);
    M(idx,idx) = M(idx,idx) + M_local;
end
end
```

### `+fem/modalAnalysis.m`
```matlab
function [Phi, omega_n] = modalAnalysis(K, M, nModes)
% Solve (K - w^2 M) Phi = 0
opts.disp = 0; opts.issym = true; opts.isreal = true;
[Phi, Omega2] = eigs(K, M, nModes, 'sm', opts);
omega_n = sqrt(diag(Omega2));
% Sort by frequency
[omega_n, idx] = sort(omega_n);
Phi = Phi(:, idx);
end
```

### `+fem/applyDirichletBCs.m`
```matlab
function [K_red, M_red, freeDOFs] = applyDirichletBCs(K, M, fixedDOFs)
allDOFs = 1:size(K,1);
freeDOFs = setdiff(allDOFs, fixedDOFs)';
K_red = K(freeDOFs, freeDOFs);
M_red = M(freeDOFs, freeDOFs);
end
```

---

## 🔹 5. AERODYNAMIC INTERFACE (`+aero/`)

### `+aero/isaAtmosphere.m`
- Standard ISA up to 11 km, then isothermal to 20 km, lapse to 32 km.
- Returns `[rho, a]` for given `h_m`.
- Implement standard formulas. Ensure vectorization.

### `+aero/pistonTheoryGAF.m`
```matlab
function [Q_k] = pistonTheoryGAF(mesh, Phi, Mach, q_inf, rho, k_vals, sweep_deg)
% 2nd order supersonic piston theory
% k = omega*b_ref / U. b_ref = mean aerodynamic chord
beta = sqrt(Mach^2 - 1);
U = Mach * sqrt(1.4*287.287); % Approximate a0*T0 correction, or pass a explicitly
b_ref = 0.22; % Placeholder mean chord
omega_vals = k_vals * U / b_ref;

nModes = size(Phi,2);
nK = length(k_vals);
Q_k = zeros(nModes, nModes, nK);

% Simplified strip integration over element faces
% For each k, compute complex aerodynamic influence
for ki = 1:nK
    k = k_vals(ki);
    % Piston theory kernel: P = (2*q/beta) * (ik * w/U + dw/dx)
    % Project onto modes: Q_ij = int(Phi_i * P(Phi_j) dA)
    % Use numerical quadrature over mesh
    Q_k(:,:,ki) = computeModeCoupling(mesh, Phi, k, Mach, q_inf, U, beta);
end
% Enforce Hermitian symmetry for stability
Q_k = (Q_k + permute(conj(Q_k), [2,1,3])) / 2;
end

% Subfunction: computeModeCoupling
function Qc = computeModeCoupling(mesh, Phi, k, Mach, q_inf, U, beta)
% Implement Gauss integration over elements
% Use local derivatives from getTransformationMatrix
% Return complex nModes x nModes matrix
% Keep it concise: assume linear pressure variation over Q4 element
% Use provided shape function derivatives
end
```

---

## 🔹 6. p-L STABILITY SOLVER (`+stability/`)

### `+stability/loewnerInterpolation.m`
```matlab
function [L, M] = loewnerInterpolation(Q_k, s_vals)
% s_vals = 1i * omega_vals
nK = length(s_vals); nModes = size(Q_k,1);
L = zeros(nModes*nK, nModes*nK);
M = zeros(nModes*nK, nModes*nK);
for i=1:nK
    for j=1:nK
        if i~=j
            ds = s_vals(i) - s_vals(j);
            L((i-1)*nModes+1:i*nModes, (j-1)*nModes+1:j*nModes) = (Q_k(:,:,i) - Q_k(:,:,j)) / ds;
            M((i-1)*nModes+1:i*nModes, (j-1)*nModes+1:j*nModes) = (s_vals(i)*Q_k(:,:,i) - s_vals(j)*Q_k(:,:,j)) / ds;
        end
    end
end
% Regularize near-zero denominators
L(abs(L)<1e-12) = 1e-12;
M(abs(M)<1e-12) = 1e-12;
end
```

### `+stability/solveFlutterPL.m`
```matlab
function [V_fl, V_div, results] = solveFlutterPL(K_red, M_red, Phi, Q_k, k_vals, flightConditions)
% Build state-space for each flight point
nF = length(flightConditions.Mach);
V_fl = zeros(nF,1); V_div = zeros(nF,1);
results = struct();

for i=1:nF
    s_vals = 1i * k_vals * flightConditions.U(i) / 0.22; % b_ref=0.22
    [L, M_loe] = loewnerInterpolation(Q_k, s_vals);
    % Assemble full aeroelastic system: (s^2*M + K - q*Q(s))x = 0
    % State-space realization: A = [0 I; -M\K 0] - q*[0 0; 0 M\Q_loe]
    % Solve generalized eigenvalue problem (A - lambda E)x = 0
    % Track max(real(lambda)) vs q_inf
    % Interpolate flutter (real=0) and divergence (imag=0)
    % Output V_fl(i), V_div(i)
end
end
```

---

## 🔹 7. MAIN ORCHESTRATOR (`mainFlutterSolver.m`)
- **Geometry:** `cr=0.300, ct=0.150, span=0.200, sweep=57.4°, nx=24, ny=12`
- **Mesh Generation:** Structured Q4 grid. Map `(xi,eta)` to physical swept coordinates.
- **BCs:** Root clamped (all 6 DOFs = 0 at `y=0` nodes)
- **Flow:**
  1. Load `lam.json` → `D_flex`, `t`, `rho_mat`
  2. Load `flight_data.csv` → filter `Mach>=1`
  3. Generate mesh → apply BCs → assemble `K, M` → modal extract `Phi` (6 modes)
  4. Loop over flight conditions → compute `Q_k(k)` → run `solveFlutterPL`
  5. Plot `V_flight vs V_flutter`, `damping vs altitude`, `mode shapes`
  6. Save `results/flutter.mat`

---

## ⚠️ 8. STRICT IMPLEMENTATION RULES
1. **NO TOOLBOXES:** Use only base MATLAB. Replace `eig`/`eigs` if unavailable (implement Arnoldi fallback if needed, but `eigs` is base).
2. **INDEXING:** MATLAB is 1-based. All DOF mappings must use `6*(node-1)+dof`.
3. **SPARSE MATRICES:** Always initialize `K = sparse(nDOF, nDOF)` and `M = sparse(...)`.
4. **COMPLEX ARITHMETIC:** `Q_k` must be `complex`. Use `conj()` and `permute` for symmetry.
5. **ERROR HANDLING:** Wrap `eigs` in `try/catch`. If it fails, fall back to `eig(full(K), full(M))`.
6. **DIMENSION CHECKS:** Add `assert(size(Phi,1) == size(K,1))` and similar checks.
7. **FILE I/O:** Use `jsondecode` and `readtable` with `CommentStyle','#'`.

---

## ✅ 9. VALIDATION CHECKPOINTS
| Step | Test | Pass Criteria |
|------|------|---------------|
| 1. Assembly | `eig(K)` on cantilever beam mesh | First 6 eigenvalues ≈ 0 (rigid body), rest > 0 |
| 2. Modal | `modalAnalysis` vs analytical clamped plate | Error < 5% in f1, f2 |
| 3. Aero | `Q_k(:,:,1)` at k=0 | Real, symmetric, matches static piston limit |
| 4. p-L | `V_fl` vs `q_inf` sweep | Single crossing of `real(λ)=0` |
| 5. Full | Run `mainFlutterSolver` | No crashes, generates `results/` folder, plots saved |
