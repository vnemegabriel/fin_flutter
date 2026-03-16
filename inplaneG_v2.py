import math
import itertools

# =============================================================================
# MATERIAL DATABASE
# Derived from manufacturer datasheets and classical micromechanics.
# T700-grade 12K carbon fiber assumed for both fabrics.
# =============================================================================

# --- Fiber (T700 carbon, typical) ---
Ef1  = 230e9   # Fiber longitudinal modulus [Pa]
Ef2  = 15e9    # Fiber transverse modulus [Pa]  (for Halpin-Tsai accuracy)
Gf12 = 27e9    # Fiber shear modulus [Pa]
nuf  = 0.20    # Fiber Poisson ratio

# --- Matrix (Epoxy LY1564 / Ampreg 22 class, Vf-compatible) ---
Em  = 3.5e9    # Matrix modulus [Pa]
num = 0.35     # Matrix Poisson ratio
Gm  = Em / (2.0 * (1.0 + num))  # Matrix shear modulus [Pa]

# --- Fin geometry ---
cr_m = 0.300   # Root chord [m]
ct_m = 0.150   # Tip chord  [m]
s_m  = 0.160   # Exposed span (height) [m]
t_m  = 0.006   # Total fin thickness [m]
sweep_deg = 57.4  # Leading-edge sweep angle [deg]

# =============================================================================
# FABRIC 1 — CARBONODB300  (±45° biaxial NCF, 300 g/m²)
# Architecture: non-crimp, stitched biaxial. Two equal layers: +45 and -45.
# Areal weight (FAW): 300 g/m²  → each ±45 ply ~150 g/m²
# Nominal ply thickness at Vf=0.50:
#   t = FAW / (Vf * rho_fiber)  ;  rho_fiber(T700) = 1800 kg/m³
# FAW_total = 0.300 kg/m²  →  t_biax_pair = 0.300 / (0.50 * 1800) = 0.333 mm per PAIR
# So each single ±45 ply: t_single = 0.167 mm
# =============================================================================
FAW_DB300 = 0.300   # kg/m²
rho_fiber = 1800.0  # kg/m³  (T700)
Vf_DB300  = 0.50
t_DB300_pair   = FAW_DB300 / (Vf_DB300 * rho_fiber)   # thickness for the full ±45 pair [m]
t_DB300_single = t_DB300_pair / 2.0                    # thickness per individual ply [m]

# =============================================================================
# FABRIC 2 — CARBONOGA90R  (UD 0°/90°, 302 g/m², 12K)
# Architecture: true woven 0/90.
# FAW_total = 302 g/m² → equal 0° and 90° tows → each direction ~151 g/m²
# t_woven = FAW / (Vf * rho_fiber) = 0.302 / (0.50 * 1800) = 0.336 mm
# =============================================================================
FAW_GA90R = 0.302   # kg/m²
Vf_GA90R  = 0.50
t_GA90R   = FAW_GA90R / (Vf_GA90R * rho_fiber)        # thickness per woven ply [m]

# =============================================================================
# MICROMECHANICS
# All moduli computed via Halpin-Tsai with validated xi values.
# =============================================================================

def halpin_tsai(Ep, Em, Vf, xi):
    """Halpin-Tsai modulus in matrix-dominated direction."""
    eta = ((Ep / Em) - 1.0) / ((Ep / Em) + xi)
    return Em * (1.0 + xi * eta * Vf) / (1.0 - eta * Vf)

def ud_properties(Ef1, Ef2, Gf12, nuf, Em, num, Gm, Vf):
    """Full UD ply elastic constants via rule-of-mixtures + Halpin-Tsai."""
    E1  = Ef1 * Vf + Em * (1.0 - Vf)                  # ROM (fiber-dominated)
    nu12 = nuf * Vf + num * (1.0 - Vf)                 # ROM
    E2  = halpin_tsai(Ef2, Em, Vf, xi=2.0)             # H-T, xi=2 (elliptical fibers)
    G12 = halpin_tsai(Gf12, Gm, Vf, xi=1.0)            # H-T, xi=1 (shear)
    return E1, E2, G12, nu12

# --- Compute UD building-block ---
E1_ud, E2_ud, G12_ud, nu12_ud = ud_properties(Ef1, Ef2, Gf12, nuf, Em, num, Gm, Vf_DB300)

# --- FABRIC 1: ±45° NCF biaxial (CARBONODB300) ---
# A ±45 NCF has zero crimp (non-crimp fabric). 
# Properties derived by exact CLPT averaging of +45 and -45 UD plies.
# For a balanced ±45 laminate: Ex = Ey = E_45, Gxy is maximized.
# Using the standard CLPT result for [+45/-45]s:
nu21_ud = nu12_ud * E2_ud / E1_ud
D_ud = 1.0 - nu12_ud * nu21_ud
Q11 = E1_ud / D_ud
Q22 = E2_ud / D_ud
Q12 = nu12_ud * E2_ud / D_ud
Q66 = G12_ud

# Transformed Qbar at 45 degrees
m45 = math.cos(math.radians(45))
n45 = math.sin(math.radians(45))
# For [+45/-45]: Ex = (Q11+Q22+2Q12+4Q66)/4 - 2*(off-axis coupling terms cancel)
# Exact formula for Ex of balanced ±45 laminate:
Ex_45lam = (Q11 * m45**4 + Q22 * n45**4 + 2*(Q12 + 2*Q66)*n45**2*m45**2 +
            Q11 * n45**4 + Q22 * m45**4 + 2*(Q12 + 2*Q66)*n45**2*m45**2) / 2.0
Gxy_45lam = ((Q11 + Q22 - 2*Q12)*n45**2*m45**2 + Q66*(n45**4 + m45**4))  # Qbar66 at 45°

E1_DB300 = Ex_45lam          # Axial stiffness of ±45 biaxial
E2_DB300 = Ex_45lam          # Balanced: same in both directions
G12_DB300 = Gxy_45lam        # Shear — this is MAXIMUM for ±45 layup
nu12_DB300 = 0.70            # ±45 laminates have high Poisson ratio
t_DB300 = t_DB300_single     # Use per-ply thickness in stacking

# --- FABRIC 2: Woven 0/90 (CARBONOGA90R) ---
# True woven: apply crimp knockdown (kc = 0.92 for plain weave carbon)
kc = 0.92
E1_GA90R  = ((E1_ud + E2_ud) / 2.0) * kc    # Averaged, crimp-reduced
E2_GA90R  = E1_GA90R                          # Balanced woven
G12_GA90R = halpin_tsai(Gf12, Gm, Vf_GA90R, xi=1.0)  # Shear: no crimp penalty
nu12_GA90R = 0.05                             # Balanced woven: low nu
t_GA90R_ply = t_GA90R                         # Per ply thickness

# =============================================================================
# PRINT MATERIAL SUMMARY
# =============================================================================
print("=" * 65)
print("  MATERIAL PROPERTIES — Derived from Specified Fabrics")
print("=" * 65)

print(f"\n  UD Building Block (T700/Epoxy, Vf={Vf_DB300:.0%})")
print(f"    E1  = {E1_ud/1e9:.2f} GPa  |  E2  = {E2_ud/1e9:.2f} GPa")
print(f"    G12 = {G12_ud/1e9:.2f} GPa  |  nu12 = {nu12_ud:.3f}")

print(f"\n  CARBONODB300 — ±45° NCF Biaxial (FAW=300 g/m²)")
print(f"    Architecture : Non-crimp, zero crimp penalty")
print(f"    t_pair = {t_DB300_pair*1e3:.3f} mm  |  t_ply = {t_DB300_single*1e3:.3f} mm")
print(f"    E1  = {E1_DB300/1e9:.2f} GPa  |  E2  = {E2_DB300/1e9:.2f} GPa")
print(f"    G12 = {G12_DB300/1e9:.2f} GPa  |  nu12 = {nu12_DB300:.2f}")

print(f"\n  CARBONOGA90R — Woven 0/90 UD (FAW=302 g/m², 12K)")
print(f"    Architecture : Plain woven, crimp kc={kc}")
print(f"    t_ply = {t_GA90R_ply*1e3:.3f} mm")
print(f"    E1  = {E1_GA90R/1e9:.2f} GPa  |  E2  = {E2_GA90R/1e9:.2f} GPa")
print(f"    G12 = {G12_GA90R/1e9:.2f} GPa  |  nu12 = {nu12_GA90R:.2f}")

# =============================================================================
# CLPT STIFFNESS MATRICES
# Plies defined as tuples: (angle_deg, E1, E2, G12, nu12, t)
# =============================================================================

def Q_matrix(E1, E2, G12, nu12):
    nu21 = nu12 * E2 / E1
    D = 1.0 - nu12 * nu21
    return E1/D, E2/D, nu12*E2/D, G12

def transform_Q(Q11, Q22, Q12, Q66, angle_deg):
    a = math.radians(angle_deg)
    m, n = math.cos(a), math.sin(a)
    m2, n2, mn = m*m, n*n, m*n
    Qbar11 = Q11*m2**2 + 2*(Q12+2*Q66)*n2*m2 + Q22*n2**2
    Qbar22 = Q11*n2**2 + 2*(Q12+2*Q66)*n2*m2 + Q22*m2**2
    Qbar12 = (Q11+Q22-4*Q66)*n2*m2 + Q12*(m2**2+n2**2)
    Qbar66 = (Q11+Q22-2*Q12-2*Q66)*n2*m2 + Q66*(n2**2+m2**2)
    Qbar16 = (Q11-Q12-2*Q66)*m2*m*n - (Q22-Q12-2*Q66)*n2*m*n
    Qbar26 = (Q11-Q12-2*Q66)*n2*m*n - (Q22-Q12-2*Q66)*m2*m*n
    return Qbar11, Qbar22, Qbar12, Qbar66, Qbar16, Qbar26

def build_ABD(ply_stack):
    """
    ply_stack: list of (angle_deg, E1, E2, G12, nu12, t)
    Returns full 6x6 ABD matrix components needed for flutter.
    """
    n_plies = len(ply_stack)
    t_total = sum(p[5] for p in ply_stack)
    z_k = -t_total / 2.0

    A = [[0.0]*6 for _ in range(6)]  # We'll track 11,22,12,66,16,26
    D = [[0.0]*6 for _ in range(6)]

    A11=A22=A12=A66=A16=A26 = 0.0
    D11=D22=D12=D66=D16=D26 = 0.0

    for (angle, e1, e2, g12, nu12, t) in ply_stack:
        q11,q22,q12,q66 = Q_matrix(e1, e2, g12, nu12)
        qb11,qb22,qb12,qb66,qb16,qb26 = transform_Q(q11,q22,q12,q66,angle)
        z0, z1 = z_k, z_k + t
        dz  = z1 - z0
        dz3 = (z1**3 - z0**3) / 3.0
        A11 += qb11*dz;  A22 += qb22*dz;  A12 += qb12*dz
        A66 += qb66*dz;  A16 += qb16*dz;  A26 += qb26*dz
        D11 += qb11*dz3; D22 += qb22*dz3; D12 += qb12*dz3
        D66 += qb66*dz3; D16 += qb16*dz3; D26 += qb26*dz3
        z_k = z1

    G_xy = A66 / t_total   # In-plane shear modulus [Pa]
    G_eff = 12.0 * D66 / t_total**3  # Bending-equivalent shear [Pa]

    return {
        'A11': A11, 'A22': A22, 'A12': A12, 'A66': A66,
        'D11': D11, 'D22': D22, 'D12': D12, 'D66': D66,
        'D16': D16, 'D26': D26,
        'G_xy': G_xy, 'G_eff': G_eff,
        't_total': t_total
    }

# =============================================================================
# LAYUP CANDIDATES
# Target total thickness: 6 mm (includes tip-to-tip, so fin skin ~3 mm each side)
# Available plies:
#   DB300 pair (+45/-45): t=0.333mm  → use as bonded pair (always ±45)
#   GA90R woven (0/90)  : t=0.336mm  → single layer gives both 0° and 90°
#
# Candidate symmetric half-layups (from midplane outward):
#   Build around a CORE of 0/90 woven + ±45 skins strategy.
#   Total target: ~6mm → each skin ~3mm → ~9 plies per side
#   
# Strategy A: Maximize D66 → ±45 dominant
# Strategy B: Maximize D11 → 0/90 dominant  
# Strategy C: Balanced (typical aerospace practice)
# =============================================================================

def make_ply(angle, fabric='GA90R'):
    if fabric == 'GA90R':
        return (angle, E1_GA90R, E2_GA90R, G12_GA90R, nu12_GA90R, t_GA90R_ply)
    else:  # DB300 individual ±45 ply
        return (angle, E1_DB300, E2_DB300, G12_DB300, nu12_DB300, t_DB300_single)

def mirror(half):
    """Create symmetric laminate from half-stack (midplane symmetry)."""
    return list(half) + list(reversed(half))

# Build candidate layups — each is a FULL symmetric laminate
# Notation: W=woven(0/90), B+=+45 NCF, B-=-45 NCF
# Woven ply at 0° represents the 0/90 fabric (both directions active)

candidates = {}

# --- Candidate A: Flutter-optimized (±45 skins, 0/90 core) ---
# Half: [B+, B-, W, W, B+, B-, W, W, B+]  ~9 plies × 0.333mm ≈ 3.0 mm
half_A = [
    make_ply(+45, 'DB300'), make_ply(-45, 'DB300'),
    make_ply(  0, 'GA90R'), make_ply(  0, 'GA90R'),
    make_ply(+45, 'DB300'), make_ply(-45, 'DB300'),
    make_ply(  0, 'GA90R'),
    make_ply(+45, 'DB300'), make_ply(-45, 'DB300'),
]
candidates['A_flutter_opt'] = mirror(half_A)

# --- Candidate B: Bending-optimized (0/90 dominant, ±45 inter-leaved) ---
half_B = [
    make_ply(  0, 'GA90R'), make_ply(  0, 'GA90R'),
    make_ply(+45, 'DB300'), make_ply(-45, 'DB300'),
    make_ply(  0, 'GA90R'), make_ply(  0, 'GA90R'),
    make_ply(+45, 'DB300'), make_ply(-45, 'DB300'),
    make_ply(  0, 'GA90R'),
]
candidates['B_bending_opt'] = mirror(half_B)

# --- Candidate C: Quasi-isotropic (aerospace standard baseline) ---
half_C = [
    make_ply(  0, 'GA90R'),
    make_ply(+45, 'DB300'), make_ply(-45, 'DB300'),
    make_ply(  0, 'GA90R'),
    make_ply(+45, 'DB300'), make_ply(-45, 'DB300'),
    make_ply(  0, 'GA90R'),
    make_ply(+45, 'DB300'), make_ply(-45, 'DB300'),
]
candidates['C_quasi_isotropic'] = mirror(half_C)

# --- Candidate D: All ±45 (maximum torsion, minimum bending — reference) ---
half_D = [
    make_ply(+45, 'DB300'), make_ply(-45, 'DB300'),
    make_ply(+45, 'DB300'), make_ply(-45, 'DB300'),
    make_ply(+45, 'DB300'), make_ply(-45, 'DB300'),
    make_ply(+45, 'DB300'), make_ply(-45, 'DB300'),
    make_ply(+45, 'DB300'), make_ply(-45, 'DB300'),
]
candidates['D_all_45'] = mirror(half_D)

# --- Candidate E: All 0/90 (maximum bending — reference) ---
half_E = [
    make_ply(0, 'GA90R'), make_ply(0, 'GA90R'), make_ply(0, 'GA90R'),
    make_ply(0, 'GA90R'), make_ply(0, 'GA90R'), make_ply(0, 'GA90R'),
    make_ply(0, 'GA90R'), make_ply(0, 'GA90R'), make_ply(0, 'GA90R'),
]
candidates['E_all_0_90'] = mirror(half_E)

# =============================================================================
# COMPUTE AND RANK
# =============================================================================
print("\n" + "=" * 65)
print("  LAMINATE STIFFNESS RESULTS — All Candidates")
print("=" * 65)
print(f"  {'Layup':<22} {'t[mm]':>6} {'D66[N·m]':>10} {'D11[N·m]':>10} "
      f"{'G_xy[GPa]':>10} {'D16[N·m]':>10}")
print("  " + "-"*62)

results = {}
for name, stack in candidates.items():
    r = build_ABD(stack)
    results[name] = r
    print(f"  {name:<22} {r['t_total']*1e3:>6.2f} {r['D66']:>10.2f} "
          f"{r['D11']:>10.2f} {r['G_xy']/1e9:>10.3f} {r['D16']:>10.3f}")

# Best for flutter = highest D66
best_flutter = max(results, key=lambda k: results[k]['D66'])
best_bending = max(results, key=lambda k: results[k]['D11'])
rec = results[best_flutter]

print(f"\n  ★  Best for flutter resistance : {best_flutter}  (D66={rec['D66']:.2f} N·m)")
print(f"  ★  Best for bending stiffness  : {best_bending}  (D11={results[best_bending]['D11']:.2f} N·m)")

# =============================================================================
# EXPORT KEY VALUES FOR FLUTTER SCRIPT
# =============================================================================
print("\n" + "=" * 65)
print("  RECOMMENDED LAYUP — Flutter-Optimized (for flutterEstimate)")
print("=" * 65)
print(f"  Layup       : {best_flutter}")
print(f"  t_total     : {rec['t_total']*1e3:.3f} mm  (target: 6.0 mm)")
print(f"  D66         : {rec['D66']:.3f} N·m   ← feed to flutterEstimate.py")
print(f"  D11         : {rec['D11']:.3f} N·m")
print(f"  G_xy (A66)  : {rec['G_xy']/1e9:.3f} GPa")
print(f"  G_eff(D66)  : {rec['G_eff']/1e9:.3f} GPa")
print(f"  D16 (bend-twist coupling) : {rec['D16']:.4f} N·m")
note = "(near-zero D16 → no bend-twist coupling — preferred)" if abs(rec['D16']) < 1.0 else "(non-zero D16 — bend-twist coupling present)"
print(f"  {note}")

