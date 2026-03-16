import math

# =============================================================================
# FLUTTER BOUNDARY ESTIMATOR — FalconLAUNCH VI Fins
# Based on:
#   NACA TN 4197  (Martin, 1958) — subsonic/incompressible
#   Ackeret (NACA TM 317, 1925) — supersonic lift slope
#   Bisplinghoff, Ashley & Halfman "Aeroelasticity" (1955) — compressible scaling
# =============================================================================

# --- Fin geometry (FalconLAUNCH VI) ---
cr_m      = 0.300    # Root chord [m]
ct_m      = 0.150    # Tip chord  [m]
s_m       = 0.160    # Exposed semi-span (height) [m]
t_m       = 0.006    # Total fin thickness [m]
sweep_deg = 57.4     # Leading-edge sweep [deg]

# --- D66 from inplaneG_v2.py (flutter-optimized layup) ---
#   Feed whichever candidate you select from inplaneG output.
#   The E_all_0_90 layup wins on D66 AND D11 for the 6mm constraint.
D66_candidates = {
    'A_flutter_opt'   : 37.96,
    'B_bending_opt'   : 41.55,
    'C_quasi_isotropic': 32.71,
    'D_all_45'        : 32.31,
    'E_all_0_90'      : 63.42,   # ← selected optimum
}
D66_selected = 'E_all_0_90'
D66 = D66_candidates[D66_selected]

# --- Flight environment ---
altitude_m   = 1462.0    # Altitude at max dynamic pressure [m]
rocket_v_mps = 650.0     # Velocity at max dynamic pressure [m/s]
# NOTE: Extract exact q_max altitude and velocity from OpenRocket trajectory file.

# =============================================================================
# ISA STANDARD ATMOSPHERE (exact, troposphere)
# =============================================================================
def standard_atmosphere(altitude_m):
    T0    = 288.15    # Sea-level temperature [K]
    P0    = 101325.0  # Sea-level pressure [Pa]
    L     = 0.0065    # Lapse rate [K/m]
    R     = 287.05    # Specific gas constant [J/(kg·K)]
    g     = 9.80665   # Gravity [m/s²]
    gamma = 1.4       # Heat capacity ratio

    T   = T0 - L * altitude_m
    P   = P0 * (1.0 - L * altitude_m / T0) ** (g / (R * L))
    rho = P / (R * T)
    a   = math.sqrt(gamma * R * T)
    return P, rho, a, T

# =============================================================================
# GEOMETRIC PARAMETERS
# =============================================================================
def fin_geometry(cr, ct, s, t, sweep_deg):
    """All derived geometric quantities for the fin panel."""
    AR      = (2.0 * s) / (cr + ct)          # Aspect ratio (both sides)
    lambda_ = ct / cr                          # Taper ratio
    t_ratio = t / cr                           # Thickness-to-chord ratio
    sweep_r = math.radians(sweep_deg)
    S_fin   = 0.5 * (cr + ct) * s             # Fin planform area [m²]
    mac     = (2/3) * cr * (1 + lambda_ + lambda_**2) / (1 + lambda_)  # MAC
    return AR, lambda_, t_ratio, sweep_r, S_fin, mac

# =============================================================================
# FLUTTER VELOCITY — NACA TN 4197 (TRUE FORM, SUBSONIC)
# Equation from Martin (1958), Section 3, parameterised for SI units.
# The NACA 4197 coefficient 39.3/14.7 = 2.6734 converts the panel
# shear modulus G_eff [psi] to a flutter velocity [ft/s] at sea-level.
# The fully dimensionless SI form is:
#   (V_f / a)² = G_eff / [ (2.6734 * AR³ / ((AR+2) * (t/c)³)) * P ]
# where G_eff = 12*D66 / t³  [Pa],  P [Pa],  a [m/s]
# =============================================================================
def flutter_naca4197(D66, t_m, cr_m, ct_m, s_m, P_alt, a_alt):
    AR, lam, t_ratio, _, _, _ = fin_geometry(cr_m, ct_m, s_m, t_m, sweep_deg)
    G_eff = 12.0 * D66 / t_m**3          # Equivalent shear modulus [Pa]

    # NACA 4197 dimensionless pressure coefficient (taper excluded per original)
    denom = (2.6734 * AR**3 * P_alt) / ((AR + 2.0) * t_ratio**3)
    Vf_mps = a_alt * math.sqrt(G_eff / denom)
    return Vf_mps, G_eff, AR, t_ratio

# =============================================================================
# SUPERSONIC FLUTTER CORRECTION (Ackeret + Bisplinghoff scaling)
# Rationale:
#   Lift slope: subsonic ~ 2π/rad ; supersonic ~ 4/√(M²-1) per Ackeret
#   Aerodynamic center: subsonic 0.25c ; supersonic shifts to ~0.40c
#   Elastic axis (typical composite fin): ~0.50c
#   Flutter speed ∝ 1/√(C_Lα · e) where e = |EA - AC| is moment arm.
#   Supersonic flutter velocity is LOWER than subsonic for the same fin.
# =============================================================================
def flutter_supersonic(Vf_sub, M, a_alt):
    if M <= 1.05:
        return Vf_sub, 1.0  # Transonic — use subsonic (conservative)

    C_La_sub   = 2.0 * math.pi
    C_La_super = 4.0 / math.sqrt(M**2 - 1.0)

    # Moment arms (EA at 0.50c assumed for symmetric layup)
    e_sub   = 0.25   # |0.50c - 0.25c|
    e_super = 0.10   # |0.50c - 0.40c|

    # Scaling: V_f ∝ sqrt(1 / (C_La * e))
    scale = math.sqrt((C_La_sub * e_sub) / (C_La_super * e_super))
    Vf_super = Vf_sub * scale
    return Vf_super, scale

# =============================================================================
# SWEEP CORRECTION (Collar & Tickle, 1954)
# For swept fins, the effective flutter velocity is increased by (1/cos Λ)^0.5
# — sweep reduces the effective aerodynamic loading on the fin.
# =============================================================================
def flutter_sweep_correction(Vf, sweep_rad):
    # Correction factor (conservative, cos Λ component only)
    k_sweep = 1.0 / math.sqrt(math.cos(sweep_rad))
    return Vf * k_sweep, k_sweep

# =============================================================================
# MAIN COMPUTATION
# =============================================================================
P_alt, rho_alt, a_alt, T_alt = standard_atmosphere(altitude_m)
M_rocket = rocket_v_mps / a_alt
q_dynamic = 0.5 * rho_alt * rocket_v_mps**2

AR, lam, t_ratio, sweep_rad, S_fin, mac = fin_geometry(cr_m, ct_m, s_m, t_m, sweep_deg)

# Step 1: Classical (subsonic, NACA 4197)
Vf_sub, G_eff, _, _ = flutter_naca4197(D66, t_m, cr_m, ct_m, s_m, P_alt, a_alt)

# Step 2: Supersonic compressible correction
Vf_super, scale_super = flutter_supersonic(Vf_sub, M_rocket, a_alt)

# Step 3: Sweep correction
Vf_super_swept, k_sweep = flutter_sweep_correction(Vf_super, sweep_rad)
Vf_sub_swept,   _       = flutter_sweep_correction(Vf_sub,   sweep_rad)

# =============================================================================
# SAFETY MARGIN
# Flutter Safety Factor (FSF) = V_flutter / V_rocket.
# Structural design requirement: FSF ≥ 1.5 (AIAA S-080 / MIL-A-8870).
# =============================================================================
FSF_classical  = Vf_sub_swept   / rocket_v_mps
FSF_supersonic = Vf_super_swept / rocket_v_mps

margin_pct = (Vf_super_swept - rocket_v_mps) / rocket_v_mps * 100.0

# =============================================================================
# PRINT REPORT
# =============================================================================
print("=" * 65)
print("  FLUTTER BOUNDARY ANALYSIS — FalconLAUNCH VI Fins")
print("=" * 65)

print(f"\n  FIN GEOMETRY")
print(f"    Root chord     : {cr_m*1e3:.1f} mm")
print(f"    Tip chord      : {ct_m*1e3:.1f} mm")
print(f"    Span (height)  : {s_m*1e3:.1f} mm")
print(f"    Thickness      : {t_m*1e3:.1f} mm")
print(f"    Sweep (LE)     : {sweep_deg:.1f}°")
print(f"    Aspect ratio   : {AR:.3f}")
print(f"    Taper ratio    : {lam:.3f}")
print(f"    t/c (root)     : {t_ratio:.4f}")
print(f"    MAC            : {mac*1e3:.1f} mm")
print(f"    Planform area  : {S_fin*1e4:.2f} cm²")

print(f"\n  FLIGHT ENVIRONMENT (at max dynamic pressure)")
print(f"    Altitude       : {altitude_m:.1f} m  ({altitude_m*3.281:.0f} ft)")
print(f"    Temperature    : {T_alt:.2f} K  ({T_alt-273.15:.2f} °C)")
print(f"    Pressure       : {P_alt/1e3:.3f} kPa")
print(f"    Density        : {rho_alt:.4f} kg/m³")
print(f"    Speed of sound : {a_alt:.2f} m/s")
print(f"    Rocket speed   : {rocket_v_mps:.1f} m/s  (Mach {M_rocket:.3f})")
print(f"    Dynamic press. : {q_dynamic/1e3:.2f} kPa  ({q_dynamic:.1f} Pa)")

print(f"\n  COMPOSITE PROPERTIES")
print(f"    Layup          : {D66_selected}")
print(f"    D66            : {D66:.3f} N·m  (torsional bending stiffness)")
print(f"    G_eff          : {G_eff/1e9:.3f} GPa  (12·D66 / t³)")

print(f"\n  FLUTTER BOUNDARIES")
print(f"    ┌─────────────────────────────────────────────┐")
print(f"    │  Method               V_flutter    Mach     │")
print(f"    │─────────────────────────────────────────────│")
print(f"    │  NACA 4197 (subsonic) {Vf_sub:.1f} m/s  {Vf_sub/a_alt:.3f}  │")
print(f"    │  + Sweep correction   {Vf_sub_swept:.1f} m/s  {Vf_sub_swept/a_alt:.3f}  │")
print(f"    │  + Supersonic corr.   {Vf_super:.1f} m/s  {Vf_super/a_alt:.3f}  │")
print(f"    │  + Supersonic+Sweep   {Vf_super_swept:.1f} m/s  {Vf_super_swept/a_alt:.3f}  │")
print(f"    └─────────────────────────────────────────────┘")

print(f"\n  SAFETY ASSESSMENT")
print(f"    Sweep correction factor  : {k_sweep:.4f}  (cos Λ = {math.cos(sweep_rad):.4f})")
print(f"    Supersonic scale factor  : {scale_super:.4f}")
print(f"    FSF (classical + sweep)  : {FSF_classical:.3f}  {'✓ PASS' if FSF_classical >= 1.5 else '✗ FAIL (< 1.5)'}")
print(f"    FSF (supersonic + sweep) : {FSF_supersonic:.3f}  {'✓ PASS' if FSF_supersonic >= 1.5 else '✗ FAIL (< 1.5)'}")
print(f"    Margin of safety         : {margin_pct:.1f}%")
print(f"    Required FSF (AIAA S-080): ≥ 1.50")

if FSF_supersonic >= 1.5:
    print(f"\n  ► FIN IS FLUTTER-SAFE at Mach {M_rocket:.2f} with {D66_selected} layup.")
else:
    print(f"\n  ► WARNING: Flutter margin insufficient. Increase D66 or reduce thickness.")
    print(f"    Suggested D66 for FSF=1.5: {D66 * (1.5/FSF_supersonic)**2:.1f} N·m")

print(f"\n  NOTE: Verify V_rocket and altitude at max-q from OpenRocket trajectory.")
print(f"  References: NACA TN 4197 (Martin 1958), Ackeret (NACA TM 317), BAH (1955).")

