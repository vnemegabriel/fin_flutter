"""
flutterEstimate_v3.py  —  FalconLAUNCH VI Fin Flutter Boundary Analysis
========================================================================
Theory  : NACA TN 4197 (Martin 1958), Ackeret (NACA TM 317, 1925),
          sweep correction from Bisplinghoff/Ashley/Halfman (1955)
Pipeline: inplaneG_v5 JSON → ISA atmosphere → subsonic → supersonic
          → sweep-corrected flutter speed → FSF vs max-q condition

Usage
-----
  python3 flutterEstimate_v3.py                          # defaults
  python3 flutterEstimate_v3.py --json laminate.json     # from inplaneG pipeline
  python3 flutterEstimate_v3.py --d66 205.96             # manual D66 [N.m]
  python3 flutterEstimate_v3.py --altitude 1462 --mach 1.942
  python3 flutterEstimate_v3.py --sweep-table            # full altitude sweep table

References
----------
[1] Martin, H.C. NACA TN 4197 (1958)
[2] Ackeret, J.   NACA TM 317  (1925)
[3] Bisplinghoff, Ashley & Halfman. Aeroelasticity (1955) Sec 5.5
[4] AIAA S-080 (1999); MIL-A-8870C (1993)
"""
import math, json, argparse, sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent / "core" / "python"))
from fin_flutter.flight_data import read_flight_data

# ─── CLI ────────────────────────────────────────────────────────────────────
def parse_args():
    p = argparse.ArgumentParser(description="FalconLAUNCH VI Flutter v3")
    p.add_argument("--json",        type=str,   default=None)
    p.add_argument("--d66",         type=float, default=None,  help="D66 [N.m]")
    p.add_argument("--t",           type=float, default=None,  help="thickness [mm]")
    p.add_argument("--altitude",    type=float, default=1462.0,help="max-q altitude [m]")
    p.add_argument("--velocity",    type=float, default=None,  help="rocket speed [m/s]")
    p.add_argument("--mach",        type=float, default=1.942, help="rocket Mach")
    p.add_argument("--span",        type=float, default=0.160, help="semi-span [m]")
    p.add_argument("--cr",          type=float, default=0.300, help="root chord [m]")
    p.add_argument("--ct",          type=float, default=0.150, help="tip chord [m]")
    p.add_argument("--sweep",       type=float, default=57.4,  help="LE sweep [deg]")
    p.add_argument("--sweep-off",   action="store_true")
    p.add_argument("--super-off",   action="store_true")
    p.add_argument("--fsf-req",     type=float, default=1.50)
    p.add_argument("--sweep-table", action="store_true")
    p.add_argument("--csv",         type=str,   default=None,
                   help="flight data CSV; sets --velocity and --altitude from max-velocity row")
    return p.parse_args()

# ─── ISA ATMOSPHERE ─────────────────────────────────────────────────────────
def isa(h):
    """ICAO ISA troposphere (0-11 km). Returns rho, p, T, a."""
    assert 0 <= h <= 11000, f"ISA troposphere only valid 0–11 km, got {h} m"  # ERRORS.md F-5
    T0, p0, L, g, R, gam = 288.15, 101325.0, 0.0065, 9.80665, 287.058, 1.4
    T   = T0 - L * h
    p   = p0 * (T / T0) ** (g / (R * L))
    rho = p / (R * T)
    a   = math.sqrt(gam * R * T)
    return rho, p, T, a

# ─── FIN GEOMETRY ───────────────────────────────────────────────────────────
def fin_ar(b, cr, ct):
    return b**2 / (0.5*(cr+ct)*b)

def fin_mac(b, cr, ct):
    lam = ct/cr
    return cr * (2/3) * (1 + lam + lam**2) / (1 + lam)

# ─── FLUTTER CHAIN ──────────────────────────────────────────────────────────
def flutter_naca4197(G_eff, rho, t, b, AR, tc):
    """
    NACA TN 4197 subsonic flutter speed.

    Physical model
    --------------
    Torsional restoring work = aerodynamic forcing work at critical speed.
    For a uniform trapezoidal panel the result is:

        V_f = sqrt( G_eff * t^3 / (rho * b^2 * K) ) * (tc_ref/tc)^0.5

    K = aerodynamic influence factor (NACA 4197 Table 2, interpolated).
    tc correction accounts for aerodynamic coupling scaling with t/c.

    Parameters
    ----------
    G_eff : float  12*D66/t^3 [Pa]  — bending-equiv. torsional modulus
    rho   : float  air density [kg/m3]
    t     : float  fin thickness [m]
    b     : float  semi-span [m]
    AR    : float  aspect ratio
    tc    : float  thickness/chord ratio

    Note
    ----
    D16 bend-twist coupling is not included. For β=20° (D16<0, washout), the
    true flutter speed is higher than computed — conservative. For β>45°
    (D16>0, washin), this approach is non-conservative. See ERRORS.md F-6.
    Reference: Weisshaar (1981) J. Aircraft 18(8):669-676
    """
    # K: single-variable linear approximation; taper ratio lambda is NOT included.
    # WARNING (ERRORS.md F-4): unvalidated against NACA TN 4197 Table 2 at
    # design point AR=0.711, lambda=0.5. Validation requires 2D (AR, lambda)
    # bilinear interpolation from Table 2 — [Martin, NACA TN 4197, 1958 Table 2].
    K  = 0.65 * (1.0 + 0.10*(AR - 1.0))
    K  = max(K, 0.40)
    V2 = G_eff * t**3 / (rho * b**2 * K)
    Vf = math.sqrt(max(V2, 0.0))
    tc_ref = 0.02
    if tc > 0:
        Vf *= (tc_ref / tc)**0.5
    return Vf

def ackeret_factor(M, x_ea=0.40):
    """
    Supersonic flutter speed amplification factor (Ackeret aerodynamics).

    Ratio of subsonic-to-supersonic aerodynamic torsional coupling, computed
    as moment arms about the fin elastic axis x_ea (fraction of chord):

        sub:   CL_alpha = 2*pi,             AC at 0.25c (subsonic thin airfoil)
        super: CL_alpha = 4/sqrt(M^2-1),    AC at 0.50c (Ackeret symmetric airfoil)

        e_sub = x_ea - 0.25   moment arm, sub  AC to EA
        e_sup = 0.50 - x_ea   moment arm, super AC to EA

        Factor = sqrt( (CL_sub * e_sub) / (CL_sup * e_sup) )

    Parameters
    ----------
    M    : float  Mach number at the flutter condition (M > 1 for supersonic branch)
    x_ea : float  elastic axis location as fraction of chord (default 0.40c)

    References
    ----------
    # Eq. 5.5 — moment arm formulation [Bisplinghoff, Ashley & Halfman, 1955 §5.5]
    # Ackeret AC at 0.50c for symmetric thin airfoil — [Ackeret, NACA TM 317, 1925]
    """
    if M <= 1.0:
        return 1.0
    e_sub = x_ea - 0.25          # moment arm: subsonic  AC (0.25c) to EA — [BAH 1955 §5.5]
    e_sup = 0.50 - x_ea          # moment arm: supersonic AC (0.50c) to EA — [Ackeret, NACA TM 317]
    cl_sub  = 2.0 * math.pi      # subsonic lift slope [/rad] — thin-airfoil theory
    cl_sup  = 4.0 / math.sqrt(M**2 - 1.0)  # Ackeret lift slope — [NACA TM 317 Eq. 16]
    return math.sqrt((cl_sub * e_sub) / (cl_sup * e_sup))

def sweep_factor(sweep_deg):
    """
    Sweep-corrected flutter speed factor: V_swept = V_unswept / cos(Lambda).

    Dynamic pressure at flutter scales as cos²(Λ), so flutter speed scales
    as cos(Λ). The correction factor is therefore 1/cos(Λ).

    # Eq. 5.5 — [Bisplinghoff, Ashley & Halfman, 1955 §5.5]
    # q_flutter ∝ cos²(Λ)  =>  Vf ∝ cos(Λ)  =>  factor = 1/cos(Λ)
    """
    if sweep_deg <= 0:
        return 1.0
    return 1.0 / math.cos(math.radians(sweep_deg))  # Eq. 5.5 — [BAH 1955 §5.5]

def compute_flutter(D66, t_mm, b, cr, ct, sweep_deg,
                    h, V_rocket=None, M_rocket=None,
                    do_sweep=True, do_super=True):
    t = t_mm * 1e-3
    rho, p, T, a = isa(h)
    if V_rocket is None:
        V_rocket = M_rocket * a
    else:
        M_rocket = V_rocket / a

    q    = 0.5 * rho * V_rocket**2
    AR   = fin_ar(b, cr, ct)
    mac  = fin_mac(b, cr, ct)
    tc   = t / mac
    Geff = 12.0 * D66 / t**3

    # Step 1 — NACA 4197 subsonic
    Vf1   = flutter_naca4197(Geff, rho, t, b, AR, tc)
    Mf1   = Vf1 / a

    # Step 2 — Ackeret supersonic correction at self-consistent flutter Mach
    # ERROR F-2 fix: evaluate ackeret_factor at the flutter Mach Mf, not at
    # M_rocket (flight Mach). Solve Mf = Mf1 * ackeret_factor(Mf) by fixed-
    # point iteration until |Mf_new - Mf| < 1e-6.
    # Reference: [Bisplinghoff, Ashley & Halfman, 1955 §5.5]
    if do_super:
        Mf_iter = Mf1
        for _ in range(50):
            f_sup_iter = ackeret_factor(Mf_iter)
            Mf_new = Mf1 * f_sup_iter
            if abs(Mf_new - Mf_iter) < 1e-6:
                break
            Mf_iter = Mf_new
        f_sup = f_sup_iter
    else:
        f_sup = 1.0
    Vf2   = Vf1 * f_sup
    Mf2   = Vf2 / a
    f_sw  = sweep_factor(sweep_deg) if do_sweep else 1.0
    Vf3   = Vf2 * f_sw
    Mf3   = Vf3 / a
    FSF   = Vf3 / V_rocket

    return dict(
        h=h, rho=rho, T=T, p=p, a=a,
        V_rocket=V_rocket, M_rocket=M_rocket, q=q,
        AR=AR, mac=mac, tc=tc, Geff=Geff,
        Vf1=Vf1, Mf1=Mf1, f_sup=f_sup,
        Vf2=Vf2, Mf2=Mf2, f_sw=f_sw,
        Vf3=Vf3, Mf3=Mf3, FSF=FSF,
    )

# ─── REPORT ─────────────────────────────────────────────────────────────────
L72 = "=" * 72

def report(r, D66, t_mm, cr, ct, b, sweep, fsf_req,
           do_sweep, do_super, src):
    print(L72)
    print("  FIN FLUTTER ANALYSIS  —  FalconLAUNCH VI  [flutterEstimate_v3.py]")
    print(L72)
    if src: print(f"  Laminate file : {src}  (tailored D66, β=20°, conservative)")
    print(f"  D66  = {D66:.4f} N.m     G_eff = {r['Geff']/1e9:.4f} GPa")
    print(f"  t    = {t_mm:.4f} mm")
    print()
    print(f"  FIN GEOMETRY")
    print(f"    cr = {cr*1e3:.1f} mm   ct = {ct*1e3:.1f} mm   b = {b*1e3:.1f} mm")
    print(f"    Λ  = {sweep:.1f} deg   AR = {r['AR']:.4f}   MAC = {r['mac']*1e3:.2f} mm")
    print(f"    t/c = {r['tc']:.4f}  ({r['tc']*100:.2f}%)")
    print()
    print(f"  ISA ATMOSPHERE  h = {r['h']:.1f} m")
    print(f"    ρ = {r['rho']:.4f} kg/m³   T = {r['T']:.2f} K   a = {r['a']:.2f} m/s")
    print(f"    p = {r['p']/1e3:.3f} kPa    q = {r['q']/1e3:.3f} kPa")
    print()
    print(f"  ROCKET  Vr = {r['V_rocket']:.1f} m/s   M = {r['M_rocket']:.4f}")
    print()
    print(f"  FLUTTER BOUNDARY CHAIN")
    print(f"    [1] NACA TN 4197 (subsonic)     : Vf = {r['Vf1']:.1f} m/s   Mf = {r['Mf1']:.3f}")
    if do_super:
        print(f"    [2] Ackeret supersonic (M={r['M_rocket']:.3f}) : "
              f"×{r['f_sup']:.4f}  → Vf = {r['Vf2']:.1f} m/s   Mf = {r['Mf2']:.3f}")
    else:
        print(f"    [2] Supersonic correction : DISABLED")
    if do_sweep:
        print(f"    [3] Sweep correction (Λ={sweep:.1f}°)  : "
              f"×{r['f_sw']:.4f}  → Vf = {r['Vf3']:.1f} m/s   Mf = {r['Mf3']:.3f}")
    else:
        print(f"    [3] Sweep correction : DISABLED")
    print()
    margin = (r['FSF']/fsf_req - 1.0)*100
    status = "PASS ✓" if r['FSF'] >= fsf_req else "FAIL ✗"
    print(f"  FLUTTER SAFETY FACTOR")
    print(f"    V_flutter = {r['Vf3']:.1f} m/s   V_rocket = {r['V_rocket']:.1f} m/s")
    print(f"    FSF = {r['FSF']:.4f}   [{status}]   required ≥ {fsf_req:.2f}   margin = {margin:+.1f}%")
    print()
    if r['FSF'] >= fsf_req:
        print(f"  ✓  Flutter boundary ({r['Vf3']:.0f} m/s) clears flight speed ({r['V_rocket']:.0f} m/s)")
        print(f"     by {margin:.1f}%.  Laminate is aeroelastically safe at max-q.")
    else:
        print(f"  ✗  Insufficient flutter margin. Increase D66 (more ±45 plies)")
        print(f"     or reduce fin span. Target: V_flutter ≥ {r['V_rocket']*fsf_req:.0f} m/s")

def altitude_table(D66, t_mm, b, cr, ct, sweep, M0,
                   do_sweep, do_super, fsf_req):
    print(f"\n{L72}")
    print(f"  ALTITUDE SWEEP  (M_rocket = {M0:.3f} constant)")
    print(L72)
    print(f"  {'h[m]':>6}  {'rho':>8}  {'a[m/s]':>7}  {'Vr[m/s]':>8}  "
          f"{'Vf[m/s]':>8}  {'Mf':>6}  {'FSF':>6}  Status")
    print(f"  {'-'*70}")
    for h in range(0, 4001, 250):
        rho, _, _, a = isa(h)
        Vr = M0 * a
        r = compute_flutter(D66, t_mm, b, cr, ct, sweep, h,
                            V_rocket=Vr, do_sweep=do_sweep, do_super=do_super)
        ok  = "PASS" if r['FSF'] >= fsf_req else "FAIL"
        tag = " ← max-q" if abs(h-1462) < 130 else ""
        print(f"  {h:>6d}  {rho:>8.4f}  {a:>7.2f}  {Vr:>8.1f}  "
              f"{r['Vf3']:>8.1f}  {r['Mf3']:>6.3f}  {r['FSF']:>6.3f}  {ok}{tag}")

# ─── MAIN ───────────────────────────────────────────────────────────────────
def main():
    args = parse_args()
    D66, t_mm, src = None, None, None

    if args.json:
        with open(args.json) as f:
            lam = json.load(f)
        D66  = lam['tailored_beta']['D66_Nm']
        t_mm = lam['tailored_beta']['t_total_mm']
        src  = args.json

    if args.d66 is not None: D66  = args.d66
    if args.t   is not None: t_mm = args.t
    if D66  is None: D66  = 205.96;  print(f"  Default D66 = {D66} N.m")
    if t_mm is None: t_mm = 5.356;   print(f"  Default t   = {t_mm} mm")

    if args.csv:
        max_vel, alt_ft = read_flight_data(args.csv)
        # CSV altitude is in ft; ISA expects metres
        if args.velocity is None:
            args.velocity = max_vel
        if args.altitude == 1462.0:          # still at default — override
            args.altitude = alt_ft * 0.3048  # ft → m

    do_sweep = not args.sweep_off
    do_super = not args.super_off
    sweep    = 0.0 if args.sweep_off else args.sweep

    r = compute_flutter(D66, t_mm, args.span, args.cr, args.ct, sweep,
                        args.altitude,
                        V_rocket=args.velocity, M_rocket=args.mach,
                        do_sweep=do_sweep, do_super=do_super)
    print()
    report(r, D66, t_mm, args.cr, args.ct, args.span, sweep,
           args.fsf_req, do_sweep, do_super, src)

    if args.sweep_table:
        altitude_table(D66, t_mm, args.span, args.cr, args.ct, sweep,
                       r['M_rocket'], do_sweep, do_super, args.fsf_req)

    print(f"\n{L72}")

if __name__ == "__main__":
    main()
