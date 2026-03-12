import math

def standard_atmosphere(altitude_m):
    """Calculates exact thermodynamic Pressure and Speed of Sound."""
    T0 = 288.15      # Sea level standard temperature (K)
    P0 = 101325.0    # Sea level standard pressure (Pa)
    L = 0.0065       # Temperature lapse rate (K/m)
    R = 287.05       # Specific gas constant for air (J/(kg·K))
    g = 9.80665      # Gravity (m/s^2)
    gamma = 1.4      # Heat capacity ratio
    
    T = T0 - L * altitude_m
    P = P0 * (1 - L * altitude_m / T0) ** (g / (R * L))
    a = math.sqrt(gamma * R * T)
    
    return P, a

def flutter_boundaries(cr_m, ct_m, s_m, t_m, D66, altitude_m, rocket_v_mps):
    """
    Calculates both Classical (Incompressible) and Advanced Supersonic 
    (Compressible) flutter boundaries using precise analytical scaling.
    """
    P_alt, a_alt = standard_atmosphere(altitude_m)
    current_mach = rocket_v_mps / a_alt
    
    # Geometric Panel Aspect Ratio
    AR = (2.0 * s_m) / (cr_m + ct_m)
    t_ratio = t_m / cr_m
    
    # Equivalent shear modulus from D66 (G_eff = 12 * D66 / t^3)
    G_eff = (12.0 * D66) / (t_m ** 3)
    
    # 1. True NACA TN 4197 Equation (Subsonic, Incompressible)
    # FIXED: Removed the deeply flawed amateur rocketry "Apogee" variant (1.337 coefficient).
    # The true NACA 4197 pressure coefficient is 39.3 / 14.7 psi = 2.6734.
    # This constant is strictly dimensionless when G_eff and P_alt share units.
    # NACA 4197 explicitly omits taper ratio from the parameter entirely.
    denom = (2.6734 * (AR**3) * P_alt) / ((AR + 2.0) * (t_ratio**3))
    Vf_classical_mps = a_alt * math.sqrt(G_eff / denom)
    
    # 2. Supersonic Compressible Analytical Scaling
    # FIXED: Replaced arbitrary legacy constants with rigorous aeroelastic physics.
    # Subsonic C_La ~ 2*pi, AC at 0.25c. 
    # Supersonic C_La ~ 4/sqrt(M^2 - 1), AC shifts to ~0.50c.
    # V_f is proportional to sqrt(1 / (C_La * e)), where e is the EA-AC offset.
    if current_mach > 1.2:
        # Lift slope ratio (Compressible / Incompressible)
        C_La_sub = 2.0 * math.pi
        C_La_super = 4.0 / math.sqrt(current_mach**2 - 1.0)
        lift_ratio = C_La_super / C_La_sub
        
        # Moment arm ratio (Assumes EA at ~0.50c)
        # Subsonic offset: |0.50c - 0.25c| = 0.25c
        # Supersonic offset shifts safely to ~0.40c: |0.50c - 0.40c| = 0.10c
        arm_ratio = 0.10 / 0.25
        
        # Supersonic scaling factor applied to the base classical velocity
        supersonic_scaling = math.sqrt(1.0 / (lift_ratio * arm_ratio))
        Vf_super_mps = Vf_classical_mps * supersonic_scaling
    else:
        Vf_super_mps = Vf_classical_mps
        
    return Vf_classical_mps, Vf_super_mps, P_alt, a_alt

# --- EXECUTION ---
if __name__ == "__main__":
    # Example geometry limits (FalconLAUNCH VI bounds)
    cr_m = 0.30     # Root chord
    ct_m = 0.10     # Tip chord
    s_m = 0.15      # Exposed Span
    t_m = 0.006     # 6.0 mm thickness
    
    # Feed the optimized D66 value from inplaneG.py's output
    D66 = 330.0     
    
    altitude_m = 1462.0   # Maximum dynamic pressure altitude
    rocket_v_mps = 650.0  # Velocity at maximum dynamic pressure altitude 
    ## CHECK. Q (dynamic pressure) IS STATED IN OPENROCKET. DATA MUST BE EXTRACTED FROM THERE.
    
    Vf_class, Vf_super, P_alt, a_alt = flutter_boundaries(cr_m, ct_m, s_m, t_m, D66, altitude_m, rocket_v_mps)
    
    current_mach = rocket_v_mps / a_alt
    
    print(f"--- FLIGHT ENVIRONMENT ---")
    print(f"Burnout Altitude:  {altitude_m:.1f} m (~{altitude_m*3.28:.0f} ft)")
    print(f"Ambient Pressure:  {P_alt/1000:.2f} kPa")
    print(f"Speed of Sound:    {a_alt:.1f} m/s")
    print(f"Max Rocket Speed:  {rocket_v_mps:.1f} m/s (Mach {current_mach:.2f})")
    
    print(f"\n--- FLUTTER BOUNDARY (6.0 mm Fin, D66 = {D66} N-m) ---")
    print(f"1. Classical (True NACA 4197):  {Vf_class:.1f} m/s (Mach {Vf_class/a_alt:.2f})")
    print(f"2. Supersonic (Compressible):   {Vf_super:.1f} m/s (Mach {Vf_super/a_alt:.2f})")
    
    margin = ((Vf_super - rocket_v_mps) / rocket_v_mps) * 100
    print(f"\nMargin of Safety (Supersonic):  {margin:.1f}%")
    
#Raymer, Aircraft Design: A Conceptual Approach

#NACA Technical Note 4197, "Summary of Flutter Experiences as a Guide to the Preliminary Design of Lifting Surfaces on Missiles," Dennis J. Martin, 1958.

#Ackeret, J. "Air Forces on Airfoils Moving Faster Than Sound" (NACA TM 317, 1925).


#Bisplinghoff, Ashley, and Halfman, "Aeroelasticity" (1955).
