import math
import itertools

# --- 1. MICROMECHANICS (WOVEN FABRIC HALPIN-TSAI) ---
def calculate_woven_ply_props(Ef1, Gf12, nuf, Em, num, Vf, t_raw, crimp_factor=0.90):
    """
    Calculates effective ply properties for a BIDIRECTIONAL WOVEN cloth.
    Models the cloth as an orthogonal cross-ply of theoretical UD tows,
    then applies an empirical crimp knockdown factor.
    """
    # 1. Calculate properties of a theoretical Unidirectional (UD) ply at the FULL Vf
    E1_ud = Ef1 * Vf + Em * (1 - Vf)
    nu12_ud = nuf * Vf + num * (1 - Vf)
    
    xi_E = 2.0 
    eta_E = ((Ef1 / Em) - 1.0) / ((Ef1 / Em) + xi_E)
    E2_ud = Em * (1.0 + xi_E * eta_E * Vf) / (1.0 - eta_E * Vf)
    
    # 2. Synthesize the Woven Ply (Cross-Ply Approximation)
    # A balanced woven ply is 50% 0-deg and 50% 90-deg by volume.
    # FIXED: We average the UD properties computed at the total Vf, 
    # preventing the double-discount penalty on stiffness.
    E1_woven = ((E1_ud + E2_ud) / 2.0) * crimp_factor
    E2_woven = E1_woven # Balanced fabric
    
    # Poisson's ratio for balanced woven is typically much lower than UD
    nu12_woven = 0.05 
    
    # Shear Modulus (Halpin-Tsai)
    # Shear is heavily matrix-dominated in 0/90 woven cloth
    Gm = Em / (2.0 * (1.0 + num))
    xi_G = 1.0
    eta_G = ((Gf12 / Gm) - 1.0) / ((Gf12 / Gm) + xi_G)
    G12_woven = Gm * (1.0 + xi_G * eta_G * Vf) / (1.0 - eta_G * Vf)
    
    # Compacted thickness
    Vf_raw = 0.50
    t_ply = t_raw * (Vf_raw / Vf)
    
    return E1_woven, E2_woven, G12_woven, nu12_woven, t_ply

# --- 2. CLPT STIFFNESS MATRIX ---
def calculate_stiffness(layup, E1, E2, G12, nu12, t_ply):
    """Calculates laminate extensional (A) and bending (D) matrices."""
    # Base Q matrix (Orthotropic)
    nu21 = nu12 * E2 / E1
    denom = 1.0 - nu12 * nu21
    Q11 = E1 / denom
    Q22 = E2 / denom
    Q12 = nu12 * E2 / denom
    Q66 = G12
    
    A66 = 0.0
    D11 = 0.0
    D66 = 0.0
    
    t_total = len(layup) * t_ply
    z_bottom = -t_total / 2.0
    
    for i, angle_deg in enumerate(layup):
        angle_rad = math.radians(angle_deg)
        m = math.cos(angle_rad)
        n = math.sin(angle_rad)
        
        # Transform Q matrix to Q_bar
        Qbar11 = Q11*m**4 + 2*(Q12 + 2*Q66)*n**2*m**2 + Q22*n**4
        Qbar66 = (Q11 + Q22 - 2*Q12 - 2*Q66)*n**2*m**2 + Q66*(n**4 + m**4)
        
        # Z coordinates for this specific ply
        z0 = z_bottom + i * t_ply
        z1 = z_bottom + (i + 1) * t_ply
        
        # Summation integration
        A66 += Qbar66 * (z1 - z0)
        term_D = (1.0/3.0) * (z1**3 - z0**3)
        D11 += Qbar11 * term_D
        D66 += Qbar66 * term_D
        
    G_xy = A66 / t_total
    return G_xy, D11, D66

# --- 3. EXECUTION & OPTIMIZATION ---
if __name__ == "__main__":
    # T700 Carbon Fiber properties
    Ef1 = 230e9
    Gf12 = 27e9 
    nuf = 0.2
    
    # BMI Resin properties
    Em = 3.5e9
    num = 0.35
    
    # Process params
    Vf = 0.50
    t_raw = 0.28e-3 
    crimp = 0.90 # Account for weave over/under undulations
    
    # The half-laminate to permute. 
    # Note: For woven cloth, "0" means fibers are at 0 and 90. 
    # "45" means fibers are at +45 and -45.
    half_plies = [0, 0, 45, 45, -45, -45]
    
    E1, E2, G12, nu12, t_ply = calculate_woven_ply_props(
        Ef1, Gf12, nuf, Em, num, Vf, t_raw, crimp_factor=crimp
    )
    
    unique_half_layups = set(itertools.permutations(half_plies))
    
    results = []
    for half in unique_half_layups:
        full_layup = list(half) + list(reversed(half))
        G_xy, D11, D66 = calculate_stiffness(full_layup, E1, E2, G12, nu12, t_ply)
        results.append((full_layup, G_xy, D11, D66))
        
    results.sort(key=lambda x: x[3], reverse=True)
    
    print("\n=== STACKING SEQUENCE OPTIMIZATION (Woven Fabric) ===")
    print(f"Base Ply E1: {E1/1e9:.2f} GPa | E2: {E2/1e9:.2f} GPa | G12: {G12/1e9:.2f} GPa")
    print("Goal: Maximize D66 (Torsional Stiffness) to prevent flutter.\n")
    
    for i in range(5):
        seq, g_xy, d11, d66 = results[i]
        print(f"Rank {i+1}: Layup {seq}")
        print(f"   D66: {d66:.2f} N-m, D11: {d11:.2f} N-m")

    # --- 4. IN-PLANE SHEAR MODULUS VERIFICATION ---
    print("\n=== IN-PLANE SHEAR MODULUS VERIFICATION ===")
    print(f"Single Woven Ply In-Plane Shear Modulus (G12): {G12/1e9:.3f} GPa")
    print(f"Full Laminate In-Plane Shear Modulus (G_xy):   {results[0][1]/1e9:.3f} GPa")