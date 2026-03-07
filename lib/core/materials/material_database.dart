import 'orthotropic_ply.dart';
import 'isotropic_material.dart';

/// Built-in material property database.
class MaterialDatabase {
  // --- Composite Prepregs ---
  static const as43501 = OrthotropicPly(
    name: 'AS4/3501-6',
    E1: 142e9,
    E2: 10.3e9,
    G12: 7.2e9,
    nu12: 0.27,
    rho: 1580.0,
    t: 125e-6,
    xt: 2280e6,
    xc: 1440e6,
    yt: 57e6,
    yc: 228e6,
    s12: 71e6,
  );

  static const t3005208 = OrthotropicPly(
    name: 'T300/5208',
    E1: 132e9,
    E2: 10.8e9,
    G12: 5.65e9,
    nu12: 0.28,
    rho: 1600.0,
    t: 125e-6,
    xt: 1500e6,
    xc: 1200e6,
    yt: 40e6,
    yc: 168e6,
    s12: 68e6,
  );

  static const im78552 = OrthotropicPly(
    name: 'IM7/8552',
    E1: 171e9,
    E2: 9.08e9,
    G12: 5.29e9,
    nu12: 0.32,
    rho: 1570.0,
    t: 131e-6,
    xt: 2326e6,
    xc: 1200e6,
    yt: 62e6,
    yc: 200e6,
    s12: 92e6,
  );

  static const eglassEpoxy = OrthotropicPly(
    name: 'E-glass/Epoxy',
    E1: 45e9,
    E2: 12e9,
    G12: 5.5e9,
    nu12: 0.28,
    rho: 2100.0,
    t: 250e-6,
    xt: 1020e6,
    xc: 620e6,
    yt: 40e6,
    yc: 140e6,
    s12: 60e6,
  );

  static const kevlarEpoxy = OrthotropicPly(
    name: 'Kevlar/Epoxy',
    E1: 76e9,
    E2: 5.5e9,
    G12: 2.3e9,
    nu12: 0.34,
    rho: 1380.0,
    t: 130e-6,
    xt: 1380e6,
    xc: 240e6,
    yt: 30e6,
    yc: 138e6,
    s12: 49e6,
  );

  static const cfrpGeneric = OrthotropicPly(
    name: 'CFRP Generic',
    E1: 120e9,
    E2: 8e9,
    G12: 4.5e9,
    nu12: 0.30,
    rho: 1550.0,
    t: 150e-6,
    xt: 1600e6,
    xc: 1100e6,
    yt: 45e6,
    yc: 180e6,
    s12: 65e6,
  );

  /// All built-in composite plies.
  static const List<OrthotropicPly> compositePlies = [
    as43501,
    t3005208,
    im78552,
    eglassEpoxy,
    kevlarEpoxy,
    cfrpGeneric,
  ];

  /// All built-in isotropic materials.
  static const List<IsotropicMaterial> isotropicMaterials = [
    IsotropicMaterial.aluminum2024T3,
    IsotropicMaterial.steel4340,
    IsotropicMaterial.titaniumTi6Al4V,
  ];

  /// Lookup composite ply by name.
  static OrthotropicPly? findCompositePly(String name) {
    for (final ply in compositePlies) {
      if (ply.name == name) return ply;
    }
    return null;
  }

  /// Lookup isotropic material by name.
  static IsotropicMaterial? findIsotropicMaterial(String name) {
    for (final mat in isotropicMaterials) {
      if (mat.name == name) return mat;
    }
    return null;
  }
}
