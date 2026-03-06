import 'dart:math';
import '../../math/matrix.dart';
import '../../models/fin_geometry.dart';
import '../../models/flight_condition.dart';
import 'horseshoe_vortex.dart';

/// A VLM panel (horseshoe vortex).
class VLMPanel {
  final int id;
  // Bound vortex endpoints (1/4-chord line)
  final Vec3 boundA;  // inboard end
  final Vec3 boundB;  // outboard end
  // Control point (3/4-chord)
  final Vec3 controlPoint;
  // Normal vector at control point
  final Vec3 normal;
  // Panel area (m²)
  final double area;
  // Spanwise width (m)
  final double spanWidth;

  const VLMPanel({
    required this.id,
    required this.boundA,
    required this.boundB,
    required this.controlPoint,
    required this.normal,
    required this.area,
    required this.spanWidth,
  });
}

/// Result of VLM solve.
class VLMResult {
  /// Vortex strengths (circulation, m²/s) for each panel.
  final List<double> gamma;

  /// Lift coefficient per panel.
  final List<double> panelCL;

  /// Total lift coefficient.
  final double CL;

  /// Total drag coefficient (induced).
  final double CDi;

  /// Aerodynamic influence coefficient matrix [A].
  final Matrix aicMatrix;

  /// All panels.
  final List<VLMPanel> panels;

  const VLMResult({
    required this.gamma,
    required this.panelCL,
    required this.CL,
    required this.CDi,
    required this.aicMatrix,
    required this.panels,
  });
}

/// Vortex Lattice Method solver for fin aerodynamics.
class VortexLattice {
  /// Solve VLM for a fin.
  VLMResult solve({
    required FinGeometry geometry,
    required FlightCondition condition,
    required int chordwisePanels,
    required int spanwisePanels,
    required double angleOfAttackRad,
  }) {
    final nx = chordwisePanels;
    final ny = spanwisePanels;
    final nPanels = nx * ny;

    // Build panels
    final panels = _buildPanels(geometry, nx, ny);

    // Build AIC matrix
    final AIC = _buildAIC(panels, condition.velocity);

    // RHS: normal velocity due to freestream
    final rhs = List<double>.generate(nPanels, (i) {
      // V_inf component normal to panel = V * sin(alpha) ≈ V * alpha for small alpha
      return condition.velocity * sin(angleOfAttackRad) * panels[i].normal.z.sign;
    });

    // Correct RHS using actual panel normal
    for (var i = 0; i < nPanels; i++) {
      final n = panels[i].normal;
      // Freestream direction: [V, 0, 0] (streamwise x)
      final vInf = Vec3(condition.velocity * cos(angleOfAttackRad),
          0, condition.velocity * sin(angleOfAttackRad));
      rhs[i] = -(n.dot(vInf)); // negative: flow-through BC
    }

    // Solve for gamma
    final gamma = AIC.solve(rhs);

    // Compute forces via Kutta-Joukowski
    final panelCL = <double>[];
    double totalLift = 0.0;
    double totalArea = 0.0;

    for (var i = 0; i < nPanels; i++) {
      final span = panels[i].spanWidth;
      final L = condition.density * condition.velocity * gamma[i] * span;
      final q = condition.dynamicPressure;
      final cl = L / (q * panels[i].area);
      panelCL.add(cl);
      totalLift += L;
      totalArea += panels[i].area;
    }

    final refArea = geometry.planformArea;
    final CL = totalLift / (condition.dynamicPressure * refArea);

    // Induced drag (simplified Trefftz plane)
    double CDi = 0.0;
    for (var i = 0; i < nPanels; i++) {
      final n = panels[i].normal;
      final vDown = _inducedNormalVelocity(panels, gamma, panels[i].controlPoint);
      CDi += 0.5 * gamma[i] * vDown / refArea;
    }

    return VLMResult(
      gamma: gamma,
      panelCL: panelCL,
      CL: CL,
      CDi: CDi.abs(),
      aicMatrix: AIC,
      panels: panels,
    );
  }

  /// Build panels discretizing the fin planform.
  List<VLMPanel> _buildPanels(FinGeometry geo, int nx, int ny) {
    final panels = <VLMPanel>[];
    var panelId = 0;

    for (var j = 0; j < ny; j++) {
      final y0 = j / ny * geo.span;
      final y1 = (j + 1) / ny * geo.span;
      final yMid = (y0 + y1) / 2.0;

      final xLE0 = geo.leadingEdgeX(y0);
      final xLE1 = geo.leadingEdgeX(y1);
      final xLEMid = geo.leadingEdgeX(yMid);
      final chordMid = geo.chordAt(yMid);

      for (var i = 0; i < nx; i++) {
        final xi0 = i / nx;
        final xi1 = (i + 1) / nx;

        final xBound0 = xLEMid + (xi0 + (xi1 - xi0) / 4.0) * chordMid;
        final xBound1 = xBound0;
        final xCtrl = xLEMid + (xi0 + (xi1 - xi0) * 3.0 / 4.0) * chordMid;

        final boundA = Vec3(xBound0, y0, 0.0);
        final boundB = Vec3(xBound1, y1, 0.0);
        final cp = Vec3(xCtrl, yMid, 0.0);
        const normal = Vec3(0.0, 0.0, 1.0); // flat plate: normal is z-direction

        final area = chordMid * (xi1 - xi0) * (y1 - y0);
        final spanWidth = y1 - y0;

        panels.add(VLMPanel(
          id: panelId++,
          boundA: boundA,
          boundB: boundB,
          controlPoint: cp,
          normal: normal,
          area: area,
          spanWidth: spanWidth,
        ));
      }
    }
    return panels;
  }

  /// Build Aerodynamic Influence Coefficient matrix.
  /// AIC[i][j] = normal velocity at control point i due to unit gamma on panel j.
  Matrix _buildAIC(List<VLMPanel> panels, double Vinf) {
    final n = panels.length;
    final AIC = Matrix(n, n);
    const wakeDir = Vec3(1.0, 0.0, 0.0); // wake trails downstream (+x)

    for (var i = 0; i < n; i++) {
      final cp = panels[i].controlPoint;
      final ni = panels[i].normal;

      for (var j = 0; j < n; j++) {
        // Unit strength horseshoe vortex on panel j
        final A = panels[j].boundA;
        final B = panels[j].boundB;

        // Bound vortex (A→B)
        final vBound = BiotSavart.finiteSegment(P: cp, A: A, B: B, gamma: 1.0);

        // Trailing vortex from A (semi-infinite, going upstream → +x direction)
        final vTrailA = BiotSavart.semiInfinite(P: cp, A: A, gamma: -1.0, direction: wakeDir);

        // Trailing vortex from B (semi-infinite, going downstream)
        final vTrailB = BiotSavart.semiInfinite(P: cp, A: B, gamma: 1.0, direction: wakeDir);

        final vTotal = vBound + vTrailA + vTrailB;

        AIC.set(i, j, ni.dot(vTotal));
      }
    }
    return AIC;
  }

  /// Compute normal induced velocity at a point from all panels.
  double _inducedNormalVelocity(
      List<VLMPanel> panels, List<double> gamma, Vec3 P) {
    const wakeDir = Vec3(1.0, 0.0, 0.0);
    Vec3 v = const Vec3(0, 0, 0);
    for (var j = 0; j < panels.length; j++) {
      final A = panels[j].boundA;
      final B = panels[j].boundB;
      final g = gamma[j];

      v = v + BiotSavart.finiteSegment(P: P, A: A, B: B, gamma: g);
      v = v + BiotSavart.semiInfinite(P: P, A: A, gamma: -g, direction: wakeDir);
      v = v + BiotSavart.semiInfinite(P: P, A: B, gamma: g, direction: wakeDir);
    }
    return v.z; // z-component = normal velocity for flat plate
  }
}
