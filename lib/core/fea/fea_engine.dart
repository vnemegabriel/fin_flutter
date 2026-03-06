import '../math/matrix.dart';
import '../materials/clt_calculator.dart';
import 'assembly/global_assembler.dart';
import 'mesh/boundary_condition.dart';
import 'mesh/mesh.dart';
import 'mesh/mesh_generator.dart';
import 'results/fea_result.dart';
import 'solver/eigenvalue_solver.dart';
import '../models/fin_geometry.dart';

/// Configuration for an FEA analysis.
class FEAConfig {
  final int spanwiseElements;
  final int chordwiseElements;
  final ElementType elementType;
  final int nModes;

  const FEAConfig({
    this.spanwiseElements = 8,
    this.chordwiseElements = 6,
    this.elementType = ElementType.kirchhoffDKQ,
    this.nModes = 6,
  });
}

/// Orchestrates FEA: mesh generation, assembly, eigenvalue solve.
class FEAEngine {
  final FEAConfig config;

  const FEAEngine({this.config = const FEAConfig()});

  /// Run full FEA analysis.
  ///
  /// [geometry] - fin geometry
  /// [abd] - laminate ABD matrices
  FEAResult analyze({
    required FinGeometry geometry,
    required LaminateABD abd,
  }) {
    // 1. Generate mesh
    final meshGen = MeshGenerator();
    final mesh = meshGen.generateFinMesh(
      geometry: geometry,
      spanwiseElements: config.spanwiseElements,
      chordwiseElements: config.chordwiseElements,
      elementType: config.elementType,
    );

    // 2. Apply boundary conditions (clamped root)
    final bc = BoundaryConditions.clampedRoot(mesh);

    // 3. Assemble global matrices
    final assembler = GlobalAssembler();
    final ({Matrix K, Matrix M}) matrices;

    if (config.elementType == ElementType.kirchhoffDKQ) {
      matrices = assembler.assembleKirchhoff(mesh: mesh, abd: abd, bc: bc);
    } else {
      matrices = assembler.assembleMindlin(mesh: mesh, abd: abd, bc: bc);
    }

    // 4. Solve eigenvalue problem
    final eigSolver = EigenvalueSolver();
    final eigResult = eigSolver.solveGeneralized(
      K: matrices.K,
      M: matrices.M,
      nModes: config.nModes,
    );

    // 5. Compute natural frequencies (filter negative eigenvalues from BC penalty)
    final validModes = <int>[];
    for (var i = 0; i < eigResult.eigenvalues.length; i++) {
      final ev = eigResult.eigenvalues[i];
      if (ev > 1.0) validModes.add(i); // filter BC penalty modes
    }

    final freqs = validModes.map((i) => eigResult.eigenvalues[i]).toList();
    // Convert omega² -> omega
    final omegaList = freqs.map((ev) => ev > 0 ? ev.abs() : 0.0).toList();

    final nModes = validModes.length;
    final modeShapes = Matrix(mesh.totalDof, nModes);
    for (var j = 0; j < nModes; j++) {
      final modeIdx = validModes[j];
      for (var i = 0; i < mesh.totalDof; i++) {
        modeShapes.set(i, j, eigResult.eigenvectors.get(i, modeIdx));
      }
    }

    // Compute generalized masses
    final genMasses = List<double>.generate(nModes, (j) {
      final phi = modeShapes.column(j);
      final Mphi = matrices.M.multiplyVector(phi);
      return _dotProduct(phi, Mphi);
    });

    return FEAResult(
      nModes: nModes,
      naturalFrequenciesRadPS: omegaList,
      modeShapes: modeShapes,
      generalizedMasses: genMasses,
    );
  }

  double _dotProduct(List<double> a, List<double> b) {
    double s = 0;
    for (var i = 0; i < a.length; i++) s += a[i] * b[i];
    return s;
  }
}
