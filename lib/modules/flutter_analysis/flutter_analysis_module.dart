import '../../core/fea/fea_engine.dart';
import '../../core/cfd/cfd_engine.dart';
import '../../core/aeroelastic/aeroelastic_coupler.dart';
import '../../core/fea/solver/flutter_solver.dart';
import '../../core/fea/solver/divergence_solver.dart';
import '../../core/math/matrix.dart';
import '../../core/models/flight_condition.dart';
import 'flutter_input.dart';
import 'flutter_output.dart';

/// Orchestrates the complete flutter and divergence analysis pipeline.
///
/// Pipeline:
/// 1. FEA: mesh + assemble + eigensolve → mode shapes
/// 2. CFD/VLM: discretize fin + solve → AIC matrix
/// 3. Aeroelastic coupling: project AIC onto mode shapes → Q_modal
/// 4. Flutter: sweep velocity, solve U-g → V_flutter
/// 5. Divergence: solve static aeroelastic eigenvalue → V_divergence
class FlutterAnalysisModule {
  /// Run full analysis. This may take several seconds and should be called
  /// in an isolate via ComputeService.
  FlutterAnalysisOutput analyze(FlutterAnalysisInput input) {
    // 1. FEA
    final feaEngine = FEAEngine(config: input.feaConfig);
    final feaResult = feaEngine.analyze(
      geometry: input.geometry,
      abd: input.laminate,
    );

    // 2. CFD
    final cfdEngine = CFDEngine(config: input.cfdConfig);
    final cfdResult = cfdEngine.analyze(
      geometry: input.geometry,
      condition: input.flightCondition,
    );

    // 3. Aeroelastic coupling
    final coupler = AeroelasticCoupler();
    final QModal = coupler.buildModalAeroMatrix(
      feaResult: feaResult,
      vlmResult: cfdResult.vlmResult,
      geometry: input.geometry,
      condition: input.flightCondition,
    );

    // Build modal stiffness and mass (diagonal matrices)
    final nModes = feaResult.nModes;
    final KModal = Matrix(nModes, nModes);
    final MModal = Matrix(nModes, nModes);
    for (var i = 0; i < nModes; i++) {
      final omega = feaResult.naturalFrequenciesRadPS[i];
      KModal.set(i, i, omega * omega * feaResult.generalizedMasses[i]);
      MModal.set(i, i, feaResult.generalizedMasses[i]);
    }

    // 4. Flutter analysis
    final flutterSolver = FlutterSolver();
    final flutterResult = flutterSolver.solveUG(
      kModal: KModal,
      mModal: MModal,
      qModal: QModal,
      rho: input.flightCondition.density,
      velocities: input.velocitySweep,
    );

    // 5. Divergence analysis
    // Build static aerodynamic stiffness matrix projected onto modal coordinates
    final nDof = feaResult.modeShapes.rows;
    final AStatic = coupler.buildStaticAeroMatrix(
      vlmResult: cfdResult.vlmResult,
      geometry: input.geometry,
      nDof: nDof,
      dofPerNode: input.feaConfig.elementType.index == 0 ? 3 : 5,
    );

    // Project onto modal space
    final Phi = feaResult.modeShapes;
    final AStaticModal = Phi.transpose() * AStatic * Phi;

    final divergenceSolver = DivergenceSolver();
    final divergenceResult = divergenceSolver.solve(
      K: KModal,
      aStatic: AStaticModal,
      rho: input.flightCondition.density,
    );

    return FlutterAnalysisOutput(
      feaResult: feaResult,
      cfdResult: cfdResult,
      flutterResult: flutterResult,
      divergenceResult: divergenceResult,
      maxFlightVelocity: input.maxFlightVelocity,
    );
  }
}
