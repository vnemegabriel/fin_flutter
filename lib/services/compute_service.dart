import 'dart:isolate';
import '../modules/flutter_analysis/flutter_analysis_module.dart';
import '../modules/flutter_analysis/flutter_input.dart';
import '../modules/flutter_analysis/flutter_output.dart';

/// Dispatches heavy computation to Dart isolates.
class ComputeService {
  /// Run a full flutter analysis in a background isolate.
  ///
  /// This prevents UI jank during long computations.
  static Future<FlutterAnalysisOutput> runFlutterAnalysis(
      FlutterAnalysisInput input) async {
    return await Isolate.run(() {
      final module = FlutterAnalysisModule();
      return module.analyze(input);
    });
  }

  /// Run any computation in a background isolate.
  static Future<R> run<R>(R Function() computation) async {
    return await Isolate.run(computation);
  }
}
