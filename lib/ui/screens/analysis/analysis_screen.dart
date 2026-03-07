import 'package:flutter/material.dart';
import 'package:go_router/go_router.dart';
import '../../../core/fea/fea_engine.dart';
import '../../../core/fea/mesh/mesh.dart';
import '../../../core/cfd/cfd_engine.dart';
import '../../../core/materials/material_database.dart';
import '../../../core/materials/clt_calculator.dart';
import '../../../core/materials/laminate.dart';
import '../../../core/models/fin_geometry.dart';
import '../../../core/models/flight_condition.dart';
import '../../../modules/flutter_analysis/flutter_analysis_module.dart';
import '../../../modules/flutter_analysis/flutter_input.dart';
import '../../../modules/flutter_analysis/flutter_output.dart';
import '../../widgets/loading_overlay.dart';

/// Analysis configuration and run screen.
class AnalysisScreen extends StatefulWidget {
  final String projectId;

  const AnalysisScreen({super.key, required this.projectId});

  @override
  State<AnalysisScreen> createState() => _AnalysisScreenState();
}

class _AnalysisScreenState extends State<AnalysisScreen> {
  int _spanwiseElements = 6;
  int _chordwiseElements = 4;
  int _chordwisePanels = 6;
  int _spanwisePanels = 8;
  int _nModes = 4;
  ElementType _elementType = ElementType.kirchhoffDKQ;
  bool _running = false;
  String? _error;
  FlutterAnalysisOutput? _result;

  Future<void> _runAnalysis() async {
    setState(() {
      _running = true;
      _error = null;
    });

    try {
      // Example default inputs (in a real app, these come from state/providers)
      final geometry = FinGeometry(
        span: 0.15,
        rootChord: 0.12,
        tipChord: 0.06,
        sweepLength: 0.04,
        thickness: 0.003,
      );
      final ply = MaterialDatabase.as43501;
      final lam = Laminate.quasiIsotropic(ply);
      final abd = lam.computeABD();
      final condition = FlightCondition.isa(altitude: 1000, mach: 0.5);

      final input = FlutterAnalysisInput(
        geometry: geometry,
        laminate: abd,
        flightCondition: condition,
        feaConfig: FEAConfig(
          spanwiseElements: _spanwiseElements,
          chordwiseElements: _chordwiseElements,
          elementType: _elementType,
          nModes: _nModes,
        ),
        cfdConfig: CFDConfig(
          chordwisePanels: _chordwisePanels,
          spanwisePanels: _spanwisePanels,
        ),
        velocityMin: 20.0,
        velocityMax: 500.0,
        velocitySteps: 40,
        maxFlightVelocity: condition.velocity,
      );

      // Run in isolate
      final result = await Future(() {
        return FlutterAnalysisModule().analyze(input);
      });

      if (mounted) {
        setState(() {
          _result = result;
          _running = false;
        });
        context.go('/project/${widget.projectId}/results');
      }
    } catch (e) {
      if (mounted) {
        setState(() {
          _error = e.toString();
          _running = false;
        });
      }
    }
  }

  @override
  Widget build(BuildContext context) {
    return Stack(
      children: [
        Scaffold(
          appBar: AppBar(
            title: const Text('Analysis Setup'),
            leading: BackButton(
                onPressed: () =>
                    context.go('/project/${widget.projectId}/flight')),
          ),
          body: SingleChildScrollView(
            padding: const EdgeInsets.all(24),
            child: ConstrainedBox(
              constraints: const BoxConstraints(maxWidth: 600),
              child: Column(
                crossAxisAlignment: CrossAxisAlignment.start,
                children: [
                  if (_error != null) ...[
                    Card(
                      color: Theme.of(context).colorScheme.errorContainer,
                      child: Padding(
                        padding: const EdgeInsets.all(12),
                        child: Text(_error!),
                      ),
                    ),
                    const SizedBox(height: 16),
                  ],
                  // FEA config
                  Text('FEA Configuration',
                      style: Theme.of(context).textTheme.titleMedium),
                  const SizedBox(height: 12),
                  DropdownButtonFormField<ElementType>(
                    value: _elementType,
                    items: const [
                      DropdownMenuItem(
                          value: ElementType.kirchhoffDKQ,
                          child: Text('Kirchhoff DKQ (thin plate)')),
                      DropdownMenuItem(
                          value: ElementType.mindlinMITC4,
                          child: Text('Mindlin MITC4 (thick plate)')),
                    ],
                    onChanged: (v) => setState(() => _elementType = v!),
                    decoration: const InputDecoration(labelText: 'Element Type'),
                  ),
                  const SizedBox(height: 12),
                  _intSlider('Spanwise Elements', _spanwiseElements, 2, 20,
                      (v) => setState(() => _spanwiseElements = v)),
                  _intSlider('Chordwise Elements', _chordwiseElements, 2, 12,
                      (v) => setState(() => _chordwiseElements = v)),
                  _intSlider('Modal Modes', _nModes, 2, 10,
                      (v) => setState(() => _nModes = v)),
                  const SizedBox(height: 24),
                  // CFD config
                  Text('CFD/VLM Configuration',
                      style: Theme.of(context).textTheme.titleMedium),
                  const SizedBox(height: 12),
                  _intSlider('Chordwise Panels', _chordwisePanels, 2, 20,
                      (v) => setState(() => _chordwisePanels = v)),
                  _intSlider('Spanwise Panels', _spanwisePanels, 2, 30,
                      (v) => setState(() => _spanwisePanels = v)),
                  const SizedBox(height: 32),
                  // Complexity estimate
                  _buildComplexityCard(),
                  const SizedBox(height: 24),
                  SizedBox(
                    width: double.infinity,
                    child: FilledButton.icon(
                      onPressed: _running ? null : _runAnalysis,
                      icon: const Icon(Icons.play_circle),
                      label: const Text('Run Analysis'),
                    ),
                  ),
                ],
              ),
            ),
          ),
        ),
        if (_running) const LoadingOverlay(message: 'Running FEA + VLM analysis...'),
      ],
    );
  }

  Widget _intSlider(String label, int value, int min, int max, void Function(int) onChanged) {
    return Column(
      crossAxisAlignment: CrossAxisAlignment.start,
      children: [
        Row(
          mainAxisAlignment: MainAxisAlignment.spaceBetween,
          children: [
            Text(label),
            Text('$value', style: const TextStyle(fontWeight: FontWeight.bold)),
          ],
        ),
        Slider(
          value: value.toDouble(),
          min: min.toDouble(),
          max: max.toDouble(),
          divisions: max - min,
          onChanged: (v) => onChanged(v.round()),
        ),
      ],
    );
  }

  Widget _buildComplexityCard() {
    final totalNodes = (_spanwiseElements + 1) * (_chordwiseElements + 1);
    final totalDof = totalNodes * 3;
    final totalPanels = _chordwisePanels * _spanwisePanels;

    return Card(
      child: Padding(
        padding: const EdgeInsets.all(12),
        child: Column(
          crossAxisAlignment: CrossAxisAlignment.start,
          children: [
            Text('Model Size Estimate',
                style: Theme.of(context).textTheme.titleSmall),
            const SizedBox(height: 8),
            _infoRow('FEA Nodes', '$totalNodes'),
            _infoRow('Structural DOF', '$totalDof'),
            _infoRow('VLM Panels', '$totalPanels'),
            _infoRow('AIC Matrix Size', '${totalPanels}×${totalPanels}'),
          ],
        ),
      ),
    );
  }

  Widget _infoRow(String label, String value) => Padding(
        padding: const EdgeInsets.symmetric(vertical: 2),
        child: Row(
          mainAxisAlignment: MainAxisAlignment.spaceBetween,
          children: [
            Text(label, style: const TextStyle(color: Colors.grey)),
            Text(value),
          ],
        ),
      );
}
