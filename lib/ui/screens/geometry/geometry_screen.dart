import 'package:flutter/material.dart';
import 'package:go_router/go_router.dart';
import '../../../core/models/fin_geometry.dart';
import '../../widgets/numerical_input_field.dart';
import '../../widgets/fin_canvas.dart';

/// Geometry input screen for fin parameters.
class GeometryScreen extends StatefulWidget {
  final String projectId;

  const GeometryScreen({super.key, required this.projectId});

  @override
  State<GeometryScreen> createState() => _GeometryScreenState();
}

class _GeometryScreenState extends State<GeometryScreen> {
  double _span = 0.15;
  double _rootChord = 0.12;
  double _tipChord = 0.06;
  double _sweepLength = 0.04;
  double _thickness = 0.003;

  FinGeometry get _geometry => FinGeometry(
        span: _span,
        rootChord: _rootChord,
        tipChord: _tipChord,
        sweepLength: _sweepLength,
        thickness: _thickness,
      );

  @override
  Widget build(BuildContext context) {
    return Scaffold(
      appBar: AppBar(
        title: const Text('Fin Geometry'),
        leading: BackButton(onPressed: () => context.go('/')),
      ),
      body: Row(
        crossAxisAlignment: CrossAxisAlignment.start,
        children: [
          // Input panel
          SizedBox(
            width: 340,
            child: SingleChildScrollView(
              padding: const EdgeInsets.all(16),
              child: Column(
                crossAxisAlignment: CrossAxisAlignment.start,
                children: [
                  Text('Planform Dimensions',
                      style: Theme.of(context).textTheme.titleMedium),
                  const SizedBox(height: 12),
                  NumericalInputField(
                    label: 'Span',
                    unit: 'm',
                    value: _span,
                    min: 0.01,
                    max: 1.0,
                    onChanged: (v) => setState(() => _span = v ?? _span),
                  ),
                  const SizedBox(height: 8),
                  NumericalInputField(
                    label: 'Root Chord',
                    unit: 'm',
                    value: _rootChord,
                    min: 0.01,
                    max: 0.5,
                    onChanged: (v) => setState(() => _rootChord = v ?? _rootChord),
                  ),
                  const SizedBox(height: 8),
                  NumericalInputField(
                    label: 'Tip Chord',
                    unit: 'm',
                    value: _tipChord,
                    min: 0.0,
                    max: 0.5,
                    onChanged: (v) => setState(() => _tipChord = v ?? _tipChord),
                  ),
                  const SizedBox(height: 8),
                  NumericalInputField(
                    label: 'Sweep Length',
                    unit: 'm',
                    value: _sweepLength,
                    min: 0.0,
                    max: 0.5,
                    onChanged: (v) => setState(() => _sweepLength = v ?? _sweepLength),
                  ),
                  const SizedBox(height: 8),
                  NumericalInputField(
                    label: 'Thickness',
                    unit: 'm',
                    value: _thickness,
                    min: 0.0005,
                    max: 0.05,
                    onChanged: (v) => setState(() => _thickness = v ?? _thickness),
                  ),
                  const SizedBox(height: 24),
                  // Derived properties
                  _buildDerivedProperties(),
                  const SizedBox(height: 24),
                  SizedBox(
                    width: double.infinity,
                    child: FilledButton.icon(
                      onPressed: () =>
                          context.go('/project/${widget.projectId}/materials'),
                      icon: const Icon(Icons.arrow_forward),
                      label: const Text('Next: Materials'),
                    ),
                  ),
                ],
              ),
            ),
          ),
          const VerticalDivider(),
          // Preview panel
          Expanded(
            child: Center(
              child: Column(
                mainAxisAlignment: MainAxisAlignment.center,
                children: [
                  Text('Fin Planform Preview',
                      style: Theme.of(context).textTheme.titleMedium),
                  const SizedBox(height: 16),
                  FinCanvas(
                    geometry: _geometry,
                    width: 350,
                    height: 250,
                  ),
                ],
              ),
            ),
          ),
        ],
      ),
    );
  }

  Widget _buildDerivedProperties() {
    final geo = _geometry;
    return Card(
      child: Padding(
        padding: const EdgeInsets.all(12),
        child: Column(
          crossAxisAlignment: CrossAxisAlignment.start,
          children: [
            Text('Derived Properties',
                style: Theme.of(context).textTheme.titleSmall),
            const SizedBox(height: 8),
            _propRow('Planform Area', '${(geo.planformArea * 1e4).toStringAsFixed(1)} cm²'),
            _propRow('Aspect Ratio', geo.aspectRatio.toStringAsFixed(2)),
            _propRow('Sweep Angle', '${geo.sweepAngleDeg.toStringAsFixed(1)}°'),
            _propRow('Taper Ratio', geo.taperRatio.toStringAsFixed(3)),
            _propRow('MAC', '${(geo.meanAerodynamicChord * 100).toStringAsFixed(1)} cm'),
          ],
        ),
      ),
    );
  }

  Widget _propRow(String label, String value) {
    return Padding(
      padding: const EdgeInsets.symmetric(vertical: 2),
      child: Row(
        mainAxisAlignment: MainAxisAlignment.spaceBetween,
        children: [
          Text(label, style: const TextStyle(color: Colors.grey)),
          Text(value, style: const TextStyle(fontWeight: FontWeight.w500)),
        ],
      ),
    );
  }
}
