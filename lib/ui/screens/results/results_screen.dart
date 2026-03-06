import 'package:flutter/material.dart';
import 'package:go_router/go_router.dart';
import '../../../core/aeroelastic/stability_margin.dart';
import '../../../modules/flutter_analysis/flutter_output.dart';
import '../../../core/fea/solver/flutter_solver.dart';
import '../../theme/color_palette.dart';

/// Results display screen.
class ResultsScreen extends StatefulWidget {
  final String projectId;
  final FlutterAnalysisOutput? output;

  const ResultsScreen({super.key, required this.projectId, this.output});

  @override
  State<ResultsScreen> createState() => _ResultsScreenState();
}

class _ResultsScreenState extends State<ResultsScreen> {
  int _tabIndex = 0;

  @override
  Widget build(BuildContext context) {
    final output = widget.output;

    if (output == null) {
      return Scaffold(
        appBar: AppBar(title: const Text('Results')),
        body: Center(
          child: Column(
            mainAxisAlignment: MainAxisAlignment.center,
            children: [
              const Icon(Icons.hourglass_empty, size: 64, color: Colors.grey),
              const SizedBox(height: 16),
              const Text('No analysis results yet.'),
              const SizedBox(height: 16),
              FilledButton.icon(
                onPressed: () =>
                    context.go('/project/${widget.projectId}/analysis'),
                icon: const Icon(Icons.play_circle),
                label: const Text('Run Analysis'),
              ),
            ],
          ),
        ),
      );
    }

    return Scaffold(
      appBar: AppBar(
        title: const Text('Analysis Results'),
        leading: BackButton(
            onPressed: () =>
                context.go('/project/${widget.projectId}/analysis')),
        actions: [
          IconButton(
            icon: const Icon(Icons.tune),
            onPressed: () =>
                context.go('/project/${widget.projectId}/optimization'),
            tooltip: 'Optimize Design',
          ),
        ],
      ),
      body: Column(
        children: [
          // Summary cards row
          Padding(
            padding: const EdgeInsets.all(16),
            child: Row(
              children: [
                Expanded(child: _FlutterSpeedCard(output: output)),
                const SizedBox(width: 16),
                Expanded(child: _DivergenceSpeedCard(output: output)),
              ],
            ),
          ),
          // Tabs
          TabBar(
            tabs: const [
              Tab(text: 'V-g Diagram'),
              Tab(text: 'Mode Shapes'),
              Tab(text: 'Summary'),
            ],
            onTap: (i) => setState(() => _tabIndex = i),
          ),
          Expanded(
            child: IndexedStack(
              index: _tabIndex,
              children: [
                _VGDiagramTab(output: output),
                _ModeShapesTab(output: output),
                _SummaryTab(output: output),
              ],
            ),
          ),
        ],
      ),
    );
  }
}

class _FlutterSpeedCard extends StatelessWidget {
  final FlutterAnalysisOutput output;

  const _FlutterSpeedCard({required this.output});

  @override
  Widget build(BuildContext context) {
    final vF = output.flutterResult.flutterSpeed;
    final margin = output.flutterMargin;
    final status = output.flutterStatus;

    final statusColor = switch (status) {
      SafetyStatus.safe => ColorPalette.safe,
      SafetyStatus.warning => ColorPalette.warning,
      SafetyStatus.marginal => ColorPalette.marginal,
      SafetyStatus.critical => ColorPalette.critical,
      null => Colors.grey,
    };

    return Card(
      child: Padding(
        padding: const EdgeInsets.all(16),
        child: Column(
          crossAxisAlignment: CrossAxisAlignment.start,
          children: [
            Row(
              mainAxisAlignment: MainAxisAlignment.spaceBetween,
              children: [
                Text('Flutter Speed',
                    style: Theme.of(context).textTheme.titleMedium),
                if (status != null)
                  Chip(
                    label: Text(status.label,
                        style: const TextStyle(color: Colors.white)),
                    backgroundColor: statusColor,
                    padding: EdgeInsets.zero,
                  ),
              ],
            ),
            const SizedBox(height: 8),
            Text(
              vF != null ? '${vF.toStringAsFixed(0)} m/s' : 'No flutter',
              style: Theme.of(context)
                  .textTheme
                  .headlineMedium
                  ?.copyWith(color: statusColor),
            ),
            if (margin != null) ...[
              const SizedBox(height: 4),
              Text('Margin: ×${margin.toStringAsFixed(2)}',
                  style: TextStyle(color: Colors.grey[600])),
            ],
          ],
        ),
      ),
    );
  }
}

class _DivergenceSpeedCard extends StatelessWidget {
  final FlutterAnalysisOutput output;

  const _DivergenceSpeedCard({required this.output});

  @override
  Widget build(BuildContext context) {
    final vD = output.divergenceResult.vDivergence;
    final margin = output.divergenceMargin;

    return Card(
      child: Padding(
        padding: const EdgeInsets.all(16),
        child: Column(
          crossAxisAlignment: CrossAxisAlignment.start,
          children: [
            Text('Divergence Speed',
                style: Theme.of(context).textTheme.titleMedium),
            const SizedBox(height: 8),
            Text(
              vD != null ? '${vD.toStringAsFixed(0)} m/s' : 'No divergence',
              style: Theme.of(context)
                  .textTheme
                  .headlineMedium
                  ?.copyWith(color: ColorPalette.primary),
            ),
            if (margin != null) ...[
              const SizedBox(height: 4),
              Text('Margin: ×${margin.toStringAsFixed(2)}',
                  style: TextStyle(color: Colors.grey[600])),
            ],
          ],
        ),
      ),
    );
  }
}

class _VGDiagramTab extends StatelessWidget {
  final FlutterAnalysisOutput output;

  const _VGDiagramTab({required this.output});

  @override
  Widget build(BuildContext context) {
    final vgCurves = output.flutterResult.vgCurves;
    if (vgCurves.isEmpty || vgCurves.every((c) => c.isEmpty)) {
      return const Center(child: Text('No V-g data available'));
    }

    return Padding(
      padding: const EdgeInsets.all(16),
      child: Column(
        crossAxisAlignment: CrossAxisAlignment.start,
        children: [
          Text('V-g Diagram (Velocity vs. Structural Damping)',
              style: Theme.of(context).textTheme.titleMedium),
          const SizedBox(height: 8),
          Text(
            'Flutter boundary: damping g crosses zero from negative to positive.',
            style: TextStyle(color: Colors.grey[600], fontSize: 12),
          ),
          const SizedBox(height: 16),
          Expanded(
            child: CustomPaint(
              painter: _VGPainter(vgCurves: vgCurves),
              size: Size.infinite,
            ),
          ),
          // Legend
          Wrap(
            spacing: 16,
            children: vgCurves.asMap().entries.map((e) {
              final color = ColorPalette.modeColors[
                  e.key % ColorPalette.modeColors.length];
              return Row(
                mainAxisSize: MainAxisSize.min,
                children: [
                  Container(width: 20, height: 3, color: color),
                  const SizedBox(width: 4),
                  Text('Mode ${e.key + 1}',
                      style: const TextStyle(fontSize: 12)),
                ],
              );
            }).toList(),
          ),
        ],
      ),
    );
  }
}

class _VGPainter extends CustomPainter {
  final List<List<VGPoint>> vgCurves;

  _VGPainter({required this.vgCurves});

  @override
  void paint(Canvas canvas, Size size) {
    const padding = 40.0;
    final plotWidth = size.width - 2 * padding;
    final plotHeight = size.height - 2 * padding;

    if (plotWidth <= 0 || plotHeight <= 0) return;

    // Find data range
    double minV = double.infinity, maxV = -double.infinity;
    double minG = double.infinity, maxG = -double.infinity;

    for (final curve in vgCurves) {
      for (final pt in curve) {
        if (pt.velocity < minV) minV = pt.velocity;
        if (pt.velocity > maxV) maxV = pt.velocity;
        if (pt.damping < minG) minG = pt.damping;
        if (pt.damping > maxG) maxG = pt.damping;
      }
    }

    if (minV == double.infinity) return;

    final vRange = (maxV - minV).clamp(1.0, double.infinity);
    final gRange = (maxG - minG).clamp(0.1, double.infinity);

    Offset toCanvas(double v, double g) => Offset(
      padding + (v - minV) / vRange * plotWidth,
      padding + plotHeight - (g - minG) / gRange * plotHeight,
    );

    // Draw axes
    final axisPaint = Paint()
      ..color = Colors.grey
      ..strokeWidth = 1;
    canvas.drawLine(Offset(padding, padding),
        Offset(padding, padding + plotHeight), axisPaint);
    canvas.drawLine(Offset(padding, padding + plotHeight),
        Offset(padding + plotWidth, padding + plotHeight), axisPaint);

    // Zero damping line
    if (minG < 0 && maxG > 0) {
      final zeroY = padding + plotHeight - (0 - minG) / gRange * plotHeight;
      canvas.drawLine(
        Offset(padding, zeroY),
        Offset(padding + plotWidth, zeroY),
        Paint()
          ..color = Colors.red.withOpacity(0.5)
          ..strokeWidth = 1
          ..style = PaintingStyle.stroke,
      );
    }

    // Draw curves
    for (var i = 0; i < vgCurves.length; i++) {
      final curve = vgCurves[i];
      if (curve.length < 2) continue;
      final color = ColorPalette.modeColors[i % ColorPalette.modeColors.length];
      final paint = Paint()
        ..color = color
        ..strokeWidth = 2
        ..style = PaintingStyle.stroke;

      final path = Path();
      path.moveTo(toCanvas(curve[0].velocity, curve[0].damping).dx,
          toCanvas(curve[0].velocity, curve[0].damping).dy);
      for (var j = 1; j < curve.length; j++) {
        final p = toCanvas(curve[j].velocity, curve[j].damping);
        path.lineTo(p.dx, p.dy);
      }
      canvas.drawPath(path, paint);
    }
  }

  @override
  bool shouldRepaint(_VGPainter old) => true;
}

class _ModeShapesTab extends StatelessWidget {
  final FlutterAnalysisOutput output;

  const _ModeShapesTab({required this.output});

  @override
  Widget build(BuildContext context) {
    final freqs = output.feaResult.naturalFrequenciesHz;
    return Padding(
      padding: const EdgeInsets.all(16),
      child: Column(
        crossAxisAlignment: CrossAxisAlignment.start,
        children: [
          Text('Natural Frequencies',
              style: Theme.of(context).textTheme.titleMedium),
          const SizedBox(height: 12),
          ...freqs.asMap().entries.map((e) => Card(
                child: ListTile(
                  leading: CircleAvatar(
                    backgroundColor: ColorPalette.modeColors[
                        e.key % ColorPalette.modeColors.length],
                    child: Text('${e.key + 1}',
                        style: const TextStyle(color: Colors.white)),
                  ),
                  title: Text('Mode ${e.key + 1}'),
                  subtitle: Text(
                      'f = ${e.value.toStringAsFixed(2)} Hz  '
                      '(ω = ${(e.value * 2 * 3.14159).toStringAsFixed(1)} rad/s)'),
                ),
              )),
        ],
      ),
    );
  }
}

class _SummaryTab extends StatelessWidget {
  final FlutterAnalysisOutput output;

  const _SummaryTab({required this.output});

  @override
  Widget build(BuildContext context) {
    return SingleChildScrollView(
      padding: const EdgeInsets.all(16),
      child: Column(
        crossAxisAlignment: CrossAxisAlignment.start,
        children: [
          _section(context, 'Structural', [
            _row('FEA Mode Count', '${output.feaResult.nModes}'),
            _row('1st Natural Freq.',
                '${output.feaResult.naturalFrequenciesHz.firstOrNull?.toStringAsFixed(2) ?? "--"} Hz'),
          ]),
          const SizedBox(height: 16),
          _section(context, 'Aerodynamic', [
            _row('Lift Coefficient CL',
                output.cfdResult.CL.toStringAsFixed(4)),
            _row('Compressibility Correction',
                output.cfdResult.compressibilityCorrectionApplied ? 'Applied' : 'Not applied'),
          ]),
          const SizedBox(height: 16),
          _section(context, 'Aeroelastic', [
            _row('Flutter Speed',
                output.flutterResult.flutterSpeed?.toStringAsFixed(1) != null
                    ? '${output.flutterResult.flutterSpeed!.toStringAsFixed(1)} m/s'
                    : 'None in sweep'),
            _row('Flutter Frequency',
                output.flutterResult.flutterFrequency?.toStringAsFixed(2) != null
                    ? '${output.flutterResult.flutterFrequency!.toStringAsFixed(2)} Hz'
                    : '--'),
            _row('Divergence Speed',
                output.divergenceResult.vDivergence?.toStringAsFixed(1) != null
                    ? '${output.divergenceResult.vDivergence!.toStringAsFixed(1)} m/s'
                    : 'None found'),
            _row('Flutter Margin',
                output.flutterMargin?.toStringAsFixed(3) ?? '--'),
          ]),
        ],
      ),
    );
  }

  Widget _section(BuildContext ctx, String title, List<Widget> rows) {
    return Card(
      child: Padding(
        padding: const EdgeInsets.all(16),
        child: Column(
          crossAxisAlignment: CrossAxisAlignment.start,
          children: [
            Text(title, style: Theme.of(ctx).textTheme.titleMedium),
            const Divider(),
            ...rows,
          ],
        ),
      ),
    );
  }

  Widget _row(String label, String value) => Padding(
        padding: const EdgeInsets.symmetric(vertical: 4),
        child: Row(
          mainAxisAlignment: MainAxisAlignment.spaceBetween,
          children: [
            Text(label, style: const TextStyle(color: Colors.grey)),
            Text(value, style: const TextStyle(fontWeight: FontWeight.w500)),
          ],
        ),
      );
}
