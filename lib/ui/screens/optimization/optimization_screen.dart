import 'package:flutter/material.dart';
import 'package:go_router/go_router.dart';
import '../../../modules/optimization/design_variables.dart';
import '../../../modules/optimization/nelder_mead.dart';
import '../../../modules/optimization/genetic_algorithm.dart';

/// Optimization screen.
class OptimizationScreen extends StatefulWidget {
  final String projectId;

  const OptimizationScreen({super.key, required this.projectId});

  @override
  State<OptimizationScreen> createState() => _OptimizationScreenState();
}

class _OptimizationScreenState extends State<OptimizationScreen> {
  String _optimizerType = 'nelder_mead';
  bool _running = false;
  final List<OptimizationIteration> _history = [];
  OptimizationIteration? _best;

  void _runOptimization() {
    setState(() {
      _running = true;
      _history.clear();
      _best = null;
    });

    // Simple demonstration: optimize Rosenbrock function
    // In real app: optimize flutter loss function
    final objective = (List<double> x) {
      // Rosenbrock: f(x,y) = (1-x)^2 + 100*(y-x^2)^2 (as demo)
      double sum = 0;
      for (var i = 0; i < x.length - 1; i++) {
        final xi = x[i] * 4 - 2;  // map [0,1] to [-2,2]
        final xi1 = x[i+1] * 4 - 2;
        sum += (1 - xi) * (1 - xi) + 100 * (xi1 - xi * xi) * (xi1 - xi * xi);
      }
      return sum;
    };

    if (_optimizerType == 'nelder_mead') {
      final optimizer = const NelderMead(maxIterations: 100);
      final iterations = optimizer.optimize(
        objective: objective,
        initialPoint: [0.5, 0.5, 0.5, 0.5],
      );

      for (final iter in iterations) {
        if (mounted) {
          setState(() {
            _history.add(iter);
            _best = iter;
          });
        }
        if (iter.converged) break;
      }
    } else {
      final optimizer = GeneticAlgorithm(maxGenerations: 50);
      final iterations = optimizer.optimize(objective: objective, nVariables: 4);

      for (final iter in iterations) {
        if (mounted) {
          setState(() {
            _history.add(iter);
            _best = iter;
          });
        }
      }
    }

    if (mounted) setState(() => _running = false);
  }

  @override
  Widget build(BuildContext context) {
    return Scaffold(
      appBar: AppBar(
        title: const Text('Design Optimization'),
        leading: BackButton(
            onPressed: () => context.go('/project/${widget.projectId}/results')),
      ),
      body: Row(
        crossAxisAlignment: CrossAxisAlignment.start,
        children: [
          // Left: configuration
          SizedBox(
            width: 300,
            child: SingleChildScrollView(
              padding: const EdgeInsets.all(16),
              child: Column(
                crossAxisAlignment: CrossAxisAlignment.start,
                children: [
                  Text('Optimizer', style: Theme.of(context).textTheme.titleMedium),
                  const SizedBox(height: 12),
                  SegmentedButton<String>(
                    segments: const [
                      ButtonSegment(value: 'nelder_mead', label: Text('Nelder-Mead')),
                      ButtonSegment(value: 'ga', label: Text('Genetic Alg.')),
                    ],
                    selected: {_optimizerType},
                    onSelectionChanged: (s) =>
                        setState(() => _optimizerType = s.first),
                  ),
                  const SizedBox(height: 24),
                  Text('Design Variables',
                      style: Theme.of(context).textTheme.titleMedium),
                  const SizedBox(height: 8),
                  const Text(
                    'Span, root chord, tip chord, sweep angle, ply angles and thicknesses will be optimized.',
                    style: TextStyle(color: Colors.grey, fontSize: 12),
                  ),
                  const SizedBox(height: 24),
                  SizedBox(
                    width: double.infinity,
                    child: FilledButton.icon(
                      onPressed: _running ? null : _runOptimization,
                      icon: _running
                          ? const SizedBox(
                              width: 16,
                              height: 16,
                              child: CircularProgressIndicator(strokeWidth: 2))
                          : const Icon(Icons.play_circle),
                      label: Text(_running ? 'Running...' : 'Run Optimization'),
                    ),
                  ),
                  if (_best != null) ...[
                    const SizedBox(height: 24),
                    Card(
                      child: Padding(
                        padding: const EdgeInsets.all(12),
                        child: Column(
                          crossAxisAlignment: CrossAxisAlignment.start,
                          children: [
                            Text('Best Result',
                                style: Theme.of(context).textTheme.titleSmall),
                            const SizedBox(height: 8),
                            Text('Objective: ${_best!.bestValue.toStringAsExponential(3)}'),
                            Text('Iterations: ${_history.length}'),
                            Text('Converged: ${_best!.converged}'),
                          ],
                        ),
                      ),
                    ),
                  ],
                ],
              ),
            ),
          ),
          const VerticalDivider(),
          // Right: convergence chart
          Expanded(
            child: Padding(
              padding: const EdgeInsets.all(16),
              child: Column(
                crossAxisAlignment: CrossAxisAlignment.start,
                children: [
                  Text('Convergence History',
                      style: Theme.of(context).textTheme.titleMedium),
                  const SizedBox(height: 8),
                  Expanded(
                    child: _history.isEmpty
                        ? const Center(
                            child: Text('Run optimization to see convergence'))
                        : CustomPaint(
                            painter: _ConvergencePainter(history: _history),
                            size: Size.infinite,
                          ),
                  ),
                ],
              ),
            ),
          ),
        ],
      ),
    );
  }
}

class _ConvergencePainter extends CustomPainter {
  final List<OptimizationIteration> history;

  _ConvergencePainter({required this.history});

  @override
  void paint(Canvas canvas, Size size) {
    if (history.isEmpty) return;
    const padding = 40.0;
    final w = size.width - 2 * padding;
    final h = size.height - 2 * padding;

    final values = history.map((i) => i.bestValue).toList();
    final minV = values.reduce((a, b) => a < b ? a : b);
    final maxV = values.reduce((a, b) => a > b ? a : b);
    final vRange = (maxV - minV).clamp(1e-10, double.infinity);

    final n = values.length.toDouble();

    final path = Path();
    for (var i = 0; i < values.length; i++) {
      final x = padding + i / n * w;
      final y = padding + h - (values[i] - minV) / vRange * h;
      if (i == 0) path.moveTo(x, y); else path.lineTo(x, y);
    }

    canvas.drawPath(
      path,
      Paint()
        ..color = Colors.blue
        ..strokeWidth = 2
        ..style = PaintingStyle.stroke,
    );

    // Axes
    final axisPaint = Paint()..color = Colors.grey..strokeWidth = 1;
    canvas.drawLine(Offset(padding, padding), Offset(padding, padding + h), axisPaint);
    canvas.drawLine(Offset(padding, padding + h), Offset(padding + w, padding + h), axisPaint);
  }

  @override
  bool shouldRepaint(_ConvergencePainter old) => old.history.length != history.length;
}
