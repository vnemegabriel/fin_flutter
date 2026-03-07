import 'package:flutter/material.dart';
import 'package:go_router/go_router.dart';
import '../../../core/materials/material_database.dart';
import '../../../core/materials/orthotropic_ply.dart';
import '../../../core/materials/clt_calculator.dart';
import '../../../core/materials/laminate.dart';

/// Material selection and laminate builder screen.
class MaterialScreen extends StatefulWidget {
  final String projectId;

  const MaterialScreen({super.key, required this.projectId});

  @override
  State<MaterialScreen> createState() => _MaterialScreenState();
}

class _MaterialScreenState extends State<MaterialScreen> {
  OrthotropicPly _selectedPly = MaterialDatabase.as43501;
  List<LaminatePly> _stack = [];
  double _newAngle = 0.0;

  Laminate get _laminate => Laminate(plies: _stack);
  LaminateABD? get _abd => _stack.isNotEmpty ? _laminate.computeABD() : null;

  void _addPly() {
    setState(() {
      _stack.add(LaminatePly(_selectedPly, _newAngle));
    });
  }

  void _removePly(int index) {
    setState(() => _stack.removeAt(index));
  }

  void _addStandardLaminate(String type) {
    setState(() {
      switch (type) {
        case 'ud':
          _stack = List.generate(8, (_) => LaminatePly(_selectedPly, 0.0));
        case 'crossply':
          _stack = [
            for (var i = 0; i < 4; i++) ...[
              LaminatePly(_selectedPly, 0.0),
              LaminatePly(_selectedPly, 90.0),
            ]
          ];
        case 'qi':
          _stack = [0.0, 45.0, -45.0, 90.0, 90.0, -45.0, 45.0, 0.0]
              .map((a) => LaminatePly(_selectedPly, a))
              .toList();
      }
    });
  }

  @override
  Widget build(BuildContext context) {
    return Scaffold(
      appBar: AppBar(
        title: const Text('Material & Laminate'),
        leading: BackButton(
            onPressed: () => context.go('/project/${widget.projectId}/geometry')),
      ),
      body: Row(
        crossAxisAlignment: CrossAxisAlignment.start,
        children: [
          // Left: controls
          SizedBox(
            width: 340,
            child: SingleChildScrollView(
              padding: const EdgeInsets.all(16),
              child: Column(
                crossAxisAlignment: CrossAxisAlignment.start,
                children: [
                  Text('Material', style: Theme.of(context).textTheme.titleMedium),
                  const SizedBox(height: 8),
                  DropdownButtonFormField<OrthotropicPly>(
                    value: _selectedPly,
                    items: MaterialDatabase.compositePlies
                        .map((p) => DropdownMenuItem(value: p, child: Text(p.name)))
                        .toList(),
                    onChanged: (p) => setState(() => _selectedPly = p!),
                    decoration:
                        const InputDecoration(labelText: 'Composite Material'),
                  ),
                  const SizedBox(height: 12),
                  _buildMaterialProperties(),
                  const SizedBox(height: 16),
                  Text('Laminate Stack',
                      style: Theme.of(context).textTheme.titleMedium),
                  const SizedBox(height: 8),
                  Row(
                    children: [
                      Expanded(
                        child: Slider(
                          value: _newAngle,
                          min: -90,
                          max: 90,
                          divisions: 36,
                          label: '${_newAngle.toStringAsFixed(0)}°',
                          onChanged: (v) => setState(() => _newAngle = v),
                        ),
                      ),
                      SizedBox(
                        width: 60,
                        child: Text('${_newAngle.toStringAsFixed(0)}°',
                            textAlign: TextAlign.center),
                      ),
                    ],
                  ),
                  Row(
                    children: [
                      Expanded(
                        child: FilledButton.icon(
                          onPressed: _addPly,
                          icon: const Icon(Icons.add),
                          label: const Text('Add Ply'),
                        ),
                      ),
                      const SizedBox(width: 8),
                      PopupMenuButton<String>(
                        onSelected: _addStandardLaminate,
                        itemBuilder: (_) => const [
                          PopupMenuItem(value: 'ud', child: Text('[0]_8')),
                          PopupMenuItem(value: 'crossply', child: Text('[0/90]_4')),
                          PopupMenuItem(
                              value: 'qi', child: Text('[0/45/-45/90]_s')),
                        ],
                        child: const Chip(label: Text('Presets')),
                      ),
                    ],
                  ),
                  const SizedBox(height: 24),
                  SizedBox(
                    width: double.infinity,
                    child: FilledButton.icon(
                      onPressed: _stack.isNotEmpty
                          ? () => context.go('/project/${widget.projectId}/flight')
                          : null,
                      icon: const Icon(Icons.arrow_forward),
                      label: const Text('Next: Flight Conditions'),
                    ),
                  ),
                ],
              ),
            ),
          ),
          const VerticalDivider(),
          // Right: ply stack + ABD display
          Expanded(
            child: SingleChildScrollView(
              padding: const EdgeInsets.all(16),
              child: Column(
                crossAxisAlignment: CrossAxisAlignment.start,
                children: [
                  Text('Ply Stack (${_stack.length} plies)',
                      style: Theme.of(context).textTheme.titleMedium),
                  const SizedBox(height: 8),
                  if (_stack.isEmpty)
                    const Padding(
                      padding: EdgeInsets.all(32),
                      child: Center(child: Text('No plies added yet')),
                    )
                  else
                    ..._stack.asMap().entries.map((e) => _PlyTile(
                          index: e.key,
                          ply: e.value,
                          onDelete: () => _removePly(e.key),
                        )),
                  if (_abd != null) ...[
                    const SizedBox(height: 16),
                    Text('Laminate Properties',
                        style: Theme.of(context).textTheme.titleMedium),
                    const SizedBox(height: 8),
                    _buildABDDisplay(_abd!),
                  ],
                ],
              ),
            ),
          ),
        ],
      ),
    );
  }

  Widget _buildMaterialProperties() {
    final ply = _selectedPly;
    return Card(
      child: Padding(
        padding: const EdgeInsets.all(12),
        child: Column(
          crossAxisAlignment: CrossAxisAlignment.start,
          children: [
            _propRow('E1', '${(ply.E1 / 1e9).toStringAsFixed(1)} GPa'),
            _propRow('E2', '${(ply.E2 / 1e9).toStringAsFixed(1)} GPa'),
            _propRow('G12', '${(ply.G12 / 1e9).toStringAsFixed(2)} GPa'),
            _propRow('ν12', ply.nu12.toStringAsFixed(3)),
            _propRow('ρ', '${ply.rho.toStringAsFixed(0)} kg/m³'),
            _propRow('t (ply)', '${(ply.t * 1000).toStringAsFixed(3)} mm'),
          ],
        ),
      ),
    );
  }

  Widget _buildABDDisplay(LaminateABD abd) {
    return Column(
      crossAxisAlignment: CrossAxisAlignment.start,
      children: [
        _matrixCard('A matrix (N/m)', abd.A),
        const SizedBox(height: 8),
        _matrixCard('D matrix (N·m)', abd.D),
        const SizedBox(height: 8),
        Card(
          child: Padding(
            padding: const EdgeInsets.all(12),
            child: Column(
              crossAxisAlignment: CrossAxisAlignment.start,
              children: [
                Text('Summary', style: Theme.of(context).textTheme.titleSmall),
                const SizedBox(height: 4),
                _propRow('Total thickness',
                    '${(abd.totalThickness * 1000).toStringAsFixed(2)} mm'),
                _propRow('Areal density',
                    '${abd.areaWeight.toStringAsFixed(2)} kg/m²'),
                _propRow('Bending stiffness D11',
                    '${(abd.D.get(0, 0)).toStringAsExponential(3)} N·m'),
              ],
            ),
          ),
        ),
      ],
    );
  }

  Widget _matrixCard(String title, dynamic matrix) {
    return Card(
      child: Padding(
        padding: const EdgeInsets.all(12),
        child: Column(
          crossAxisAlignment: CrossAxisAlignment.start,
          children: [
            Text(title, style: Theme.of(context).textTheme.titleSmall),
            const SizedBox(height: 8),
            for (var i = 0; i < 3; i++)
              Padding(
                padding: const EdgeInsets.symmetric(vertical: 2),
                child: Row(
                  children: [
                    for (var j = 0; j < 3; j++)
                      Expanded(
                        child: Text(
                          matrix.get(i, j).toStringAsExponential(2),
                          textAlign: TextAlign.right,
                          style: const TextStyle(
                              fontFamily: 'monospace', fontSize: 11),
                        ),
                      ),
                  ],
                ),
              ),
          ],
        ),
      ),
    );
  }

  Widget _propRow(String label, String value) => Padding(
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

class _PlyTile extends StatelessWidget {
  final int index;
  final LaminatePly ply;
  final VoidCallback onDelete;

  const _PlyTile({
    required this.index,
    required this.ply,
    required this.onDelete,
  });

  @override
  Widget build(BuildContext context) {
    return Card(
      margin: const EdgeInsets.only(bottom: 4),
      child: ListTile(
        leading: CircleAvatar(
          radius: 16,
          child: Text('${index + 1}', style: const TextStyle(fontSize: 12)),
        ),
        title: Text(ply.ply.name),
        subtitle: Text('${ply.angleDegrees.toStringAsFixed(0)}° · '
            't=${(ply.ply.t * 1000).toStringAsFixed(3)}mm'),
        trailing: IconButton(
          icon: const Icon(Icons.delete, size: 20),
          onPressed: onDelete,
        ),
        dense: true,
      ),
    );
  }
}
