import 'package:flutter/material.dart';
import 'package:go_router/go_router.dart';
import '../../../core/models/flight_condition.dart';
import '../../widgets/numerical_input_field.dart';

/// Flight condition input screen.
class FlightConditionScreen extends StatefulWidget {
  final String projectId;

  const FlightConditionScreen({super.key, required this.projectId});

  @override
  State<FlightConditionScreen> createState() => _FlightConditionScreenState();
}

class _FlightConditionScreenState extends State<FlightConditionScreen> {
  double _altitude = 1000.0;
  double _mach = 0.5;
  bool _useISA = true;

  FlightCondition get _condition =>
      FlightCondition.isa(altitude: _altitude, mach: _mach);

  @override
  Widget build(BuildContext context) {
    final cond = _condition;
    return Scaffold(
      appBar: AppBar(
        title: const Text('Flight Conditions'),
        leading: BackButton(
            onPressed: () => context.go('/project/${widget.projectId}/materials')),
      ),
      body: SingleChildScrollView(
        padding: const EdgeInsets.all(24),
        child: ConstrainedBox(
          constraints: const BoxConstraints(maxWidth: 600),
          child: Column(
            crossAxisAlignment: CrossAxisAlignment.start,
            children: [
              Text('Atmospheric Conditions',
                  style: Theme.of(context).textTheme.titleLarge),
              const SizedBox(height: 16),
              SwitchListTile(
                title: const Text('Use ISA Standard Atmosphere'),
                subtitle:
                    const Text('Automatically compute density from altitude'),
                value: _useISA,
                onChanged: (v) => setState(() => _useISA = v),
              ),
              const SizedBox(height: 16),
              NumericalInputField(
                label: 'Altitude',
                unit: 'm',
                value: _altitude,
                min: 0,
                max: 50000,
                onChanged: (v) => setState(() => _altitude = v ?? _altitude),
              ),
              const SizedBox(height: 12),
              NumericalInputField(
                label: 'Mach Number',
                unit: 'M',
                value: _mach,
                min: 0.01,
                max: 0.85,
                hint: 'Must be < 0.85 for P-G correction',
                onChanged: (v) => setState(() => _mach = v ?? _mach),
              ),
              const SizedBox(height: 24),
              // Derived conditions
              Card(
                child: Padding(
                  padding: const EdgeInsets.all(16),
                  child: Column(
                    crossAxisAlignment: CrossAxisAlignment.start,
                    children: [
                      Text('ISA Computed Conditions',
                          style: Theme.of(context).textTheme.titleMedium),
                      const SizedBox(height: 12),
                      _row('Velocity', '${cond.velocity.toStringAsFixed(1)} m/s'),
                      _row('Air Density', '${cond.density.toStringAsFixed(4)} kg/m³'),
                      _row('Dynamic Pressure q',
                          '${cond.dynamicPressure.toStringAsFixed(0)} Pa'),
                      _row('Speed of Sound',
                          '${cond.speedOfSound.toStringAsFixed(1)} m/s'),
                    ],
                  ),
                ),
              ),
              if (_mach >= 0.85)
                const Padding(
                  padding: EdgeInsets.only(top: 8),
                  child: Card(
                    color: Colors.amber,
                    child: Padding(
                      padding: EdgeInsets.all(12),
                      child: Row(
                        children: [
                          Icon(Icons.warning),
                          SizedBox(width: 8),
                          Expanded(
                            child: Text(
                              'Mach ≥ 0.85: Prandtl-Glauert correction not valid. '
                              'Results may be inaccurate in transonic regime.',
                            ),
                          ),
                        ],
                      ),
                    ),
                  ),
                ),
              const SizedBox(height: 32),
              SizedBox(
                width: double.infinity,
                child: FilledButton.icon(
                  onPressed: () =>
                      context.go('/project/${widget.projectId}/analysis'),
                  icon: const Icon(Icons.arrow_forward),
                  label: const Text('Next: Analysis Setup'),
                ),
              ),
            ],
          ),
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
