import 'package:flutter/material.dart';
import 'package:flutter/services.dart';

/// A text field for numerical input with unit label.
class NumericalInputField extends StatefulWidget {
  final String label;
  final String? unit;
  final double? value;
  final void Function(double?) onChanged;
  final double? min;
  final double? max;
  final String? hint;
  final int decimalPlaces;

  const NumericalInputField({
    super.key,
    required this.label,
    this.unit,
    this.value,
    required this.onChanged,
    this.min,
    this.max,
    this.hint,
    this.decimalPlaces = 3,
  });

  @override
  State<NumericalInputField> createState() => _NumericalInputFieldState();
}

class _NumericalInputFieldState extends State<NumericalInputField> {
  late TextEditingController _controller;
  String? _errorText;

  @override
  void initState() {
    super.initState();
    _controller = TextEditingController(
      text: widget.value?.toStringAsFixed(widget.decimalPlaces) ?? '',
    );
  }

  @override
  void didUpdateWidget(NumericalInputField oldWidget) {
    super.didUpdateWidget(oldWidget);
    if (oldWidget.value != widget.value && widget.value != null) {
      final newText = widget.value!.toStringAsFixed(widget.decimalPlaces);
      if (_controller.text != newText) {
        _controller.text = newText;
      }
    }
  }

  @override
  void dispose() {
    _controller.dispose();
    super.dispose();
  }

  void _onChanged(String text) {
    if (text.isEmpty) {
      setState(() => _errorText = null);
      widget.onChanged(null);
      return;
    }
    final value = double.tryParse(text);
    if (value == null) {
      setState(() => _errorText = 'Invalid number');
      widget.onChanged(null);
      return;
    }
    if (widget.min != null && value < widget.min!) {
      setState(() => _errorText = 'Min: ${widget.min}');
    } else if (widget.max != null && value > widget.max!) {
      setState(() => _errorText = 'Max: ${widget.max}');
    } else {
      setState(() => _errorText = null);
    }
    widget.onChanged(value);
  }

  @override
  Widget build(BuildContext context) {
    return TextFormField(
      controller: _controller,
      keyboardType: const TextInputType.numberWithOptions(decimal: true, signed: true),
      inputFormatters: [
        FilteringTextInputFormatter.allow(RegExp(r'[-+]?[0-9]*\.?[0-9]*')),
      ],
      onChanged: _onChanged,
      decoration: InputDecoration(
        labelText: widget.label,
        hintText: widget.hint,
        errorText: _errorText,
        suffixText: widget.unit,
      ),
    );
  }
}
