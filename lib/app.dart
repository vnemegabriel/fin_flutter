import 'package:flutter/material.dart';
import 'ui/router/app_router.dart';
import 'ui/theme/app_theme.dart';

/// Root application widget.
class FinFlutterApp extends StatelessWidget {
  const FinFlutterApp({super.key});

  @override
  Widget build(BuildContext context) {
    return MaterialApp.router(
      title: 'Fin Flutter Analyzer',
      theme: AppTheme.light,
      darkTheme: AppTheme.dark,
      routerConfig: AppRouter.router,
      debugShowCheckedModeBanner: false,
    );
  }
}
