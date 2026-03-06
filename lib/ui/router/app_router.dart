import 'package:flutter/material.dart';
import 'package:go_router/go_router.dart';
import '../screens/home/home_screen.dart';
import '../screens/geometry/geometry_screen.dart';
import '../screens/materials/material_screen.dart';
import '../screens/flight_conditions/flight_condition_screen.dart';
import '../screens/analysis/analysis_screen.dart';
import '../screens/results/results_screen.dart';
import '../screens/optimization/optimization_screen.dart';
import '../screens/project/project_screen.dart';

/// Application router using go_router.
class AppRouter {
  static final GoRouter router = GoRouter(
    initialLocation: '/',
    routes: [
      GoRoute(
        path: '/',
        builder: (context, state) => const HomeScreen(),
      ),
      GoRoute(
        path: '/project/:id',
        builder: (context, state) {
          final id = state.pathParameters['id']!;
          return ProjectScreen(projectId: id);
        },
        routes: [
          GoRoute(
            path: 'geometry',
            builder: (context, state) {
              final id = state.pathParameters['id']!;
              return GeometryScreen(projectId: id);
            },
          ),
          GoRoute(
            path: 'materials',
            builder: (context, state) {
              final id = state.pathParameters['id']!;
              return MaterialScreen(projectId: id);
            },
          ),
          GoRoute(
            path: 'flight',
            builder: (context, state) {
              final id = state.pathParameters['id']!;
              return FlightConditionScreen(projectId: id);
            },
          ),
          GoRoute(
            path: 'analysis',
            builder: (context, state) {
              final id = state.pathParameters['id']!;
              return AnalysisScreen(projectId: id);
            },
          ),
          GoRoute(
            path: 'results',
            builder: (context, state) {
              final id = state.pathParameters['id']!;
              return ResultsScreen(projectId: id);
            },
          ),
          GoRoute(
            path: 'optimization',
            builder: (context, state) {
              final id = state.pathParameters['id']!;
              return OptimizationScreen(projectId: id);
            },
          ),
        ],
      ),
    ],
  );
}
