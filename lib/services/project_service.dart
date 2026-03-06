import 'dart:convert';
import 'package:shared_preferences/shared_preferences.dart';
import '../models/project.dart';
import '../core/models/fin_geometry.dart';
import '../core/models/flight_condition.dart';
import '../core/materials/laminate.dart';
import '../core/materials/clt_calculator.dart';
import '../core/materials/orthotropic_ply.dart';
import '../core/materials/material_database.dart';
import '../models/analysis_session.dart';

/// CRUD service for projects with JSON persistence.
class ProjectService {
  static const String _storageKey = 'fin_flutter_projects';
  SharedPreferences? _prefs;

  Future<void> init() async {
    _prefs = await SharedPreferences.getInstance();
  }

  /// Load all projects from storage.
  Future<List<Project>> loadProjects() async {
    await _ensureInit();
    final jsonStr = _prefs!.getString(_storageKey);
    if (jsonStr == null) return [];
    try {
      final List<dynamic> data = jsonDecode(jsonStr) as List<dynamic>;
      return data.map((d) => _projectFromJson(d as Map<String, dynamic>)).toList();
    } catch (e) {
      return [];
    }
  }

  /// Save all projects to storage.
  Future<void> saveProjects(List<Project> projects) async {
    await _ensureInit();
    final data = projects.map(_projectToJson).toList();
    await _prefs!.setString(_storageKey, jsonEncode(data));
  }

  /// Create a new project.
  Future<Project> createProject(String name, {String description = ''}) async {
    final project = Project(
      id: DateTime.now().millisecondsSinceEpoch.toString(),
      name: name,
      description: description,
      sessions: [],
      createdAt: DateTime.now(),
      updatedAt: DateTime.now(),
    );
    final projects = await loadProjects();
    projects.add(project);
    await saveProjects(projects);
    return project;
  }

  /// Update an existing project.
  Future<void> updateProject(Project project) async {
    final projects = await loadProjects();
    final idx = projects.indexWhere((p) => p.id == project.id);
    if (idx >= 0) {
      projects[idx] = project;
    } else {
      projects.add(project);
    }
    await saveProjects(projects);
  }

  /// Delete a project by ID.
  Future<void> deleteProject(String id) async {
    final projects = await loadProjects();
    projects.removeWhere((p) => p.id == id);
    await saveProjects(projects);
  }

  Future<void> _ensureInit() async {
    _prefs ??= await SharedPreferences.getInstance();
  }

  // --- JSON serialization ---

  Map<String, dynamic> _projectToJson(Project p) => {
        'id': p.id,
        'name': p.name,
        'description': p.description,
        'sessions': p.sessions.map(_sessionToJson).toList(),
        'createdAt': p.createdAt.toIso8601String(),
        'updatedAt': p.updatedAt.toIso8601String(),
      };

  Project _projectFromJson(Map<String, dynamic> json) => Project(
        id: json['id'] as String,
        name: json['name'] as String,
        description: json['description'] as String? ?? '',
        sessions: (json['sessions'] as List<dynamic>?)
                ?.map((s) => _sessionFromJson(s as Map<String, dynamic>))
                .toList() ??
            [],
        createdAt: DateTime.parse(json['createdAt'] as String),
        updatedAt: DateTime.parse(json['updatedAt'] as String),
      );

  Map<String, dynamic> _sessionToJson(AnalysisSession s) => {
        'id': s.id,
        'name': s.name,
        'geometry': s.geometry.toJson(),
        'flightCondition': s.flightCondition.toJson(),
        'laminateName': s.laminate.name,
        'plyCount': s.laminate.plyCount,
        'createdAt': s.createdAt.toIso8601String(),
        'updatedAt': s.updatedAt.toIso8601String(),
      };

  AnalysisSession _sessionFromJson(Map<String, dynamic> json) {
    // Simplified reconstruction - real app would serialize full laminate
    final ply = MaterialDatabase.as43501;
    final laminate = Laminate.unidirectional(ply, json['plyCount'] as int? ?? 4);

    return AnalysisSession(
      id: json['id'] as String,
      name: json['name'] as String,
      geometry: FinGeometry.fromJson(json['geometry'] as Map<String, dynamic>),
      laminate: laminate,
      flightCondition:
          FlightCondition.fromJson(json['flightCondition'] as Map<String, dynamic>),
      results: null, // Results are not persisted (re-computed on demand)
      createdAt: DateTime.parse(json['createdAt'] as String),
      updatedAt: DateTime.parse(json['updatedAt'] as String),
    );
  }
}
