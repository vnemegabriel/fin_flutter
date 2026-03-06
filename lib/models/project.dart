import 'analysis_session.dart';

/// A project containing one or more analysis sessions.
class Project {
  final String id;
  final String name;
  final String description;
  final List<AnalysisSession> sessions;
  final DateTime createdAt;
  final DateTime updatedAt;

  const Project({
    required this.id,
    required this.name,
    this.description = '',
    this.sessions = const [],
    required this.createdAt,
    required this.updatedAt,
  });

  Project copyWith({
    String? name,
    String? description,
    List<AnalysisSession>? sessions,
  }) =>
      Project(
        id: id,
        name: name ?? this.name,
        description: description ?? this.description,
        sessions: sessions ?? this.sessions,
        createdAt: createdAt,
        updatedAt: DateTime.now(),
      );

  AnalysisSession? get activeSession =>
      sessions.isNotEmpty ? sessions.last : null;

  @override
  String toString() => 'Project($id: $name, ${sessions.length} sessions)';
}
