import 'package:flutter/material.dart';
import 'package:go_router/go_router.dart';
import '../../../models/project.dart';
import '../../../services/project_service.dart';

/// Home screen: project list.
class HomeScreen extends StatefulWidget {
  const HomeScreen({super.key});

  @override
  State<HomeScreen> createState() => _HomeScreenState();
}

class _HomeScreenState extends State<HomeScreen> {
  final _projectService = ProjectService();
  List<Project> _projects = [];
  bool _loading = true;

  @override
  void initState() {
    super.initState();
    _loadProjects();
  }

  Future<void> _loadProjects() async {
    await _projectService.init();
    final projects = await _projectService.loadProjects();
    if (mounted) {
      setState(() {
        _projects = projects;
        _loading = false;
      });
    }
  }

  Future<void> _createProject() async {
    final name = await _showNameDialog(context);
    if (name == null || name.isEmpty) return;
    final project = await _projectService.createProject(name);
    if (mounted) {
      setState(() => _projects.add(project));
      context.go('/project/${project.id}/geometry');
    }
  }

  Future<void> _deleteProject(Project project) async {
    final confirm = await showDialog<bool>(
      context: context,
      builder: (ctx) => AlertDialog(
        title: const Text('Delete Project'),
        content: Text('Delete "${project.name}"? This cannot be undone.'),
        actions: [
          TextButton(onPressed: () => ctx.pop(false), child: const Text('Cancel')),
          FilledButton(
            onPressed: () => ctx.pop(true),
            style: FilledButton.styleFrom(backgroundColor: Colors.red),
            child: const Text('Delete'),
          ),
        ],
      ),
    );
    if (confirm == true) {
      await _projectService.deleteProject(project.id);
      if (mounted) setState(() => _projects.removeWhere((p) => p.id == project.id));
    }
  }

  @override
  Widget build(BuildContext context) {
    return Scaffold(
      appBar: AppBar(
        title: const Text('Fin Flutter Analyzer'),
        actions: [
          IconButton(
            icon: const Icon(Icons.refresh),
            onPressed: _loadProjects,
            tooltip: 'Refresh',
          ),
        ],
      ),
      body: _loading
          ? const Center(child: CircularProgressIndicator())
          : _projects.isEmpty
              ? _buildEmptyState(context)
              : _buildProjectList(),
      floatingActionButton: FloatingActionButton.extended(
        onPressed: _createProject,
        icon: const Icon(Icons.add),
        label: const Text('New Project'),
      ),
    );
  }

  Widget _buildEmptyState(BuildContext context) {
    return Center(
      child: Column(
        mainAxisAlignment: MainAxisAlignment.center,
        children: [
          Icon(Icons.rocket_launch, size: 80, color: Colors.grey[300]),
          const SizedBox(height: 16),
          Text(
            'No Projects Yet',
            style: Theme.of(context).textTheme.headlineSmall,
          ),
          const SizedBox(height: 8),
          Text(
            'Create a new project to start analyzing fin flutter and divergence.',
            style: TextStyle(color: Colors.grey[600]),
            textAlign: TextAlign.center,
          ),
          const SizedBox(height: 24),
          FilledButton.icon(
            onPressed: _createProject,
            icon: const Icon(Icons.add),
            label: const Text('Create First Project'),
          ),
        ],
      ),
    );
  }

  Widget _buildProjectList() {
    return ListView.separated(
      padding: const EdgeInsets.all(16),
      itemCount: _projects.length,
      separatorBuilder: (_, __) => const SizedBox(height: 8),
      itemBuilder: (context, index) {
        final project = _projects[index];
        return _ProjectCard(
          project: project,
          onTap: () => context.go('/project/${project.id}/geometry'),
          onDelete: () => _deleteProject(project),
        );
      },
    );
  }
}

class _ProjectCard extends StatelessWidget {
  final Project project;
  final VoidCallback onTap;
  final VoidCallback onDelete;

  const _ProjectCard({
    required this.project,
    required this.onTap,
    required this.onDelete,
  });

  @override
  Widget build(BuildContext context) {
    return Card(
      child: ListTile(
        leading: const CircleAvatar(child: Icon(Icons.rocket_launch)),
        title: Text(project.name),
        subtitle: Text(
          '${project.sessions.length} sessions · '
          'Updated ${_formatDate(project.updatedAt)}',
        ),
        trailing: PopupMenuButton<String>(
          onSelected: (value) {
            if (value == 'delete') onDelete();
          },
          itemBuilder: (_) => [
            const PopupMenuItem(
              value: 'delete',
              child: Row(
                children: [
                  Icon(Icons.delete, color: Colors.red),
                  SizedBox(width: 8),
                  Text('Delete'),
                ],
              ),
            ),
          ],
        ),
        onTap: onTap,
      ),
    );
  }

  String _formatDate(DateTime dt) {
    final now = DateTime.now();
    final diff = now.difference(dt);
    if (diff.inDays > 0) return '${diff.inDays}d ago';
    if (diff.inHours > 0) return '${diff.inHours}h ago';
    return 'just now';
  }
}

Future<String?> _showNameDialog(BuildContext context) {
  final controller = TextEditingController();
  return showDialog<String>(
    context: context,
    builder: (ctx) => AlertDialog(
      title: const Text('Project Name'),
      content: TextField(
        controller: controller,
        decoration: const InputDecoration(hintText: 'e.g., L2 Rocket Fins'),
        autofocus: true,
      ),
      actions: [
        TextButton(onPressed: () => ctx.pop(null), child: const Text('Cancel')),
        FilledButton(
          onPressed: () => ctx.pop(controller.text.trim()),
          child: const Text('Create'),
        ),
      ],
    ),
  );
}
