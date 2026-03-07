import 'package:flutter/material.dart';
import 'package:go_router/go_router.dart';

/// Project shell screen with navigation rail.
class ProjectScreen extends StatefulWidget {
  final String projectId;

  const ProjectScreen({super.key, required this.projectId});

  @override
  State<ProjectScreen> createState() => _ProjectScreenState();
}

class _ProjectScreenState extends State<ProjectScreen> {
  int _selectedIndex = 0;

  static const List<_NavItem> _navItems = [
    _NavItem(icon: Icons.straighten, label: 'Geometry', path: 'geometry'),
    _NavItem(icon: Icons.layers, label: 'Materials', path: 'materials'),
    _NavItem(icon: Icons.air, label: 'Flight', path: 'flight'),
    _NavItem(icon: Icons.play_circle, label: 'Analysis', path: 'analysis'),
    _NavItem(icon: Icons.bar_chart, label: 'Results', path: 'results'),
    _NavItem(icon: Icons.tune, label: 'Optimize', path: 'optimization'),
  ];

  void _navigateTo(int index) {
    setState(() => _selectedIndex = index);
    context.go('/project/${widget.projectId}/${_navItems[index].path}');
  }

  @override
  Widget build(BuildContext context) {
    final isWide = MediaQuery.sizeOf(context).width > 600;

    return Scaffold(
      appBar: AppBar(
        title: const Text('Fin Flutter Analysis'),
        leading: BackButton(onPressed: () => context.go('/')),
      ),
      body: Row(
        children: [
          if (isWide)
            NavigationRail(
              selectedIndex: _selectedIndex,
              onDestinationSelected: _navigateTo,
              labelType: NavigationRailLabelType.all,
              destinations: _navItems
                  .map((item) => NavigationRailDestination(
                        icon: Icon(item.icon),
                        label: Text(item.label),
                      ))
                  .toList(),
            ),
          const VerticalDivider(thickness: 1, width: 1),
          Expanded(
            child: Center(
              child: Text(
                'Select a section from the navigation',
                style: Theme.of(context).textTheme.bodyLarge,
              ),
            ),
          ),
        ],
      ),
      bottomNavigationBar: isWide
          ? null
          : NavigationBar(
              selectedIndex: _selectedIndex,
              onDestinationSelected: _navigateTo,
              destinations: _navItems
                  .map((item) => NavigationDestination(
                        icon: Icon(item.icon),
                        label: item.label,
                      ))
                  .toList(),
            ),
    );
  }
}

class _NavItem {
  final IconData icon;
  final String label;
  final String path;
  const _NavItem({required this.icon, required this.label, required this.path});
}
