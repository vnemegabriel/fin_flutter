/// A node in the FEA mesh.
class Node {
  final int id;
  final double x;  // physical x-coordinate (m)
  final double y;  // physical y-coordinate (m)
  final double z;  // physical z-coordinate (m), 0 for 2D plate

  const Node({required this.id, required this.x, required this.y, this.z = 0.0});

  @override
  String toString() => 'Node($id: [$x, $y, $z])';
}

/// Element type enumeration.
enum ElementType { kirchhoffDKQ, mindlinMITC4 }

/// A quadrilateral FEA plate element with 4 corner nodes.
class QuadElement {
  final int id;
  final List<int> nodeIds; // 4 corner node IDs, counter-clockwise

  const QuadElement({required this.id, required this.nodeIds})
      : assert(nodeIds.length == 4);

  @override
  String toString() => 'QuadElement($id: $nodeIds)';
}

/// FEA mesh: collection of nodes and elements.
class Mesh {
  final List<Node> nodes;
  final List<QuadElement> elements;
  final ElementType elementType;

  const Mesh({
    required this.nodes,
    required this.elements,
    required this.elementType,
  });

  int get nodeCount => nodes.length;
  int get elementCount => elements.length;

  /// DOF per node for the element type.
  int get dofPerNode => switch (elementType) {
        ElementType.kirchhoffDKQ => 3,  // w, theta_x, theta_y
        ElementType.mindlinMITC4 => 5,  // w, theta_x, theta_y, u, v
      };

  /// Total DOF count.
  int get totalDof => nodeCount * dofPerNode;

  /// Get node by ID.
  Node nodeById(int id) => nodes.firstWhere((n) => n.id == id);

  @override
  String toString() =>
      'Mesh(${nodeCount} nodes, ${elementCount} elements, $elementType)';
}
