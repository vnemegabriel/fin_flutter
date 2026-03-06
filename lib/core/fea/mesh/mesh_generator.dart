import '../../models/fin_geometry.dart';
import 'mesh.dart';

/// Generates structured quadrilateral meshes over a trapezoidal fin planform.
class MeshGenerator {
  /// Generate a structured quad mesh.
  ///
  /// [geometry] - fin geometry (span, chords, sweep)
  /// [spanwiseElements] - ny: number of elements in the spanwise direction
  /// [chordwiseElements] - nx: number of elements in the chordwise direction
  /// [elementType] - element formulation
  Mesh generateFinMesh({
    required FinGeometry geometry,
    required int spanwiseElements,
    required int chordwiseElements,
    required ElementType elementType,
  }) {
    final ny = spanwiseElements;
    final nx = chordwiseElements;

    // Create nodes: (ny+1) x (nx+1) grid
    final nodes = <Node>[];
    final nodeIndex = (int i, int j) => i * (nx + 1) + j;

    for (var i = 0; i <= ny; i++) {
      final s = i / ny; // normalized spanwise coordinate [0,1]
      final y = s * geometry.span;
      final xLE = geometry.leadingEdgeX(y);
      final chord = geometry.chordAt(y);

      for (var j = 0; j <= nx; j++) {
        final t = j / nx; // normalized chordwise coordinate [0,1]
        final x = xLE + t * chord;
        nodes.add(Node(id: nodeIndex(i, j), x: x, y: y));
      }
    }

    // Create elements
    final elements = <QuadElement>[];
    var elemId = 0;
    for (var i = 0; i < ny; i++) {
      for (var j = 0; j < nx; j++) {
        // Counter-clockwise ordering:
        // n0(i,j) -> n1(i,j+1) -> n2(i+1,j+1) -> n3(i+1,j)
        elements.add(QuadElement(
          id: elemId++,
          nodeIds: [
            nodeIndex(i, j),
            nodeIndex(i, j + 1),
            nodeIndex(i + 1, j + 1),
            nodeIndex(i + 1, j),
          ],
        ));
      }
    }

    return Mesh(
      nodes: nodes,
      elements: elements,
      elementType: elementType,
    );
  }
}
