import 'dart:typed_data';
import 'package:archive/archive.dart';
import 'package:xml/xml.dart';
import 'ork_geometry.dart';

/// Parser for OpenRocket .ork files.
///
/// .ork files are ZIP archives containing an XML document.
/// This parser extracts fin geometry and simulation data.
class OrkParser {
  /// Parse an .ork file from its raw bytes.
  Future<OrkDocument> parse(Uint8List fileBytes) async {
    // Step 1: Extract ZIP archive
    Archive archive;
    try {
      archive = ZipDecoder().decodeBytes(fileBytes);
    } catch (e) {
      throw OrkParseException('Failed to decode .ork ZIP archive: $e');
    }

    // Step 2: Find the XML document (usually named 'rocket.ork' or similar)
    ArchiveFile? xmlFile;
    for (final file in archive) {
      if (file.isFile && (file.name.endsWith('.ork') || file.name == 'rocket.ork')) {
        xmlFile = file;
        break;
      }
      // Some .ork files just have the XML at root level in the ZIP
      if (file.isFile && file.name.toLowerCase().contains('xml')) {
        xmlFile = file;
        break;
      }
    }

    // If no file found, try the first file
    xmlFile ??= archive.files.where((f) => f.isFile).firstOrNull;

    if (xmlFile == null) {
      throw OrkParseException('No XML document found in .ork archive');
    }

    // Step 3: Parse XML
    final xmlBytes = xmlFile.content as List<int>;
    final xmlString = String.fromCharCodes(xmlBytes);
    XmlDocument xmlDoc;
    try {
      xmlDoc = XmlDocument.parse(xmlString);
    } catch (e) {
      throw OrkParseException('Failed to parse XML in .ork file: $e');
    }

    // Step 4: Extract data
    return _extractDocument(xmlDoc);
  }

  OrkDocument _extractDocument(XmlDocument doc) {
    final root = doc.rootElement;

    // Rocket name
    final rocketName = _findText(root, 'name') ?? 'Unknown Rocket';

    // Find fin sets
    final finSets = <OrkFinGeometry>[];
    for (final finNode in root.findAllElements('trapezoidfinset')
        .followedBy(root.findAllElements('freeformfinset'))
        .followedBy(root.findAllElements('ellipticalfinset'))) {
      final geo = _extractFinGeometry(finNode);
      if (geo != null) finSets.add(geo);
    }

    // Find simulation data
    final flightProfiles = <OrkFlightProfile>[];
    for (final simNode in root.findAllElements('simulation')) {
      final profile = _extractFlightProfile(simNode);
      if (profile != null && profile.dataPoints.isNotEmpty) {
        flightProfiles.add(profile);
      }
    }

    return OrkDocument(
      rocketName: rocketName,
      finSets: finSets,
      flightProfiles: flightProfiles,
    );
  }

  OrkFinGeometry? _extractFinGeometry(XmlElement finNode) {
    try {
      // Parse geometry fields - OpenRocket stores in SI units (meters)
      final span = _findDouble(finNode, 'span') ??
          _findDouble(finNode, 'height') ?? 0.1;
      final rootChord = _findDouble(finNode, 'rootchord') ??
          _findDouble(finNode, 'length') ?? 0.1;
      final tipChord = _findDouble(finNode, 'tipchord') ??
          _findDouble(finNode, 'rootchord') ?? rootChord * 0.5;
      final sweepLength = _findDouble(finNode, 'sweeplength') ??
          _findDouble(finNode, 'sweep') ?? 0.0;
      final thickness = _findDouble(finNode, 'thickness') ?? 0.003;
      final crossSection = _findText(finNode, 'crosssection') ?? 'square';

      final finCountStr = _findText(finNode, 'fincount');
      final finCount = int.tryParse(finCountStr ?? '4') ?? 4;

      // Material
      final materialNode = finNode.findElements('material').firstOrNull;
      final materialName = materialNode?.getAttribute('name') ??
          _findText(finNode, 'material') ?? 'Unknown';

      // Mount diameter (look at parent rocket body tube)
      final bodyTubeNode = finNode.parent?.parent;
      final mountDiam = _findDouble(bodyTubeNode as XmlElement?, 'outerradius');

      return OrkFinGeometry(
        span: span,
        rootChord: rootChord,
        tipChord: tipChord,
        sweepLength: sweepLength,
        thickness: thickness,
        crossSection: crossSection,
        materialName: materialName,
        finCount: finCount,
        mountDiameter: (mountDiam ?? 0.05) * 2.0,
      );
    } catch (e) {
      return null;
    }
  }

  OrkFlightProfile? _extractFlightProfile(XmlElement simNode) {
    final dataPoints = <OrkFlightDataPoint>[];

    try {
      final flightDataNode = simNode.findElements('flightdata').firstOrNull;
      if (flightDataNode == null) return null;

      for (final dp in flightDataNode.findElements('datapoint')) {
        final time = _findDouble(dp, 't') ?? _findDouble(dp, 'time') ?? 0;
        final altitude = _findDouble(dp, 'altitude') ?? _findDouble(dp, 'z') ?? 0;
        final velocity = _findDouble(dp, 'velocity') ??
            _findDouble(dp, 'totalvelocity') ?? 0;
        final mach = _findDouble(dp, 'mach') ??
            (velocity / 343.0); // estimate if not provided

        dataPoints.add(OrkFlightDataPoint(
          time: time,
          altitude: altitude,
          velocity: velocity,
          mach: mach,
        ));
      }
    } catch (e) {
      // Continue with what we have
    }

    return OrkFlightProfile(dataPoints: dataPoints);
  }

  String? _findText(XmlElement? parent, String tagName) {
    if (parent == null) return null;
    try {
      return parent.findElements(tagName).firstOrNull?.innerText.trim();
    } catch (_) {
      return null;
    }
  }

  double? _findDouble(XmlElement? parent, String tagName) {
    final text = _findText(parent, tagName);
    if (text == null) return null;
    return double.tryParse(text);
  }
}

/// Exception thrown when .ork parsing fails.
class OrkParseException implements Exception {
  final String message;
  const OrkParseException(this.message);

  @override
  String toString() => 'OrkParseException: $message';
}
