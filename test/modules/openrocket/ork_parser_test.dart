import 'dart:typed_data';
import 'package:archive/archive.dart';
import 'package:flutter_test/flutter_test.dart';
import 'package:fin_flutter/modules/openrocket/ork_parser.dart';

/// Build a ZIP archive in memory containing a single XML file named 'rocket.ork'.
Uint8List _buildOrkZip(String xmlContent) {
  final archive = Archive();
  final bytes = xmlContent.codeUnits;
  final file = ArchiveFile('rocket.ork', bytes.length, bytes);
  archive.addFile(file);
  final zipBytes = ZipEncoder().encode(archive)!;
  return Uint8List.fromList(zipBytes);
}

const _validOrkXml = '''<?xml version="1.0" encoding="UTF-8"?>
<openrocket version="1.7">
  <rocket>
    <name>TestRocket</name>
    <subcomponents>
      <stage>
        <subcomponents>
          <bodytube>
            <outerradius>0.025</outerradius>
            <subcomponents>
              <trapezoidfinset>
                <name>TestFin</name>
                <rootchord>0.15</rootchord>
                <tipchord>0.08</tipchord>
                <span>0.20</span>
                <sweeplength>0.05</sweeplength>
                <thickness>0.003</thickness>
                <fincount>4</fincount>
              </trapezoidfinset>
            </subcomponents>
          </bodytube>
        </subcomponents>
      </stage>
    </subcomponents>
  </rocket>
</openrocket>''';

const _missingTipChordXml = '''<?xml version="1.0" encoding="UTF-8"?>
<openrocket version="1.7">
  <rocket>
    <subcomponents>
      <stage>
        <subcomponents>
          <bodytube>
            <subcomponents>
              <trapezoidfinset>
                <rootchord>0.15</rootchord>
                <span>0.20</span>
              </trapezoidfinset>
            </subcomponents>
          </bodytube>
        </subcomponents>
      </stage>
    </subcomponents>
  </rocket>
</openrocket>''';

const _withFlightDataXml = '''<?xml version="1.0" encoding="UTF-8"?>
<openrocket version="1.7">
  <rocket>
    <subcomponents>
      <stage>
        <subcomponents>
          <bodytube>
            <subcomponents>
              <trapezoidfinset>
                <rootchord>0.1</rootchord>
                <tipchord>0.06</tipchord>
                <span>0.15</span>
                <sweeplength>0.03</sweeplength>
                <thickness>0.003</thickness>
              </trapezoidfinset>
            </subcomponents>
          </bodytube>
        </subcomponents>
      </stage>
    </subcomponents>
    <simulations>
      <simulation>
        <flightdata>
          <datapoint><t>0.0</t><altitude>0.0</altitude><velocity>0.0</velocity></datapoint>
          <datapoint><t>1.0</t><altitude>50.0</altitude><velocity>100.0</velocity></datapoint>
          <datapoint><t>2.0</t><altitude>200.0</altitude><velocity>250.0</velocity></datapoint>
        </flightdata>
      </simulation>
    </simulations>
  </rocket>
</openrocket>''';

const _emptyFinSetXml = '''<?xml version="1.0" encoding="UTF-8"?>
<openrocket version="1.7">
  <rocket>
    <subcomponents>
      <stage>
        <subcomponents>
          <bodytube/>
        </subcomponents>
      </stage>
    </subcomponents>
  </rocket>
</openrocket>''';

void main() {
  group('OrkParser', () {
    final parser = OrkParser();

    test('parses valid synthetic .ork file — correct geometry values', () async {
      final bytes = _buildOrkZip(_validOrkXml);
      final doc = await parser.parse(bytes);

      expect(doc.hasFinSets, isTrue);
      final fin = doc.primaryFinSet!;
      expect(fin.span, closeTo(0.20, 0.001));
      expect(fin.rootChord, closeTo(0.15, 0.001));
      expect(fin.tipChord, closeTo(0.08, 0.001));
      expect(fin.sweepLength, closeTo(0.05, 0.001));
    });

    test('missing tipchord falls back to rootChord (no crash)', () async {
      final bytes = _buildOrkZip(_missingTipChordXml);
      final doc = await parser.parse(bytes);

      expect(doc.hasFinSets, isTrue);
      final fin = doc.primaryFinSet!;
      // Parser falls back to rootchord when tipchord is missing
      expect(fin.tipChord, greaterThan(0.0));
      expect(fin.rootChord, closeTo(0.15, 0.001));
    });

    test('invalid ZIP bytes throw OrkParseException', () async {
      final garbage = Uint8List.fromList([0, 1, 2, 3, 4, 5, 6, 7, 8, 9]);
      expect(
        () async => await parser.parse(garbage),
        throwsA(isA<OrkParseException>()),
      );
    });

    test('empty fin set returns empty finSets list (no crash)', () async {
      final bytes = _buildOrkZip(_emptyFinSetXml);
      final doc = await parser.parse(bytes);

      expect(doc.finSets, isEmpty);
      expect(doc.hasFinSets, isFalse);
      expect(doc.primaryFinSet, isNull);
    });

    test('OrkFinGeometry.toFinGeometry() produces valid FinGeometry', () async {
      final bytes = _buildOrkZip(_validOrkXml);
      final doc = await parser.parse(bytes);

      final orkFin = doc.primaryFinSet!;
      final geo = orkFin.toFinGeometry();

      expect(geo.span, closeTo(0.20, 0.001));
      expect(geo.rootChord, closeTo(0.15, 0.001));
      expect(geo.tipChord, closeTo(0.08, 0.001));
      expect(geo.planformArea, greaterThan(0.0));
      expect(geo.aspectRatio, greaterThan(0.0));
    });
  });
}
