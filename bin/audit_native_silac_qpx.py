#!/usr/bin/env python3
"""Audit a retained, single-fraction-group SILAC bundle without reprocessing it.

Requires qpx, pyarrow and pyopenms. Writes only to a NEW --outdir, outside the
input bundle. This is an offline validation increment, not a measurement export.
"""

import argparse
from collections import Counter, defaultdict
from dataclasses import asdict
from functools import lru_cache
import hashlib
import importlib.metadata
import json
import math
from pathlib import Path
import re
import shutil
import statistics
import struct
import subprocess
import sys
import xml.etree.ElementTree as ET

import export_silac_ratios as exporter


def sha256(path):
    digest = hashlib.sha256()
    with path.open('rb') as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()


def cv(row):
    return {item['cv_name']: item['cv_value'] for item in row.get('cv_params') or []}


def f32(value):
    return struct.unpack('f', struct.pack('f', float(value)))[0]


def valid_ratio(numerator, denominator):
    if (numerator is None or denominator is None or not math.isfinite(numerator)
            or not math.isfinite(denominator) or numerator <= 0 or denominator <= 0):
        return None
    ratio = numerator / denominator
    return ratio if math.isfinite(ratio) and ratio > 0 else None


def unique_index(rows, key):
    result = {}
    for row in rows:
        identity = key(row)
        if identity in result:
            raise ValueError(f'Ambiguous or duplicate identity: {identity!r}')
        result[identity] = row
    return result


def reference_errors(psms, features):
    """Check both directions and reciprocity; IDs are local to this bundle."""
    errors = []
    for pid, psm in psms.items():
        fid = psm['feature_id']
        if fid is not None and (fid not in features or pid not in features[fid]['psm_ids']):
            errors.append(('psm_to_feature', pid, fid))
    for fid, feature in features.items():
        for pid in feature['psm_ids']:
            if pid not in psms or psms[pid]['feature_id'] != fid:
                errors.append(('feature_to_psm', fid, pid))
            elif psms[pid]['run_file_name'] != feature['id_run_file_name']:
                errors.append(('identification_run', fid, pid))
    return errors


class Checks:
    def __init__(self):
        self.rows = {}

    def check(self, name, passed, detail=None):
        row = self.rows.setdefault(name, {'checked': 0, 'failed': 0, 'examples': []})
        row['checked'] += 1
        if not passed:
            row['failed'] += 1
            if len(row['examples']) < 5:
                row['examples'].append(detail)

    def number(self, name, actual, expected, *, rel_tol=3e-7, abs_tol=1e-6):
        passed = (actual is not None and expected is not None and
                  math.isfinite(actual) and math.isfinite(expected) and
                  math.isclose(actual, expected, rel_tol=rel_tol, abs_tol=abs_tol))
        safe = lambda value: str(value) if value is not None and not math.isfinite(value) else value
        self.check(name, passed, {'actual': safe(actual), 'expected': safe(expected)})
        row = self.rows[name]
        if actual is not None and expected is not None and math.isfinite(actual) and math.isfinite(expected):
            row['max_absolute_difference'] = max(row.get('max_absolute_difference', 0), abs(actual - expected))
            if expected:
                row['max_relative_difference'] = max(row.get('max_relative_difference', 0), abs((actual - expected) / expected))


def array(params, name, convert=str):
    return exporter.array(params.get(name), convert)


@lru_cache(maxsize=None)
def proforma(sequence):
    from pyopenms import AASequence
    from qpx.converters.openms_consensus.feature_adapter import to_proforma
    return to_proforma(AASequence.fromString(sequence)) if sequence else ''


def read_metadata(sdrf, design, root, checks):
    import csv
    with sdrf.open() as handle:
        physical = list(csv.reader(handle, delimiter='\t'))
    header, records = physical[0], physical[1:]
    # Preserve physical repeated modification columns; do not round-trip DictReader.
    checks.check('sdrf_row_widths', all(len(row) == len(header) for row in records))
    required = ['comment[data file]', 'comment[label]', 'comment[fraction identifier]',
                'characteristics[biological unit]', 'factor value[treatment]']
    indexes = {name: header.index(name) for name in required}
    records = [{name: row[index] for name, index in indexes.items()} for row in records]
    sdrf_by_key = unique_index(records, lambda row: (Path(row['comment[data file]']).stem, row['comment[label]']))
    sections = design.read_text().split('\n\n')
    design_rows = list(csv.DictReader(sections[0].splitlines(), delimiter='\t'))
    unique_index(design_rows, lambda row: (row['Spectra_Filepath'], row['Label']))
    groups, samples = exporter.read_design(design)
    if len(groups) != 1:
        raise ValueError('This first audit supports one fraction group; no cross-group inference is attempted.')
    exporter.check_design_maps(root, groups)
    reference, chemistry, normalized = exporter.read_settings(root)
    # Names come from the recorded design/map chemistry, reconciled with the SDRF.
    maps = {node.get('id'): dict(node.attrib, **exporter.user_params(node))
            for node in root.findall('mapList/map')}
    names = {'': 'SILAC light', 'Lys4,Arg6': 'SILAC medium', 'Lys8,Arg10': 'SILAC heavy'}
    labels = {index + 1: names[value] for index, value in enumerate(chemistry)}
    expected, actual = set(), set()
    manifest = []
    for row in design_rows:
        run, channel = Path(row['Spectra_Filepath']).stem, int(row['Label'])
        key = run, labels[channel]
        source = sdrf_by_key.get(key)
        checks.check('design_sdrf_run_channel', source is not None, key)
        if source is not None:
            checks.check('design_sdrf_fraction', row['Fraction'] == source['comment[fraction identifier]'], key)
            condition = samples[row['Sample']]['MSstats_Condition']
            checks.check('design_sdrf_condition', condition.endswith(source['factor value[treatment]']), key)
            manifest.append({'run': run, 'channel': channel, 'label': labels[channel], 'fraction': int(row['Fraction']),
                             'fraction_group': int(row['Fraction_Group']), 'sample': row['Sample'],
                             'biological_unit': source['characteristics[biological unit]'], 'condition': condition})
        expected.add((run, channel, int(row['Fraction_Group']), int(row['Fraction']), row['Sample']))
    for row in maps.values():
        actual.add((Path(row['name']).stem, int(row['channel_id']) + 1,
                    int(row['fraction_group']), int(row['fraction']), row['sample_name']))
    checks.check('design_xml_map_identity', expected == actual and len(expected) == len(maps))
    checks.check('sdrf_design_exact_population', set(sdrf_by_key) == {(m['run'], m['label']) for m in manifest})
    checks.check('one_biological_unit', len({m['biological_unit'] for m in manifest}) == 1)
    return maps, labels, reference, normalized, manifest


def audit(bundle, sdrf, out):
    import pyarrow as pa
    import pyarrow.parquet as pq
    import qpx

    checks = Checks()
    sources = [bundle / 'experiment.consensusXML', bundle / 'experiment.mzTab', sdrf]
    designs = list(bundle.glob('*openms_design.tsv'))
    if len(designs) != 1:
        raise ValueError('Expected exactly one published experimental design')
    design = designs[0]
    sources.append(design)
    configs = list(bundle.glob('*config.tsv'))
    if len(configs) != 1:
        raise ValueError('Expected exactly one published per-file OpenMS configuration')
    sources.extend(configs)
    sources.extend(bundle / 'qpx' / f'quantms.{kind}.parquet' for kind in ('psm', 'feature', 'pg'))
    inputs = [{'path': str(path.resolve()), 'bytes': path.stat().st_size, 'sha256': sha256(path)} for path in sources]
    root = ET.parse(sources[0]).getroot()
    maps, labels, reference, normalized, manifest = read_metadata(sdrf, design, root, checks)
    if normalized:
        raise ValueError('This audit targets raw pilot ratios; normalized input needs a separate reconciliation.')
    runs = {row['run'] for row in manifest}
    import csv
    with configs[0].open() as handle:
        config_rows = list(csv.DictReader(handle, delimiter='\t'))
    config_by_run = unique_index(config_rows, lambda row: Path(row['Filename']).stem)
    checks.check('configuration_exact_run_population', set(config_by_run) == runs)
    tables, rows, inventory = {}, {}, {}
    for kind in ('psm', 'feature', 'pg'):
        path = bundle / 'qpx' / f'quantms.{kind}.parquet'
        table = pq.read_table(path)
        tables[kind], rows[kind] = table, table.to_pylist()
        inventory[kind] = {'rows': table.num_rows,
                           'fields': [{'name': f.name, 'type': str(f.type), 'nullable': f.nullable,
                                       'null_count': table[f.name].null_count} for f in table.schema],
                           'metadata': {k.decode(): v.decode() for k, v in (table.schema.metadata or {}).items()}}
        (out / f'{kind}.arrow.txt').write_text(table.schema.to_string(truncate_metadata=False) + '\n')
    with qpx.open_dataset(bundle / 'qpx') as dataset:
        validation = dataset.validate(strict=True)
        formal = {kind: {'valid': result.is_valid, 'summary': result.summary,
                         'issues': [asdict(issue) for issue in result.issues]} for kind, result in validation.items()}
    for kind in tables:
        checks.check(f'canonical_schema_{kind}', kind in formal and formal[kind]['valid'])
        unique_index(rows[kind], lambda row, k=kind: row[f'{k}_id'])
    psms = unique_index(rows['psm'], lambda row: row['psm_id'])
    features = unique_index(rows['feature'], lambda row: row['feature_id'])
    for error in reference_errors(psms, features):
        checks.check('bidirectional_psm_feature_references', False, error)
    checks.check('bidirectional_psm_feature_references', not reference_errors(psms, features))

    protein_identification = root.find('.//ProteinIdentification')
    protein_hits = {hit.get('id'): hit for hit in protein_identification.findall('ProteinHit')}
    pmeta = exporter.user_params(protein_identification)
    source_runs = [Path(value).stem for value in array(pmeta, 'spectra_data')]
    xml_psms = {}
    xml_psm_nodes = root.findall('.//PeptideIdentification') + root.findall('UnassignedPeptideIdentification')
    for node in xml_psm_nodes:
        if node.find('PeptideHit') is None:
            continue
        meta = exporter.user_params(node)
        run = source_runs[int(meta['id_merge_index'])]
        scan_match = re.search(r'(?:^|\s)scan=(\d+)(?:\s|$)', node.get('spectrum_reference', ''))
        if scan_match is None:
            raise ValueError('Expected explicit Thermo scan references in retained pilot XML')
        for hit in node.findall('PeptideHit'):
            hmeta = exporter.user_params(hit)
            form = proforma(hmeta.get('MS1Label:labeled_sequence', hit.get('sequence')))
            key = run, (int(scan_match[1]),), form, int(hit.get('charge'))
            if key in xml_psms:
                raise ValueError(f'Duplicate XML PSM key: {key}')
            xml_psms[key] = (node, hit, hmeta)
    qpx_psms = unique_index(rows['psm'], lambda row: (row['run_file_name'], tuple(row['scan']), row['peptidoform'], row['charge']))
    checks.check('psm_xml_exact_identity_population', set(qpx_psms) == set(xml_psms),
                 {'qpx_only': len(set(qpx_psms) - set(xml_psms)), 'xml_only': len(set(xml_psms) - set(qpx_psms))})
    for key in qpx_psms.keys() & xml_psms.keys():
        row, (node, hit, hmeta) = qpx_psms[key], xml_psms[key]
        checks.number('psm_pep_xml', row['posterior_error_probability'], float(hmeta['Posterior Error Probability_score']), rel_tol=1e-12)
        checks.number('psm_rt_xml', row['rt'], float(node.get('RT')))
        checks.number('psm_observed_mz_xml', row['observed_mz'], float(node.get('MZ')))
        checks.check('psm_label_annotation', cv(row).get('MS1Label:channel') == hmeta.get('MS1Label:channel'), key)
        checks.check('psm_target_flag', row['is_decoy'] == (hmeta.get('target_decoy') == 'decoy'), key)
        checks.check('psm_protein_evidence_xml', set(row['protein_accessions']) ==
                     {protein_hits[ref].get('accession') for ref in hit.get('protein_refs', '').split()}, key)
    xml_features = {}
    for node in root.findall('consensusElementList/consensusElement'):
        subelements = node.findall('groupedElementList/element')
        feature_runs = {Path(maps[element.get('map')]['name']).stem for element in subelements}
        if len(feature_runs) != 1:
            raise ValueError('This pilot audit requires each feature to belong to exactly one acquisition.')
        key = next(iter(feature_runs)), int(node.get('charge')), f32(node.find('centroid').get('rt'))
        if key in xml_features:
            raise ValueError(f'Ambiguous XML feature key: {key}')
        xml_features[key] = node
    qpx_features = unique_index(rows['feature'], lambda row: (row['run_file_name'], row['charge'], row['rt']))
    checks.check('feature_xml_exact_identity_population', set(qpx_features) == set(xml_features),
                 {'qpx_only': len(set(qpx_features) - set(xml_features)), 'xml_only': len(set(xml_features) - set(qpx_features))})
    peptide_evidence, native_peptides = defaultdict(list), {}
    for key in qpx_features.keys() & xml_features.keys():
        row, node = qpx_features[key], xml_features[key]
        hit = node.find('PeptideIdentification/PeptideHit')
        checks.check('feature_peptidoform_xml', row['peptidoform'] == proforma(hit.get('sequence') if hit is not None else ''), key)
        checks.check('feature_protein_evidence_xml', {p['accession'] for p in row['pg_accessions']} ==
                     ({protein_hits[ref].get('accession') for ref in hit.get('protein_refs', '').split()}
                      if hit is not None else set()), key)
        checks.number('feature_mz_xml', row['observed_mz'], float(node.find('centroid').get('mz')))
        intensity = unique_index(row['intensities'], lambda item: item['label'])
        checks.check('feature_channels', set(intensity) == set(labels.values()), key)
        for element in node.findall('groupedElementList/element'):
            label = labels[int(maps[element.get('map')]['channel_id']) + 1]
            if label in intensity:
                checks.number('feature_channel_intensity_xml', intensity[label]['intensity'], float(element.get('it')))
        for pid in row['psm_ids']:
            if pid in psms:
                checks.check('feature_psm_run_scan_charge',
                             psms[pid]['run_file_name'] == key[0] and psms[pid]['charge'] == row['charge']
                             and psms[pid]['scan'] == row['scan'], row['feature_id'])
        params = exporter.user_params(node)
        evidence = dict(zip(array(params, 'MS1Label:evidence_ratio_channel', int), array(params, 'MS1Label:evidence_ratio', float)))
        peptide_arrays = dict(zip(array(params, 'MS1Label:peptide_ratio_channel', int),
                                  zip(array(params, 'MS1Label:peptide_ratio', float), array(params, 'MS1Label:peptide_ratio_count', int))))
        denominator = intensity.get(labels[reference], {}).get('intensity')
        for channel, label in labels.items():
            ratio = valid_ratio(intensity.get(label, {}).get('intensity'), denominator)
            checks.check('valid_feature_ratio_denominator', ratio is not None, key)
            if channel in evidence:
                checks.number('feature_evidence_ratio_xml', ratio, evidence[channel])
            if hit is not None and ratio is not None:
                form_key = row['peptidoform'], channel
                peptide_evidence[form_key].append(ratio)
                stored = peptide_arrays.get(channel)
                checks.check('feature_native_peptide_ratio_present', stored is not None, key)
                if stored is not None:
                    checks.check('repeated_peptide_aggregate_consistent', form_key not in native_peptides or native_peptides[form_key] == stored, form_key)
                    native_peptides[form_key] = stored
    for key, values in peptide_evidence.items():
        if key in native_peptides:
            expected, count = native_peptides[key]
            checks.number('peptide_ratio_from_qpx_evidence', statistics.median(values), expected)
            checks.check('peptide_support_xml', len(values) == count, key)

    hits = {hit.get('id'): hit for hit in protein_identification.findall('ProteinHit')}
    native_groups = {}
    for name, value in pmeta.items():
        if re.fullmatch(r'indistinguishable_proteins_\d+', name):
            score, *members = value.split(',')
            identity = tuple(sorted(hits[member].get('accession') for member in members))
            native_groups[identity] = (name, float(score))
    qgroups = {tuple(sorted(row['pg_accessions'])) for row in rows['pg']}
    checks.check('protein_group_exact_membership_population', qgroups == set(native_groups))
    group_sets = {frozenset(group) for group in qgroups}
    mapping_counts = Counter()
    mapping_examples = []
    for row in rows['feature']:
        accessions = frozenset(item['accession'] for item in row['pg_accessions'])
        if not accessions:
            mapping_counts['identification_free'] += 1
        elif accessions in group_sets:
            mapping_counts['exact_inferred_group'] += 1
        else:
            candidates = [sorted(group) for group in group_sets if group <= accessions]
            mapping_counts['multiple_inferred_groups'] += 1
            checks.check('nonexact_feature_group_mapping_is_shared', len(candidates) > 1 and row['unique'] is False, row['feature_id'])
            if len(mapping_examples) < 5:
                mapping_examples.append({'feature_id': row['feature_id'], 'protein_evidence': sorted(accessions),
                                         'compatible_groups': sorted(candidates), 'anchor': row['anchor_protein']})
    pg_keys = unique_index(rows['pg'], lambda row: (tuple(sorted(row['pg_accessions'])), tuple(sorted(row['grouped_runs'])), row['label']))
    for (identity, grouped_runs, label), row in pg_keys.items():
        checks.check('protein_group_run_membership', set(grouped_runs) == runs, row['pg_id'])
        checks.check('protein_group_label', label is None or label in labels.values(), row['pg_id'])
        if identity not in native_groups:
            continue
        prefix, score = native_groups[identity]
        checks.number('protein_group_qvalue_xml', row['global_qvalue'], score, rel_tol=1e-12)
        if label is None:
            checks.check('identification_only_group_no_quantity', row['intensity'] is None and not row['additional_intensities'], row['pg_id'])
            continue
        channel = next(c for c, name in labels.items() if name == label)
        abundance = dict(zip(array(pmeta, prefix + '_fraction_group_level_label', int), array(pmeta, prefix + '_fraction_group_level_abundance', float)))
        checks.number('protein_abundance_xml', row['intensity'], abundance.get(channel))
        ratios = dict(zip(array(pmeta, prefix + '_fraction_group_level_ratio_label', int), array(pmeta, prefix + '_fraction_group_level_ratio', float)))
        counts = dict(zip(array(pmeta, prefix + '_fraction_group_level_ratio_label', int), array(pmeta, prefix + '_fraction_group_level_ratio_count', int)))
        additional = {(item['label'], value['intensity_name']): value['intensity_value']
                      for item in row['additional_intensities'] or [] for value in item['intensities']}
        checks.check('protein_ratio_presence_xml', ((label, 'ratio') in additional) == (channel in ratios), row['pg_id'])
        if channel in ratios:
            checks.number('protein_ratio_xml', additional.get((label, 'ratio')), ratios[channel])
            checks.check('protein_ratio_support_xml', cv(row).get('ratio_count') == str(counts[channel]), row['pg_id'])
    peptide_export, protein_export = exporter.export_ratios(bundle / 'experiment.mzTab', bundle / 'experiment.consensusXML', design)
    for row in peptide_export:
        if row['ratio_status'] != 'reported':
            continue
        key = proforma(row['peptidoform']), row['numerator_channel']
        checks.check('mztab_peptide_ratio_has_qpx_evidence', key in peptide_evidence, key)
        if key in peptide_evidence:
            checks.number('peptide_ratio_qpx_to_mztab', statistics.median(peptide_evidence[key]), row['ratio'])
            checks.check('peptide_support_qpx_to_mztab', len(peptide_evidence[key]) == row['ratio_count'], key)
    counts = {'psms': len(psms), 'features': len(features),
              'identified_features': sum(bool(row['peptidoform']) for row in features.values()),
              'identification_free_features': sum(not row['peptidoform'] for row in features.values()),
              'peptide_forms': len({key[0] for key in peptide_evidence}), 'protein_rows': len(rows['pg']),
              'protein_groups': len(qgroups),
              'peptide_nonreference_ratio_cells': sum(row['ratio_status'] == 'reported' for row in peptide_export),
              'protein_nonreference_ratio_cells': sum(row['ratio_status'] == 'reported' for row in protein_export),
              'protein_label_counts': dict(Counter(str(row['label']) for row in rows['pg'])),
              'psm_runs': dict(sorted(Counter(row['run_file_name'] for row in psms.values()).items())),
              'psms_with_feature_reference': sum(row['feature_id'] is not None for row in psms.values())}
    schema_dir = Path(qpx.__file__).parent / 'core/data/schemas'
    (out / 'canonical_schemas').mkdir()
    for name in ('psm.yaml', 'feature.yaml', 'pg.yaml', 'types.yaml'):
        shutil.copyfile(schema_dir / name, out / 'canonical_schemas' / name)
    report = {'schema_version': 1, 'inputs': inputs, 'counts': counts, 'inventory': inventory,
              'official_validation': formal, 'checks': checks.rows,
              'feature_protein_mapping': {'counts': dict(mapping_counts), 'examples': mapping_examples,
                  'interpretation': 'Native feature.pg_accessions carries peptide protein evidence. For shared peptides it spans multiple inferred groups; never join by anchor or treat the entire list as one inferred group.'},
              'manifest': manifest, 'reference_channel': reference, 'normalization': 'raw',
              'versions': {name: importlib.metadata.version(name) for name in ('qpx', 'pyarrow', 'pyopenms', 'duckdb')},
              'python': sys.version, 'command': sys.argv,
              'git_head': subprocess.check_output(['git', 'rev-parse', 'HEAD'], text=True).strip(),
              'script_sha256': sha256(Path(__file__)),
              'input_hashes_unchanged': all(sha256(Path(item['path'])) == item['sha256'] for item in inputs)}
    checks.check('scientific_inputs_unchanged', report['input_hashes_unchanged'])
    # Do not confuse schema validity with all semantic checks passing.
    report['cross_artifact_failures'] = {name: row for name, row in checks.rows.items() if row['failed']}
    (out / 'validation.json').write_text(json.dumps(report, indent=2, allow_nan=False) + '\n')
    shutil.copyfile(__file__, out / 'audit_native_silac_qpx.py')
    (out / 'summary.md').write_text(
        '# PXD001101 native QPX validation\n\n'
        f"QPX consumer: {report['versions']}. Raw reference channel: {reference}.\n\n"
        + '\n'.join(f'- {kind}: {value["summary"]}' for kind, value in formal.items())
        + f'\n\nCross-artifact failed checks: {len(report["cross_artifact_failures"])}.\n\n'
        + '```json\n' + json.dumps(counts, indent=2) + '\n```\n\n'
        + 'Detailed checks, discrepancies, source hashes, schema inventory and acquisition/channel manifest are in validation.json.\n'
        + '\nThis audit does not certify arbitrary consumers, nested semantic validity, site localization, biological replication, '
        'detection sensitivity, or quantitative accuracy. Peptide aggregate ratios/support remain XML/mzTab-sourced; '
        'QPX intensities are used only for reconciliation. Empty-sequence features remain QC evidence. '
        'Null-label protein rows require explicit handling. PEP is not PSM q-value; group and member scores remain distinct.\n')
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--bundle', required=True, type=Path)
    parser.add_argument('--sdrf', required=True, type=Path)
    parser.add_argument('--outdir', required=True, type=Path)
    args = parser.parse_args()
    out = args.outdir.resolve()
    if out.is_relative_to(args.bundle.resolve()):
        parser.error('--outdir must be outside the original bundle')
    out.mkdir(parents=True, exist_ok=False)
    report = audit(args.bundle.resolve(), args.sdrf.resolve(), out)
    print(json.dumps({'outdir': str(out), 'counts': report['counts'], 'failed_checks': report['cross_artifact_failures']}, indent=2))
    return 1 if report['cross_artifact_failures'] else 0


if __name__ == '__main__':
    sys.exit(main())
