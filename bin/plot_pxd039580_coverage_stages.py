#!/usr/bin/env python3
"""Plot the retained six-injection HeLa baselines and two-isotope trial."""
import argparse
import collections
import csv
import json
from pathlib import Path
import re
import shutil
import xml.etree.ElementTree as ET

from plot_pxd003327_coverage_stages import extract_stages, write_tsv
from compare_pxd003327_runs import ROOT, sha, plt
from matplotlib.ticker import StrMethodFormatter


ANALYSES = {
    'original': 'results/pxd039580_raw_all6_20261006_210459/analysis',
    'upgrade': 'results/openms37_validation/hela_review/analysis',
    'isotopes2': 'results/openms37_tuning/hela_detector_isotopes2/review/analysis',
}
LOGS = {
    'upgrade': 'results/openms37_validation/hela_review/ms1.command.out',
    'isotopes2': 'work/69/eb8b7a4ae49cc82267e806579a1247/.command.out',
}
LABELS = ['3.6.0\nOriginal', '3.7.0\nBaseline', '3.7.0\n2-isotope trial']
COLORS = ['#7b8797', '#d09236', '#137f88']
METRICS = [
    ('identified_target_psms', 'Identified', 'Target PSMs in the final identification output'),
    ('mapped_multiplets', 'Mapped', 'MS1 multiplets with mapped identifications'),
    ('retained_multiplets', 'Retained', 'Multiplets kept after label-consistency checks'),
    ('quantified_peptide_forms', 'Quantified', 'Distinct modified peptides with usable native ratios'),
]


def table(path):
    with path.open() as handle:
        yield from csv.DictReader(handle, delimiter='\t')


def count_psms(path):
    names, counts = {}, collections.Counter()
    with path.open() as handle:
        for fields in csv.reader(handle, delimiter='\t'):
            if fields and fields[0] == 'MTD' and re.fullmatch(r'ms_run\[\d+\]-location', fields[1]):
                names[re.search(r'\d+', fields[1])[0]] = Path(fields[2]).stem
            elif fields and fields[0] == 'PSH':
                spectrum = fields.index('spectra_ref')
                decoy = fields.index('opt_global_cv_MS:1002217_decoy_peptide')
            elif fields and fields[0] == 'PSM':
                assert fields[decoy] == '0'
                counts[names[re.search(r'ms_run\[(\d+)\]', fields[spectrum])[1]]] += 1
    return counts


def count_exported(path, parameters):
    """Recount run/multiplet pairs while bounding memory for large consensusXML."""
    maps, counts = {}, collections.Counter()
    iterator = ET.iterparse(path, events=('start', 'end'))
    _, root = next(iterator)
    verified_settings = False
    for event, node in iterator:
        if event != 'end':
            continue
        if node.tag == 'dataProcessing':
            software = node.find('software')
            if software is not None and software.get('name') == 'MS1LabeledWorkflow':
                actual = {p.get('name').removeprefix('parameter: '): p.get('value')
                          for p in node.findall('UserParam') if p.get('name').startswith('parameter: ')}
                assert actual == parameters
                verified_settings = True
            node.clear()
        elif node.tag == 'map':
            maps[node.get('id')] = Path(node.get('name')).stem
        elif node.tag == 'consensusElement':
            counts.update({maps[e.get('map')] for e in node.findall('./groupedElementList/element')})
            node.clear()
        elif node.tag == 'UnassignedPeptideIdentification':
            node.clear()
            root.remove(node)
        elif node.tag in ('IdentificationRun', 'consensusElementList'):
            node.clear()
    assert verified_settings
    return counts


def ratio_counts(path, level):
    fields = ('sequence', 'modifications') if level == 'peptide' else ('protein_group',)
    entities, per_run, cells = set(), collections.defaultdict(set), 0
    for row in table(path):
        assert row['normalization'] == 'raw'
        if row['ratio_status'] != 'reported':
            continue
        assert int(row['numerator_channel']) in (2, 3) and int(row['reference_channel']) == 1
        identity = tuple(row[k] for k in fields)
        entities.add(identity)
        per_run[row['runs']].add(identity)
        cells += 1
    return len(entities), {k: len(v) for k, v in per_run.items()}, cells


def draw(rows):
    fig, axes = plt.subplots(2, 2, figsize=(12.8, 9.8))
    fig.subplots_adjust(left=.10, right=.97, top=.79, bottom=.28, hspace=.95, wspace=.32)
    fig.text(.07, .95, 'PXD039580 HeLa | coverage across completed runs', fontsize=21, weight='bold')
    fig.text(.07, .904, 'Original 3.6 build  →  3.7 baseline  →  3.7 + isotope requirement 3:6 to 2:6',
             fontsize=12.5, color='#415164')
    fig.text(.07, .868, 'Six injections: two triplex mixtures × three technical injections | Comet 0.02 Da throughout',
             fontsize=10.8, color='#566473')
    for ax, (key, heading, definition) in zip(axes.flat, METRICS):
        values = [row[key] for row in rows]
        ax.set_title(heading, loc='left', fontsize=16, weight='bold', pad=34)
        ax.text(0, 1.08, definition, transform=ax.transAxes, fontsize=10, color='#566473')
        for i, value in enumerate(values):
            if value is None:
                ax.text(i, max(v for v in values if v is not None) * .10, 'N/A\nlog unavailable',
                        ha='center', fontsize=10, color='#566473')
            else:
                bars = ax.bar(i, value, color=COLORS[i], width=.60, zorder=3)
                ax.bar_label(bars, labels=[f'{value:,}'], padding=6, fontsize=12, weight='bold')
        ax.set_xticks(range(3), LABELS, fontsize=10.5)
        ax.set_ylim(0, max(v for v in values if v is not None) * 1.28)
        ax.set_xlim(-.65, 2.65)
        ax.set_axisbelow(True)
        ax.grid(axis='y', color='#e5e9ef')
        ax.yaxis.set_major_formatter(StrMethodFormatter('{x:,.0f}'))
        ax.tick_params(axis='x', length=0, pad=9)
        gains = []
        for i, label in enumerate(('Upgrade', 'Isotope trial')):
            a, b = values[i:i + 2]
            if a is None:
                change = 'N/A'
            elif b == a:
                change = '0.0%'
            elif key == 'identified_target_psms' and abs(100 * (b / a - 1)) < .05:
                change = f'{b - a:+,} PSMs'
            else:
                change = f'{100 * (b / a - 1):+.1f}%'
            gains.append(f'{label} {change}')
        ax.text(0, -.40, '    |    '.join(gains), transform=ax.transAxes, fontsize=10,
                color='#1f606a', weight='bold')
    proteins = ' → '.join(f"{r['quantified_protein_groups']:,}" for r in rows)
    fig.text(.07, .155, f'Quantified protein groups:  {proteins}', fontsize=11, color='#415164')
    fig.text(.07, .116, 'The 2-isotope bar is an experimental sensitivity trial; 3:6 remains the controlled baseline.',
             fontsize=10.5, color='#8b5a17', weight='bold')
    fig.text(.07, .038,
             'Panels use different units. Mapped/retained counts sum over six injections; peptide/protein counts are experiment-wide unions.\n'
             'All runs: protein ratio support ≥2, MBR off, light reference, raw ratios. Only native M/L and H/L count toward coverage.\n'
             'Pinned prerelease builds; full-pipeline upgrade comparison. Original mapping count is unavailable, not zero.',
             fontsize=9.2, color='#566473', linespacing=1.5)
    return fig


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--outdir', type=Path, default=ROOT / 'results/pxd039580_coverage_stages_20261022')
    args = parser.parse_args()
    out = args.outdir.resolve()
    out.mkdir(parents=True, exist_ok=True)
    hashes, records, rows, injections, bundles = {}, {}, [], [], {}

    def fingerprint(path):
        digest = sha(path)
        hashes[str(path.relative_to(ROOT))] = digest
        return digest

    for name, analysis in ANALYSES.items():
        analysis = ROOT / analysis
        metrics = analysis / 'metrics.json'
        fingerprint(metrics)
        bundle, = json.loads(metrics.read_text())['bundles']
        assert bundle['status'] == 'analyzed' and not bundle['diagnostics']
        bundles[name] = bundle
        for source in bundle['inputs'].values():
            assert fingerprint(Path(source['path'])) == source['sha256']
        settings = bundle['processing']['parameters']
        assert settings['ratios:min_ratio_count'] == '2'
        assert settings['match_between_runs'] == settings['ratios:normalize'] == 'false'
        psms = count_psms(Path(bundle['inputs']['mztab']['path']))
        exported = count_exported(Path(bundle['inputs']['consensusxml']['path']), settings)
        assert dict(exported) == bundle['processing']['final_features_per_run']
        assert sum(psms.values()) == bundle['identification']['unique_resolved_psms']
        stages = None
        if name in LOGS:
            log = ROOT / LOGS[name]
            stages = extract_stages(log.read_text())
            fingerprint(log)
            shutil.copyfile(log, out / f'{name}.ms1.log')
            assert stages.keys() == psms.keys() == exported.keys()
            assert all(stages[run]['retained_multiplets'] == exported[run] for run in stages)
        row = dict(configuration=name, identified_target_psms=sum(psms.values()),
                   detected_multiplets=sum(r['detected_multiplets'] for r in stages.values()) if stages else None,
                   mapped_multiplets=sum(r['mapped_multiplets'] for r in stages.values()) if stages else None,
                   retained_multiplets=sum(exported.values()))
        run_ratios = {}
        for level, metric in [('peptide', 'quantified_peptide_forms'), ('protein', 'quantified_protein_groups')]:
            path = analysis / bundle['output_folder'] / f'{level}_ratios.tsv'
            fingerprint(path)
            row[metric], run_ratios[metric], row[f'{level}_ratio_cells'] = ratio_counts(path, level)
        rows.append(row)
        for run in sorted(psms):
            injections.append(dict(configuration=name, injection=run,
                mixture='Mix1' if 'A3R1' in run else 'Mix2', identified_target_psms=psms[run],
                mapped_multiplets=stages[run]['mapped_multiplets'] if stages else None,
                retained_multiplets=exported[run],
                **{k: v.get(run, 0) for k, v in run_ratios.items()}))
        records[name] = {'analysis': str(analysis.relative_to(ROOT)), 'inputs': bundle['inputs'],
                         'software': bundle['processing']['software'], 'ms1_parameters': settings,
                         'ms1_log': LOGS.get(name)}
        print(f"Verified {name}: {row}", flush=True)
    a, b = [bundles[n]['processing']['parameters'] for n in ('upgrade', 'isotopes2')]
    changes = {k: [a.get(k), b.get(k)] for k in a.keys() | b.keys() if a.get(k) != b.get(k)}
    assert changes == {'algorithm:isotopes_per_peptide': ['3:6', '2:6']}
    assert [(r['identified_target_psms'], r['quantified_peptide_forms'], r['quantified_protein_groups'])
            for r in rows] == [(272604, 6130, 694), (272600, 16940, 1538), (272600, 19964, 1737)]
    write_tsv(out / 'coverage_by_configuration.tsv', rows)
    write_tsv(out / 'coverage_by_injection.tsv', injections)
    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 11, 'axes.spines.top': False,
                         'axes.spines.right': False, 'pdf.fonttype': 42})
    fig = draw(rows)
    fig.savefig(out / 'hela_coverage.png', dpi=170, facecolor='white')
    fig.savefig(out / 'PXD039580_coverage_stages.pdf', facecolor='white')
    plt.close(fig)
    notes = ['# PXD039580 HeLa: coverage across completed runs', '',
             '[PNG](hela_coverage.png) · [PDF](PXD039580_coverage_stages.pdf)', '',
             'Same six injections throughout: two triplex mixtures, three technical injections each. '
             'Plots use the retained 3.6 original, 3.7 baseline and completed 2-isotope sensitivity trial. '
             'No scientific processing was launched.', '',
             '| Configuration | Target PSMs | Mapped multiplets | Retained multiplets | Quantified peptide forms | Quantified protein groups |',
             '|---|---:|---:|---:|---:|---:|']
    for r in rows:
        notes.append('| ' + ' | '.join('N/A' if r[k] is None else str(r[k]) for k in
            ('configuration', 'identified_target_psms', 'mapped_multiplets', 'retained_multiplets',
             'quantified_peptide_forms', 'quantified_protein_groups')) + ' |')
    notes += ['', '## Settings and interpretation', '',
        '- Original: OpenMS 3.6.0-pre-HEAD-2026-09-30; knock_out=false; precursor m/z mapping; requantify absent from the parameter record (not evidence of false).',
        '- Upgrade: OpenMS 3.7.0-pre-HEAD-2026-10-19; knock_out=true; requantify=true; peptide m/z mapping. Full-pipeline comparison also changes identification storage, input ordering and threads; individual causal contributions are not isolated.',
        '- Trial: same 3.7 build and parameters as the baseline, except algorithm:isotopes_per_peptide=2:6 instead of 3:6. Identification inputs were reused. Coverage rises, but Mix1 accuracy/precision tradeoffs mean 3:6 remains the controlled baseline.',
        '- Comet fragment setting is already 0.02 Da throughout; precursor tolerance is 4.5 ppm. This is different from the E. coli fragment-tolerance tuning comparison.',
        '- All plotted runs: protein ratio support 2, MBR off, normalization off, light reference. Native positive M/L and H/L establish coverage; derived H/M is excluded.', '',
        '## Units and validation', '',
        'Identified counts are final target PSM rows recounted from mzTab. Mapped counts come from matching MS1 task logs. '
        'Retained counts are exported run/multiplet pairs independently recounted from consensusXML and checked against logs where available. '
        'The original task directory and mapping log are unavailable; its mapped count is N/A, never estimated from another run. '
        'Its retained count comes from its exact saved consensusXML.', '',
        'Peptide/protein counts are distinct experiment-wide identities with a reported native raw ratio, recounted from the existing analysis exports. '
        'Peptide identities distinguish modifications. All analysis input hashes are checked against the actual mzTab, consensusXML and design files. '
        'Recorded XML settings and per-run retained counts are independently checked. The four panels have different units and do not form a numerical funnel.', '',
        'The isotope trial is not a strict superset: the existing review reports 4,415 gained and 1,391 lost peptide forms, '
        'with mixture-dependent quality changes. See [baseline review](../openms37_validation/hela_review/report.md) '
        'and [isotope trial review](../openms37_tuning/hela_detector_isotopes2/review/report.md). '
        'The two-file Thermo-native pilot differs in injection subset and settings and is excluded from this six-file comparison.', '',
        'Other completed branches are separate interventions: intensity 500 left quantitative evidence unchanged; '
        'support 1 increased protein coverage; MBR increased transfer-based coverage. They are not combined into the third bar.', '',
        'Reproduce from the repository root:', '', '```bash',
        '/home/philipp/miniconda3/envs/silac/bin/python bin/plot_pxd039580_coverage_stages.py', '```', '']
    (out / 'README.md').write_text('\n'.join(notes))
    fingerprint(Path(__file__).resolve())
    (out / 'provenance.json').write_text(json.dumps({'source_hashes': hashes, 'configurations': records,
        'trial_parameter_changes': changes, 'checks': 'Input hashes; target PSM recount; XML settings; exported run/multiplet recount; log retained counts; native ratio coverage recount.'}, indent=2) + '\n')
    print(f'Output: {out}', flush=True)


if __name__ == '__main__':
    main()
