#!/usr/bin/env python3
"""Make a concise coverage comparison from retained outputs and MS1 task logs."""
import argparse
import collections
import csv
import json
from pathlib import Path
import re
import shutil
import xml.etree.ElementTree as ET

from compare_pxd003327_runs import ROOT, export_ratios, read_mztab, sha, plt, PdfPages
from matplotlib.ticker import StrMethodFormatter

COLORS = ['#7b8797', '#d09236', '#137f88']
STAGES = ['old', 'upgrade', 'current']
LABELS = ['3.6.0\nOriginal', '3.7.0\nUpgrade', '3.7.0\n+ tuning']


def write_tsv(path, rows):
    with path.open('w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter='\t')
        writer.writeheader()
        writer.writerows(rows)


def extract_stages(text):
    result = {}
    for run, block in re.findall(r'Processing (\S+) \(fraction group.*?\)(.*?)(?=Processing |Progress of \'aligning maps\'|Building graph)', text, re.S):
        def number(pattern):
            match = re.search(pattern, block)
            if not match:
                raise ValueError(f'Missing log measurement: {run}: {pattern}')
            return int(match[1])
        result[Path(run).stem] = dict(
            detected_multiplets=number(r'Detected (\d+) peptide multiplet'),
            mapping_input_multiplets=number(r'Mapped identifications onto \d+ of (\d+)'),
            mapped_multiplets=number(r'Mapped identifications onto (\d+)'),
            retained_multiplets=number(r'Resolved (\d+) multiplet'),
            label_conflicts=number(r'; (\d+) identified multiplet'),
            input_psms_including_decoys=number(r'Mapping (\d+) peptides to'),
            unassigned_psms=number(r'Unassigned peptides: (\d+)'),
        )
        row = result[Path(run).stem]
        assert row['mapped_multiplets'] == row['retained_multiplets'] + row['label_conflicts']
    assert result
    return result


def draw(plex, rows):
    selected = [next(r for r in rows if r['configuration'] == f'{plex}_{s}') for s in STAGES]
    fig, axes = plt.subplots(2, 2, figsize=(12.8, 9.2))
    fig.subplots_adjust(left=.09, right=.96, top=.79, bottom=.26, hspace=.82, wspace=.30)
    fig.text(.07, .95, f'PXD003327 {plex} | what improved?', fontsize=22, weight='bold')
    fig.text(.07, .902, 'Original 3.6 build  →  3.7 upgrade  →  3.7 + 0.02 Da fragment tolerance', fontsize=13, color='#415164')
    subset = 'Two injections: E. coli 4:1 and 1:1' if plex == 'duplex' else 'One injection: E. coli 1:2:4'
    fig.text(.07, .865, subset, fontsize=11, color='#566473')
    metrics = [
        ('identified_target_psms', 'Identified', 'Target PSMs in the final identification output'),
        ('mapped_multiplets', 'Mapped', 'MS1 multiplets with mapped identifications'),
        ('retained_multiplets', 'Retained', 'Multiplets kept after label-consistency checks'),
        ('quantified_peptide_forms', 'Quantified', 'Distinct modified peptides with usable quantification'),
    ]
    for ax, (metric, heading, definition) in zip(axes.flat, metrics):
        values = [row[metric] for row in selected]
        ax.set_title(heading, loc='left', fontsize=16, weight='bold', pad=34)
        ax.text(0, 1.08, definition, transform=ax.transAxes, fontsize=10, color='#566473')
        bars = ax.bar(range(3), values, color=COLORS, width=.60, zorder=3)
        ax.bar_label(bars, labels=[f'{n:,}' for n in values], padding=6, fontsize=13, weight='bold')
        ax.set_xticks(range(3), LABELS, fontsize=11)
        ax.set_ylim(0, max(values) * 1.28)
        ax.set_xlim(-.6, 2.6)
        ax.set_axisbelow(True)
        ax.grid(axis='y', color='#e5e9ef')
        ax.yaxis.set_major_formatter(StrMethodFormatter('{x:,.0f}'))
        ax.tick_params(axis='x', length=0, pad=9)
        gains = [100 * (values[1] / values[0] - 1), 100 * (values[2] / values[1] - 1)]
        ax.text(0, -.34, f'Upgrade {gains[0]:+.1f}%    |    Further tuning {gains[1]:+.1f}%',
                transform=ax.transAxes, fontsize=10.5, color='#1f606a', weight='bold')
    detectors = ' → '.join(f"{row['detected_multiplets']:,}" for row in selected)
    proteins = ' → '.join(str(row['quantified_protein_groups']) for row in selected)
    fig.text(.07, .13, f'Detected MS1 multiplets:  {detectors}', fontsize=11, color='#415164')
    fig.text(.55, .13, f'Quantified protein groups:  {proteins}', fontsize=11, color='#415164')
    fig.text(.07, .055,
             'Each panel has its own unit; these bars are not a single attrition funnel. Multiplets group isotope channels.\n'
             'Mapped/retained counts are summed over injections. Peptide/protein counts are distinct identities across the subset.\n'
             'Pinned pre-release builds; retained completed runs. All stages: protein support ≥2, MBR off, automatic OpenMS picking.',
             fontsize=9.5, color='#566473', linespacing=1.5)
    return fig


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', type=Path, default=ROOT / 'results/pxd003327_version_comparison_20261021')
    parser.add_argument('--outdir', type=Path)
    args = parser.parse_args()
    source = args.source.resolve()
    out = args.outdir.resolve() if args.outdir else source / 'coverage_stages'
    out.mkdir(parents=True, exist_ok=True)
    manifest = json.loads((source / 'comparison.json').read_text())
    with (source / 'tasks.tsv').open() as handle:
        tasks = {r['configuration']: r for r in csv.DictReader(handle, delimiter='\t') if r['task'] == 'SILAC:MS1_LABELED_WORKFLOW'}
    rows, injection_rows, hashes = [], [], {}
    for name, record in manifest['runs'].items():
        base = ROOT / record['path']
        log = ROOT / tasks[name]['workdir'] / '.command.out'
        text = log.read_text()
        stage_counts = extract_stages(text)
        sections = read_mztab(base / 'experiment.mzTab')
        xmlpath = base / 'experiment.consensusXML'
        xml = ET.parse(xmlpath).getroot()
        maps = {m.get('id'): Path(m.get('name')).stem for m in xml.findall('./mapList/map')}
        exported = collections.Counter()
        for node in xml.findall('./consensusElementList/consensusElement'):
            exported.update({maps[e.get('map')] for e in node.findall('./groupedElementList/element')})
        design, = base.glob('*.sdrf_openms_design.tsv')
        peptides, proteins = export_ratios(base / 'experiment.mzTab', xmlpath, design)
        quantified_peptides = [r for r in peptides if r['ratio_status'] == 'reported']
        quantified_proteins = [r for r in proteins if r['ratio_status'] == 'reported']
        for run, counts in stage_counts.items():
            assert exported[run] == counts['retained_multiplets'], (name, run, exported[run], counts)
            psms = [r for r in sections['PSM'] if r['run'] == run]
            injection_rows.append(dict(configuration=name, injection=run, identified_target_psms=len(psms),
                identified_sequences=len({r['sequence'] for r in psms}), **counts,
                exported_multiplets=exported[run],
                quantified_peptide_forms=len({(r['sequence'], r['modifications']) for r in quantified_peptides if r['runs'] == run}),
                quantified_protein_groups=len({r['protein_group'] for r in quantified_proteins if r['runs'] == run})))
        row = dict(configuration=name, identified_target_psms=len(sections['PSM']),
                   **{k: sum(v[k] for v in stage_counts.values()) for k in next(iter(stage_counts.values()))},
                   exported_multiplet_run_pairs=sum(exported.values()),
                   quantified_peptide_forms=len({(r['sequence'], r['modifications']) for r in quantified_peptides}),
                   quantified_protein_groups=len({r['protein_group'] for r in quantified_proteins}))
        assert row['quantified_peptide_forms'] == record['coverage']['peptide_entities']
        assert row['quantified_protein_groups'] == record['coverage']['protein_entities']
        rows.append(row)
        for p in (log, base / 'experiment.mzTab', xmlpath, design):
            hashes[str(p.relative_to(ROOT))] = sha(p)
        shutil.copyfile(log, out / f'{name}.ms1.log')
    write_tsv(out / 'coverage_by_configuration.tsv', rows)
    write_tsv(out / 'coverage_by_injection.tsv', injection_rows)
    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 11, 'axes.spines.top': False,
                         'axes.spines.right': False, 'pdf.fonttype': 42})
    with PdfPages(out / 'PXD003327_coverage_stages.pdf') as pdf:
        pdf.infodict().update(Title='PXD003327: identification, mapping, retention and quantification')
        for plex in ('duplex', 'triplex'):
            fig = draw(plex, rows)
            fig.savefig(out / f'{plex}_coverage.png', dpi=170, facecolor='white')
            fig.savefig(out / f'{plex}_coverage.pdf', facecolor='white')
            pdf.savefig(fig, facecolor='white')
            plt.close(fig)
    notes = ['# PXD003327: simple coverage comparison', '',
             '[Two-page PDF](PXD003327_coverage_stages.pdf) · [Duplex PNG](duplex_coverage.png) · [Triplex PNG](triplex_coverage.png)', '',
             'Three completed stages: old 3.6.0 build / 0.8 Da; 3.7.0 build / 0.8 Da; 3.7.0 build / 0.02 Da. The last contrast shows the latest adopted E. coli search adjustment. No new processing or input data is needed for these plots.', '',
             '| Stage | Identified target PSMs | Mapped multiplets | Retained multiplets | Quantified peptide forms | Quantified proteins |',
             '|---|---:|---:|---:|---:|---:|']
    for r in rows:
        notes.append('| ' + ' | '.join(str(r[k]) for k in ('configuration', 'identified_target_psms', 'mapped_multiplets', 'retained_multiplets', 'quantified_peptide_forms', 'quantified_protein_groups')) + ' |')
    notes += ['', '## Definitions', '',
              '- **Identified:** final target PSM rows, recomputed from mzTab.',
              '- **Mapped:** detected MS1 multiplets with mapped identification(s), read from the MS1 log. This is a multiplet count, not a PSM count.',
              '- **Retained:** label-consistent multiplets kept by the resolver, read from the log and independently matched to exported consensusXML multiplet/run pairs.',
              '- **Quantified:** distinct modified peptide identities with at least one valid stored channel comparison. This counts coverage only; no ratio values or accuracy distributions are plotted.',
              '- **Detected:** detector-reported multiplets before partial-multiplet consolidation. The TSV also gives the post-consolidation mapping-input count. Detection alone does not establish identification or usable quantification.', '',
              'Duplex mapped/retained counts sum over the two physical injections. Quantified peptide and protein counts are unions across them. The four metrics use different units and must not be interpreted as one numerical funnel. `coverage_by_injection.tsv` supplies all individual injection counts, mapping-input counts, and label-conflict counts.', '',
              'The old counts come from the actual saved task logs matching the selected outputs. The older October 5 duplex checkpoint has tiny detector/retention differences (e.g. 434 rather than 433 retained multiplets for 4:1); those historical report counts are not substituted for the matching-run logs.', '',
              'The 3.7 upgrade includes partial-multiplet detection, requantification and peptide-m/z mapping behavior. The later 0.02 Da runs hold those settings fixed. Both current E. coli profiles retain automatic OpenMS picking, protein support 2 and MBR off. Other experimental branches (MS2-only picking, MBR, support 1) are not combined into an unexecuted configuration.', '',
              'These are retained completed runs, not new executions of today’s source. The previous report documents later code/configuration edits, including the now-explicit Comet maximum-variable-modification argument. See [full parameter record](../changes_and_parameters.md).', '',
              'The plots use OpenMS `3.6.0-pre-HEAD-2026-09-30` and `3.7.0-pre-HEAD-2026-10-19`. The exact source paths, SHA-256 hashes and copied MS1 logs are retained here; the parent report records immutable image digests and all task commands.', '',
              'Reproduce from the repository root:', '', '```bash',
              '/home/philipp/miniconda3/envs/silac/bin/python bin/plot_pxd003327_coverage_stages.py', '```', '']
    (out / 'README.md').write_text('\n'.join(notes))
    (out / 'provenance.json').write_text(json.dumps({'source_hashes': hashes, 'source_report': str(source),
        'checks': 'Mapped = retained + label conflicts; retained = exported run/multiplet pairs; quantified counts match original output audit.',
        'configurations': manifest['runs']}, indent=2) + '\n')
    print(json.dumps({'output': str(out), 'coverage': rows}, indent=2))


if __name__ == '__main__':
    main()
