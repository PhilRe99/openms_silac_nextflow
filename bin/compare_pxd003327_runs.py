#!/usr/bin/env python3
"""Plot retained PXD003327 runs and archive their executed settings (read-only inputs)."""
import argparse
import collections
import csv
import difflib
import hashlib
import json
import math
import os
from pathlib import Path
import re
import shutil
import sys
import textwrap
import xml.etree.ElementTree as ET

os.environ.setdefault('MPLCONFIGDIR', '/tmp/pxd003327-comparison-mpl')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'bin'))
from export_silac_ratios import export_ratios

STAGES = ('old', 'upgrade', 'current')
LABELS = ('3.6 build\n0.8 Da', '3.7 build\n0.8 Da', '3.7 build\n0.02 Da')
COLORS = ('#728296', '#d3933c', '#177e89')
OLD = 'sha256:824bcdf017b96f3f4509f68395efe1ce19014be7cc2b07263ae3cb2e652dac02'
NEW = 'sha256:734197c186b23a98bdb5cb25071d28e10b125480f8d67d8aa59c76c3632c5ac0'
PATHS = {
    'duplex_old': 'results/openms_version_tests/duplex_baseline/PXD003327_duplex',
    'triplex_old': 'results/id_format_native_triplex/PXD003327_triplex',
    'duplex_upgrade': 'results/openms37_validation/duplex/PXD003327_duplex',
    'triplex_upgrade': 'results/openms37_validation/triplex/PXD003327_triplex',
    'duplex_current': 'results/openms37_validation/duplex_002da/PXD003327_duplex_002da',
    'triplex_current': 'results/openms37_validation/triplex_002da/PXD003327_triplex_002da',
}
NAMES = dict(zip(PATHS, ('small_galileo', 'ecstatic_meucci', 'naughty_legentil',
                        'naughty_linnaeus', 'compassionate_swanson', 'lonely_darwin')))


def sha(path):
    h = hashlib.sha256()
    with path.open('rb') as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def tsv(path, rows):
    if not rows:
        return
    with path.open('w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter='\t')
        writer.writeheader()
        writer.writerows(rows)


def read_mztab(path):
    headers, sections, names = {}, collections.defaultdict(list), {}
    with path.open() as handle:
        for fields in csv.reader(handle, delimiter='\t'):
            if fields and fields[0] == 'MTD' and re.fullmatch(r'ms_run\[\d+\]-location', fields[1]):
                names[re.search(r'\d+', fields[1])[0]] = Path(fields[2]).stem
            elif fields and fields[0] in ('PSH', 'PEH', 'PRH'):
                headers[{'PSH': 'PSM', 'PEH': 'PEP', 'PRH': 'PRT'}[fields[0]]] = fields[1:]
            elif fields and fields[0] in ('PSM', 'PEP', 'PRT'):
                assert len(fields[1:]) == len(headers[fields[0]])
                sections[fields[0]].append(dict(zip(headers[fields[0]], fields[1:])))
    for row in sections['PSM']:
        assert row['opt_global_cv_MS:1002217_decoy_peptide'] == '0'
        row['run'] = names[re.search(r'ms_run\[(\d+)\]', row['spectra_ref'])[1]]
    return sections


def ratio_key(row, level):
    fields = ('sequence', 'modifications') if level == 'peptide' else ('protein_group',)
    return tuple(row[f] for f in fields + ('fraction_group', 'numerator_channel', 'reference_channel', 'normalization'))


def target(row):
    seq = row['sequence']
    aa = 'K' if seq.count('K') == 1 and seq.count('R') == 0 else 'R' if seq.count('R') == 1 and seq.count('K') == 0 else None
    if aa is None:
        return None
    run, ch = row['runs'], row['numerator_channel']
    expected = (1. if run == 'Chris_Ecoli_1-1' else (.25 if aa == 'K' else 4.)) if ch == 2 and run != 'Chris_Ecoli_1-2-4' else ({2: 2., 3: 4.} if aa == 'K' else {2: .5, 3: .25})[ch]
    return run, ch, aa, expected


def accuracy(rows):
    errors = np.array([math.log2(r['ratio'] / target(r)[3]) for r in rows])
    return dict(n=len(rows), median_ratio=float(np.median([r['ratio'] for r in rows])) if rows else None,
                median_log2_error=float(np.median(errors)) if rows else None,
                median_absolute_log2_error=float(np.median(abs(errors))) if rows else None,
                iqr_log2_error=float(np.percentile(errors, 75) - np.percentile(errors, 25)) if rows else None,
                over_twofold_error_fraction=float(np.mean(abs(errors) > 1)) if rows else None)


def archived_tasks(name, tolerance):
    if name.endswith('_current') or name.endswith('_upgrade'):
        key = name.replace('_current', '_002').replace('_upgrade', '_08')
        return tolerance['runs'][key]['tasks']
    if name == 'duplex_old':
        with (ROOT / 'results/openms_version_tests/duplex_baseline_trace.tsv').open() as handle:
            rows = list(csv.DictReader(handle, delimiter='\t'))
        tasks = []
        for row in rows:
            shard, prefix = row['hash'].split('/')
            directory, = (ROOT / 'work/openms_version_tests' / shard).glob(prefix + '*')
            assert row['status'] == 'COMPLETED' and row['exit'] == '0'
            tasks.append(dict(name=row['name'], workdir=str(directory.relative_to(ROOT)), status=row['status']))
        return tasks
    log = (ROOT / 'results/id_format_native_comparison/triplex.native.nextflow.log').read_text()
    tasks = [dict(name=n, workdir=str(Path(p).relative_to(ROOT)), status=s) for n, s, code, p in re.findall(
        r'Task completed > TaskHandler\[id: \d+; name: (.*?); status: (.*?); exit: (.*?); error: .*?; workDir: (.*?)\]', log) if code == '0']
    for shard, prefix, task in re.findall(r'\[([0-9a-f]{2})/([0-9a-f]+)\] Cached process > (.*)', log):
        directory, = (ROOT / 'work/id_format_baseline' / shard).glob(prefix + '*')
        tasks.append(dict(name=task.strip(), workdir=str(directory.relative_to(ROOT)), status='CACHED'))
    return tasks


def title(fig, heading, subheading):
    fig.suptitle(heading, x=.07, y=.97, ha='left', fontsize=18, fontweight='bold')
    fig.text(.07, .918, subheading, fontsize=10, color='#465464')


def footer(fig, text):
    fig.text(.07, .035, text, fontsize=9, color='#465464', va='bottom', linespacing=1.45)


def text_page(pdf, heading, paragraphs):
    fig = plt.figure(figsize=(11.7, 8.3))
    title(fig, heading, 'PXD003327 | retained local runs | report date: 2026-10-21')
    y = .85
    for paragraph in paragraphs:
        lines = '\n'.join(textwrap.wrap(paragraph, width=112, break_long_words=False))
        fig.text(.07, y, lines, fontsize=10.5, va='top', linespacing=1.35)
        y -= .026 * (lines.count('\n') + 1) + .025
    assert y > .04, f'Text page overflow: {heading}'
    if heading == 'Coverage and SILAC ratio comparison':
        fig.savefig(Path(pdf._filename).parent / 'report_overview.png', dpi=130)
    pdf.savefig(fig)
    plt.close(fig)


def make_plots(out, runs, coverage, comparisons):
    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 10, 'axes.spines.top': False,
                         'axes.spines.right': False, 'pdf.fonttype': 42})
    with PdfPages(out / 'PXD003327_comparison.pdf') as pdf:
        pdf.infodict().update(Title='PXD003327: OpenMS 3.6 versus current profile settings', Author='NF-SILAC analysis')
        text_page(pdf, 'Coverage and SILAC ratio comparison', [
            'Three stages: old OpenMS 3.6 build with 0.8 Da fragments; pinned 3.7 development build with the same tolerance; then 3.7 with 0.02 Da fragments. Duplex contains two different mixtures; triplex contains one injection. These are not biological replicates.',
            'Current means the E. coli test profiles read on 2026-10-21: OpenMS picking enabled, automatic MS levels, partial detection and requantification enabled, peptide m/z mapping, MBR off, minimum protein ratio support 2, light reference and no ratio normalization.',
            'The 0.02 Da outputs are completed retained runs with those effective MS1 settings. They are not a fresh execution of today\'s source. Today\'s Comet configuration also explicitly sets max_variable_mods_in_peptide=5; the retained commands leave it implicit. The archived artifacts do not independently establish that old default.',
            'Duplex peptide forms with ratios: 701 to 963 (+37.4%); protein groups: 117 to 139 (+18.8%). Triplex peptide forms: 162 to 432 (+166.7%); protein groups: 38 to 78 (+105.3%). Greater coverage does not remove the persistent residue-specific ratio bias.',
            'Ratios are stored positive finite native raw ratios, without imputation or reconstruction from summed intensities. Benchmark classes contain exactly one K and no R, or one R and no K. A pooled protein expected ratio is inappropriate for these mixtures.',
            'The accompanying changes_and_parameters.md, settings TSVs, command snapshots, source hashes and comparison JSON contain the exact evidence. Original scientific results and workflow configuration were read without modification.'
        ])
        for plex in ('duplex', 'triplex'):
            fig, axes = plt.subplots(2, 3, figsize=(13.5, 8.5))
            fig.subplots_adjust(top=.80, bottom=.16, left=.08, right=.97, hspace=.66, wspace=.4)
            title(fig, f'{plex.title()} | identification and ratio coverage', 'All three stages use native idparquet identification bundles; counts are recomputed from the final outputs.')
            fields = [('target_psms', 'Target PSM rows'), ('identified_forms', 'Identified modified forms'),
                      ('peptide_entities', 'Peptide forms with ratios'), ('protein_entities', 'Protein groups with ratios'),
                      ('peptide_cells', 'Peptide ratio cells'), ('protein_cells', 'Protein ratio cells')]
            for ax, (field, label) in zip(axes.flat, fields):
                values = [runs[f'{plex}_{s}']['coverage'][field] for s in STAGES]
                bars = ax.bar(range(3), values, color=COLORS, width=.65)
                ax.bar_label(bars, labels=[f'{v:,}' for v in values], padding=4, fontsize=10)
                ax.set_xticks(range(3), LABELS, fontsize=9)
                ax.set_ylim(0, max(values) * 1.27)
                ax.set_title(label, loc='left', fontsize=11, pad=10)
                ax.set_axisbelow(True)
                ax.grid(axis='y', alpha=.2)
                ax.text(.98, .94, f'{100 * (values[-1] / values[0] - 1):+.1f}% vs old', transform=ax.transAxes,
                        ha='right', va='top', color=COLORS[-1], fontsize=9)
            footer(fig, 'Forms distinguish modifications. Entity counts are unions across injections/channels; ratio cells include fraction group and channel.\nCoverage denotes identification or usable-ratio yield, not amino-acid sequence coverage or a fixed-denominator sensitivity estimate.')
            fig.savefig(out / f'{plex}_coverage.png', dpi=150)
            pdf.savefig(fig)
            plt.close(fig)
        for plex in ('duplex', 'triplex'):
            classes = sorted(runs[f'{plex}_old']['classes'])
            fig, axes = plt.subplots(2, 2, figsize=(13.5, 9))
            fig.subplots_adjust(top=.80, bottom=.16, left=.09, right=.96, hspace=.6, wspace=.3)
            title(fig, f'{plex.title()} | known-mixture peptide ratio distributions', 'Same log2 observed/expected scale in every panel; 0 = expected, +1 = twice expected, -1 = half expected.')
            for ax, cls in zip(axes.flat, classes):
                samples = [[math.log2(r['ratio'] / cls[3]) for r in runs[f'{plex}_{s}']['classes'][cls]] for s in STAGES]
                boxes = ax.boxplot(samples, patch_artist=True, showfliers=True, widths=.5,
                                   flierprops=dict(marker='.', markersize=3, alpha=.35),
                                   medianprops=dict(color='black', linewidth=1.4))
                for patch, color in zip(boxes['boxes'], COLORS):
                    patch.set_facecolor(color)
                    patch.set_alpha(.75)
                ax.axhline(0, color='#284b35', linewidth=1)
                ax.axhspan(-1, 1, color='#75aa86', alpha=.08)
                run, ch, aa, expected = cls
                ax.set_title(f"{run.replace('Chris_Ecoli_', '')} | single-{aa} | {'M/L' if ch == 2 and plex == 'triplex' else 'H/L'} = {expected:g}", loc='left', fontsize=11)
                ax.set_xticks([1, 2, 3], [f'{label}\nn={len(v)}; med={np.median([r["ratio"] for r in runs[f"{plex}_{s}"]["classes"][cls]]):.3g}' for s, label, v in zip(STAGES, LABELS, samples)], fontsize=8)
                ax.set_ylim(-4, 4)
                ax.set_ylabel('log2(observed / expected)')
                ax.grid(axis='y', alpha=.2)
                tails = sum(sum(abs(x) > 4 for x in v) for v in samples)
                if tails:
                    ax.text(.02, .97, f'{tails} points beyond axis; retained in statistics', transform=ax.transAxes, va='top', fontsize=8)
            footer(fig, 'Boxes: median and interquartile range; whiskers: 1.5 x IQR; dots: outliers. Each stage includes all eligible reported peptide forms.\nPopulation composition changes between stages. Exact benchmark summaries and shared-only populations are supplied as TSV files.')
            fig.savefig(out / f'{plex}_ratio_distributions.png', dpi=150)
            pdf.savefig(fig)
            plt.close(fig)
        for level in ('peptide', 'protein'):
            fig, axes = plt.subplots(1, 2, figsize=(13.5, 7.5))
            fig.subplots_adjust(top=.79, bottom=.23, left=.08, right=.97, wspace=.28)
            title(fig, f'Shared {level} ratios | direct old versus current comparison', 'Each point is the identical entity / fraction group / numerator / reference ratio cell in both outputs.')
            for plex, ax in zip(('duplex', 'triplex'), axes):
                a, b = [runs[f'{plex}_{s}']['maps'][level] for s in ('old', 'current')]
                keys = sorted(a.keys() & b.keys())
                x = np.array([math.log2(a[k]['ratio']) for k in keys])
                y = np.array([math.log2(b[k]['ratio']) for k in keys])
                ax.scatter(x, y, s=10, color=COLORS[-1], alpha=.4, rasterized=True)
                bounds = [math.floor(min(x.min(), y.min())) - .5, math.ceil(max(x.max(), y.max())) + .5]
                ax.plot(bounds, bounds, color='#3d4752', linewidth=1)
                ax.set(xlim=bounds, ylim=bounds, xlabel='Old log2 ratio', ylabel='Current log2 ratio', title=plex.title())
                ax.set_aspect('equal', adjustable='box')
                c = comparisons[f'{plex}_old_to_current'][level]
                ax.text(.03, .97, f"Shared {c['shared']:,}\nGained {c['gained']:,}; lost {c['lost']:,}\nShared changes >2-fold: {c['over_twofold_change']}", transform=ax.transAxes, va='top', fontsize=9)
                ax.grid(alpha=.15)
            footer(fig, 'Gained/lost refer to exact ratio-cell identities, including exact protein-group membership. Protein-group splits/merges count as turnover.\nProtein ratios can shift when their supporting K/R peptide composition changes; the diagonal tests agreement, not biological accuracy.')
            fig.savefig(out / f'{level}_paired_ratios.png', dpi=150)
            pdf.savefig(fig)
            plt.close(fig)
        text_page(pdf, 'Executed changes | version and MS1 behavior', [
            'Old image: OpenMS 3.6.0-pre-HEAD-2026-09-30, revision 5d5cbff. Upgraded image: 3.7.0-pre-HEAD-2026-10-19, revision 842adef. Full immutable container digests are listed in the companion document and manifest.',
            'Version upgrade: algorithm:knock_out false -> true; id_mapping:mz_reference precursor -> peptide; requantify was absent from the old parameter record and is true in the new record. Absence does not establish a recorded old false value. All other recorded MS1 scientific parameter values agree.',
            'Search tuning: Comet fragment_mass_tolerance 0.8 -> 0.02 Da. Precursor tolerance remains 10 ppm. The 3.7 runs before and after tolerance tuning have byte-identical experimental designs. The saved resolved configurations and executed tool commands isolate this tolerance change and metadata filenames.',
            'The version contrast includes the changed tool build and its inherited defaults. It does not identify the independent contribution of each algorithm change. The 0.8 versus 0.02 Da contrast isolates the search setting but still includes ordinary identification/rescoring rerun variability.',
            'Unchanged MS1 settings include reference_channel=1, min_ratio_count=2, normalize=false, match_between_runs=false, proteinFDR=0.01, picked_proteinFDR=true, protein_inference=aggregation and protein_quantification=unique_peptides.',
            'Percolator applies score:fdr=0.01 with score_type=pep and the external subprocess backend. MS1 records psmFDR=1.0; that must not be mislabeled as the upstream Percolator cutoff. Configured thresholds are not independently measured error rates.'
        ])
        text_page(pdf, 'Current configuration | scope and interpretation', [
            'The current profiles explicitly enable OpenMS PeakPickerHiRes for these supplied profile-mzML E. coli fixtures. peakpicking_ms_levels=null omits the MS-level override, matching the retained automatic-picking route. The global openms_peakpicking=false default does not override the test profiles.',
            'No Thermo RAW conversion was used for these E. coli results. The newer native Thermo centroiding path, the MS2-only experiment, support-1 runs, MBR trials and isotope-2 trials are not included as current results because they do not match the selected E. coli profiles.',
            'Identification is native idparquet at all three selected stages, so storage format is held constant. Earlier idXML migration results are documented separately; the accepted legacy duplex also has 701 peptide forms and 117 protein groups with ratios.',
            'Full per-run MS1 settings, search configurations, exact task commands and command diffs are archived with the report. Today\'s configuration and relevant workflow modules are snapshotted for auditing. New explicit Comet max-variable-modification configuration has not been retrospectively certified for the retained outputs.',
            'The benchmark remains a small local regression subset. The local PXD003327 label has a documented PXD003845 accession provenance ambiguity. Expectations here are tied to the actual local Chris_Ecoli file names and curated mixture definitions.',
            'No new mass-spectrometry processing was launched. This PDF reports existing measurements, not a fresh run of today\'s source. Use changes_and_parameters.md for the parameter table, source paths, exact invocation history and reproduction command.'
        ])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--outdir', type=Path, default=ROOT / 'results/pxd003327_version_comparison_20261021')
    args = parser.parse_args()
    out = args.outdir.resolve()
    out.mkdir(parents=True, exist_ok=True)
    tolerance = json.loads((ROOT / 'results/openms37_validation/tolerance_review/comparison.json').read_text())
    runs, coverage, per_injection, benchmark, settings, task_rows = {}, [], [], [], [], []
    manifest = {'analysis_date': '2026-10-21', 'runs': {}, 'source_hashes': {}, 'current_config_snapshot': {}}
    # Freeze the current source, without modifying user configuration or original results.
    sources = ['nextflow.config', 'conf/modules.config', 'workflows/silac.nf', 'assets/silac_chemistry.json',
               'modules/local/openms/comet/main.nf', 'modules/local/openms/openms_peak_picker/main.nf',
               'bin/export_silac_ratios.py', 'bin/compare_pxd003327_runs.py']
    sources += [f'testdata/PXD003327_{plex}{suffix}.sdrf.tsv' for plex in ('duplex', 'triplex') for suffix in ('', '_002da')]
    for relative in sources:
        source = ROOT / relative
        dest = out / 'source_snapshot' / relative
        dest.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(source, dest)
        manifest['current_config_snapshot'][relative] = sha(source)
    history = (ROOT / '.nextflow/history').read_text().splitlines()
    (out / 'invocations.tsv').write_text('\n'.join(line for line in history if any(f'\t{name}\t' in line for name in NAMES.values())) + '\n')
    for name, relative in PATHS.items():
        base = ROOT / relative
        plex, stage = name.split('_')
        design, = base.glob('*.sdrf_openms_design.tsv')
        xml = ET.parse(base / 'experiment.consensusXML').getroot()
        proc, = [p for p in xml.findall('dataProcessing') if p.find('software').get('name') == 'MS1LabeledWorkflow']
        params = {p.get('name').removeprefix('parameter: '): p.get('value') for p in proc.findall('UserParam') if p.get('name').startswith('parameter: ')}
        peptide, protein = export_ratios(base / 'experiment.mzTab', base / 'experiment.consensusXML', design)
        section = read_mztab(base / 'experiment.mzTab')
        maps = {}
        for level, rows in [('peptide', peptide), ('protein', protein)]:
            maps[level] = {ratio_key(r, level): r for r in rows if r['ratio_status'] == 'reported'}
            assert len(maps[level]) == sum(r['ratio_status'] == 'reported' for r in rows)
            assert all(r['normalization'] == 'raw' for r in rows)
            assert all(math.isfinite(r['ratio']) and r['ratio'] > 0 for r in maps[level].values())
            tsv(out / f'{name}.{level}_ratios.tsv', rows)
        c = dict(configuration=name, run_name=NAMES[name], target_psms=len(section['PSM']),
                 identified_forms=len({(r['sequence'], r['modifications']) for r in section['PSM']}))
        for level in maps:
            n = 2 if level == 'peptide' else 1
            c[f'{level}_entities'] = len({k[:n] for k in maps[level]})
            c[f'{level}_cells'] = len(maps[level])
        coverage.append(c)
        classes = collections.defaultdict(list)
        for row in maps['peptide'].values():
            if (cls := target(row)):
                classes[cls].append(row)
        for cls, rows in sorted(classes.items()):
            benchmark.append(dict(configuration=name, run=cls[0], channel=cls[1], residue_class=cls[2], expected=cls[3], **accuracy(rows)))
        for run in sorted({r['run'] for r in section['PSM']}):
            psms = [r for r in section['PSM'] if r['run'] == run]
            per_injection.append(dict(configuration=name, injection=run, target_psms=len(psms),
                identified_forms=len({(r['sequence'], r['modifications']) for r in psms}),
                peptide_ratio_cells=sum(r['runs'] == run for r in maps['peptide'].values()),
                protein_ratio_cells=sum(r['runs'] == run for r in maps['protein'].values())))
        task_list = archived_tasks(name, tolerance)
        assert len(task_list) == (15 if plex == 'duplex' else 9)
        tasks = {}
        for task in task_list:
            directory = ROOT / task['workdir']
            command = directory / '.command.sh'
            text = command.read_text()
            is_openms = not any(s in task['name'] for s in ('SDRF_PARSING', 'RESOLVE_SILAC_CONFIG'))
            if is_openms:
                assert (OLD if stage == 'old' else NEW) in (directory / '.command.run').read_text()
            command_name = re.sub(r'[^a-zA-Z0-9_.-]', '_', task['name']) + '.sh'
            destination = out / 'commands' / name / command_name
            destination.parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(command, destination)
            tasks[task['name']] = text
            task_rows.append(dict(configuration=name, task=task['name'], status=task['status'], workdir=task['workdir'], command_sha256=sha(command)))
        configs = {p.name: json.loads(p.read_text()) for p in base.glob('configuration/*.json')}
        if not configs:  # The original duplex publishes these under configuration as well.
            raise ValueError(f'No resolved search configs: {base}')
        for path in [base / 'experiment.mzTab', base / 'experiment.consensusXML', design, *base.glob('configuration/*.json')]:
            manifest['source_hashes'][str(path.relative_to(ROOT))] = sha(path)
        for parameter, value in params.items():
            settings.append(dict(configuration=name, parameter=parameter, value=value))
        runs[name] = dict(coverage=c, maps=maps, classes=classes, params=params, configs=configs, tasks=tasks, design=design)
        manifest['runs'][name] = dict(path=relative, run_name=NAMES[name], software=proc.find('software').get('version'),
                                    image=OLD if stage == 'old' else NEW, coverage=c, ms1_settings=params, resolved_search_configs=configs)
    comparisons, populations, deltas = {}, [], []
    ignored = {'design', 'in', 'ids', 'out', 'out_cxml', 'out_qpx', 'log', 'debug', 'threads', 'no_progress', 'force', 'test'}
    for plex in ('duplex', 'triplex'):
        for first, second in [('old', 'upgrade'), ('upgrade', 'current'), ('old', 'current')]:
            a, b = [runs[f'{plex}_{s}'] for s in (first, second)]
            pair = f'{plex}_{first}_to_{second}'
            comparison = {}
            for param in sorted(a['params'].keys() | b['params'].keys()):
                av, bv = a['params'].get(param), b['params'].get(param)
                if av != bv:
                    deltas.append(dict(comparison=pair, parameter=param, before=av if av is not None else '[absent]',
                                       after=bv if bv is not None else '[absent]', kind='execution/path' if param in ignored else 'scientific'))
            assert a['design'].read_bytes() == b['design'].read_bytes()
            for level in ('peptide', 'protein'):
                x, y = a['maps'][level], b['maps'][level]
                common = x.keys() & y.keys()
                comparison[level] = dict(shared=len(common), gained=len(y.keys() - x.keys()), lost=len(x.keys() - y.keys()),
                    unchanged=sum(x[k]['ratio'] == y[k]['ratio'] for k in common),
                    over_twofold_change=sum(abs(math.log2(y[k]['ratio'] / x[k]['ratio'])) > 1 for k in common))
                if level == 'peptide':
                    for population, rows in [('shared_before', [x[k] for k in common]), ('shared_after', [y[k] for k in common]),
                                             ('gained', [y[k] for k in y.keys() - x.keys()]), ('lost', [x[k] for k in x.keys() - y.keys()])]:
                        classes = collections.defaultdict(list)
                        for row in rows:
                            if (cls := target(row)):
                                classes[cls].append(row)
                        for cls, items in sorted(classes.items()):
                            populations.append(dict(comparison=pair, population=population, run=cls[0], channel=cls[1], residue_class=cls[2], expected=cls[3], **accuracy(items)))
            comparisons[pair] = comparison
            diff = []
            for task in sorted(a['tasks'].keys() & b['tasks'].keys()):
                diff.extend(difflib.unified_diff(a['tasks'][task].splitlines(True), b['tasks'][task].splitlines(True),
                                               fromfile=f'{first}/{task}', tofile=f'{second}/{task}'))
            (out / f'{pair}.commands.diff').write_text(''.join(diff))
    for filename, rows in [('coverage.tsv', coverage), ('coverage_by_injection.tsv', per_injection), ('benchmark.tsv', benchmark),
                           ('population_accuracy.tsv', populations), ('ms1_parameters.tsv', settings), ('parameter_changes.tsv', deltas), ('tasks.tsv', task_rows)]:
        tsv(out / filename, rows)
    manifest['comparisons'] = comparisons
    (out / 'comparison.json').write_text(json.dumps(manifest, indent=2, allow_nan=False) + '\n')
    make_plots(out, runs, coverage, comparisons)
    write_document(out, runs, coverage, benchmark, deltas)
    print(json.dumps({'output': str(out), 'coverage': coverage, 'comparisons': comparisons}, indent=2))


def write_document(out, runs, coverage, benchmark, deltas):
    lines = ['# PXD003327: OpenMS 3.6 build versus current E. coli profile settings', '',
             'Generated 2026-10-21 from retained outputs. No pipeline processing was launched; source results and configuration were not modified.', '',
             '## Scope and definition of current', '',
             '“Current” uses the completed 0.02 Da runs whose effective MS1 settings match the current `test_duplex` and `test_triplex` profiles: support 2, MBR off, partial detection/requantification on, peptide m/z mapping, automatic OpenMS picking, and raw light-reference ratios. These outputs are **not a fresh execution of today’s source**. The current Comet configuration explicitly supplies `-max_variable_mods_in_peptide 5`; the retained commands omit this option. Its earlier implicit value is not independently established here. Source snapshots preserve this distinction.', '',
             'Both profiles set `openms_peakpicking=true`, overriding the global false default. `peakpicking_ms_levels=null` omits the MS-level override. No Thermo RAW conversion was used. Support-1, MBR-on, MS2-only and isotope-2 trials do not represent these selected profile settings.', '',
             '## Selected runs', '', '| Configuration | Nextflow run | Result directory |', '|---|---|---|']
    lines += [f'| {name} | {NAMES[name]} | `{path}` |' for name, path in PATHS.items()]
    lines += ['', 'The old triplex invocation reused four cached preparation/metadata tasks and executed five identification/MS1 tasks. The other selected invocations completed all 15 duplex or 9 triplex tasks without caching. `invocations.tsv` preserves the exact launch commands, duration, session and status from Nextflow history.', '',
              '## Exact software and changed scientific settings', '',
              f'- Old: `3.6.0-pre-HEAD-2026-09-30`, revision `5d5cbff`; image `ghcr.io/openms/openms-tools-thirdparty@{OLD}`.',
              f'- Upgrade/current: `3.7.0-pre-HEAD-2026-10-19`, revision `842adef`; image `ghcr.io/openms/openms-tools-thirdparty@{NEW}`.', '',
              '| Setting | Old 3.6 build | Upgrade, 0.8 Da | Current-profile retained run |', '|---|---|---|---|',
              '| Comet fragment tolerance | 0.8 Da | 0.8 Da | 0.02 Da |',
              '| `algorithm:knock_out` (partial multiplets) | false | true | true |',
              '| `requantify` | absent from parameter record | true | true |',
              '| `id_mapping:mz_reference` | precursor | peptide | peptide |', '',
              'An absent old `requantify` entry is not an observed false value. Other recorded MS1 scientific parameters agree. `parameter_changes.tsv` contains every changed/added/removed MS1 key including output/path metadata. `ms1_parameters.tsv` contains **all** recorded keys and values for all six runs. `comparison.json` also contains every resolved search configuration, including exact modification lists and binary-modification vectors.', '',
              'The version step changes the whole OpenMS build and inherited MS1 defaults. Its gain cannot be attributed to one algorithm option. The tolerance step retains the same pinned image and MS1 settings. Existing command/configuration comparison evidence isolates the fragment tolerance and renamed SDRF/design filenames; experimental design contents remain byte-identical.', '',
              '## Shared settings and exact parameter values', '',
              '| Parameter / operation | Value in the retained comparisons |', '|---|---|',
              '| Precursor tolerance | 10 ppm |', '| Enzyme / missed cleavages | Trypsin / 2 |',
              '| Instrument / isotope error | high_res / 0/1/2/3 |',
              '| Fixed modification | Carbamidomethyl (C) |',
              '| Non-label variable modifications | Acetyl (Protein N-term); Deamidated (N); Deamidated (Q); Oxidation (M) |',
              '| Duplex labels | `[][Lys8,Arg10]` |', '| Triplex labels | `[][Lys4,Arg6][Lys8,Arg10]` |',
              '| Decoys | reverse; `DECOY_` prefix |', '| Identification storage | complete native idparquet bundles at all selected stages |',
              '| Percolator | external subprocess; `score_type=pep`; `score:fdr=0.01` |',
              '| MS1 `psmFDR` / `FDR_type` | 1.0 / PSM (separate from upstream 0.01 Percolator cutoff) |',
              '| Protein FDR / picked | 0.01 / true |', '| Inference / quantification | aggregation / unique_peptides |',
              '| `max_nr_labelled_aas` | 2 |', '| Reference channel / normalize | 1 (light) / false |',
              '| `ratios:min_ratio_count` | 2 |', '| `match_between_runs` | false |',
              '| Fraction aggregation | sum |', '| Detector charge / isotopes | 1:4 / 3:6 |',
              '| Detector m/z tolerance | 6 ppm |', '| Intensity cutoff | 1000 |',
              '| Peptide / averagine similarity / scaling | 0.5 / 0.4 / 0.95 |',
              '| ID mapping RT / m/z tolerance | 5 / 20 ppm |', '| Threads | 2 for configured multithreaded tasks |', '',
              'Label variable modifications: duplex K8 = `Label:13C(6)15N(2) (K)` and R10 = `Label:13C(6)15N(4) (R)`; triplex additionally K4 = `Label:2H(4) (K)` and R6 = `Label:13C(6) (R)`. Exact ordered lists and binary groups are preserved in `comparison.json` and task command files. Search target database: sibling `openms_silac_bash/db/ecoli_k12_uniprot_reviewed.fasta`; decoy FASTA is produced by DecoyDatabase.', '',
              '## Recomputed coverage', '', '| Configuration | Target PSM rows | Identified forms | Peptide forms with ratios | Protein groups with ratios | Peptide cells | Protein cells |', '|---|---:|---:|---:|---:|---:|---:|']
    for row in coverage:
        lines.append('| ' + ' | '.join(str(row[k]) for k in ('configuration', 'target_psms', 'identified_forms', 'peptide_entities', 'protein_entities', 'peptide_cells', 'protein_cells')) + ' |')
    lines += ['', 'Peptide forms distinguish modifications; entities are unions across fraction groups/channels. Ratio cells identify entity, fraction group, numerator, reference and normalization. These are yield counts, not protein amino-acid coverage. Protein group membership is matched exactly; splits and merges count as changed identities.', '',
              '## Known-mixture ratio accuracy', '',
              'Only exactly one K/no R or one R/no K. Expected H/L in duplex 4:1: K=0.25, R=4; duplex 1:1: both 1. Triplex M/L: K=2, R=0.5; H/L: K=4, R=0.25. Expectations are the project’s curated local mixture definitions, not newly inferred public-dataset claims.', '',
              '| Configuration | Injection | Channel | Class | Expected | n | Median observed | Median absolute log2 error |', '|---|---|---:|---|---:|---:|---:|---:|']
    for row in benchmark:
        lines.append(f"| {row['configuration']} | {row['run']} | {row['channel']} | {row['residue_class']} | {row['expected']:g} | {row['n']} | {row['median_ratio']:.5f} | {row['median_absolute_log2_error']:.5f} |")
    lines += ['', 'Triplex K M/L remains about 23% low despite higher coverage. Shared-only, gained and lost peptide populations are summarized separately in `population_accuracy.tsv`. No statistical significance is claimed from these few injections; mixtures are not replicates. Protein ratios have no universal expected value here because K and R have different expectations.', '',
              '## Configuration changes after these completed runs', '',
              '- `conf/modules.config` now explicitly forwards partial detection, requantification, ID m/z reference and MBR; their current values equal the recorded new-build defaults in the selected outputs.',
              '- Test profiles now select `*_002da.sdrf.tsv`. The optional global `comet_fragment_mass_tolerance` override remains null; the 0.02 Da value comes from metadata.',
              '- OpenMS picking is optional globally, with explicit true in the E. coli profiles; automatic MS-level selection preserves the evaluated route.',
              '- Comet now explicitly forwards `comet_max_variable_mods=5`, absent as an explicit argument in retained commands. Its previous implicit value is unresolved in this analysis.',
              '- RAW preparation now supports native Thermo centroiding and optional standalone OpenMS picking. This code path did not produce these mzML-based E. coli results.', '',
              '## Artifacts and reproduction', '',
              '- `PXD003327_comparison.pdf`: nine pages of coverage, ratio distributions, paired peptide/protein ratios, and parameter explanations.',
              '- `coverage.tsv`, `coverage_by_injection.tsv`, `benchmark.tsv`, `population_accuracy.tsv`: numerical plot/statistics tables.',
              '- `*.peptide_ratios.tsv`, `*.protein_ratios.tsv`: source-level exported observations including missingness/status.',
              '- `comparison.json`: exact software, input hashes, full settings, resolved search configurations, coverage and turnover.',
              '- `commands/`, `*.commands.diff`, `tasks.tsv`: verbatim executed task scripts, their pairwise diffs, source task paths/status and SHA-256 hashes.',
              '- `source_snapshot/`: today’s configuration, relevant modules, metadata, chemistry and analysis scripts.', '',
              'Run from the repository root:', '', '```bash',
              'MPLCONFIGDIR=/tmp/pxd003327-comparison-mpl /home/philipp/miniconda3/envs/silac/bin/python bin/compare_pxd003327_runs.py', '```', '',
              'The source output/command directories must remain available. Results are stored under Git-ignored `results/`; preserve this report directory separately if needed. The local accession label PXD003327 has a documented PXD003845 provenance ambiguity; results here refer specifically to the named local Chris_Ecoli files.', '']
    (out / 'changes_and_parameters.md').write_text('\n'.join(lines))


if __name__ == '__main__':
    main()
