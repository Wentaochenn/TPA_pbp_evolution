#!/usr/bin/env python3
"""Read-only audit of direct protein numbering between two local references."""
import argparse
import hashlib
from pathlib import Path

from common import sha256, table, write_json

# Corrected original project annotations, 1-based inclusive chromosome intervals.
COORDINATES = {
    'mrcA': {'SS14': (772313, 774967, '-', 884),
             'Nichols': (772433, 775087, '-', 884)},
    'pbp1': {'SS14': (535746, 537620, '+', 624),
             'Nichols': (535881, 537755, '+', 624)},
    'pbp2': {'SS14': (823707, 825665, '+', 652),
             'Nichols': (823827, 825785, '+', 652)},
    'tp47': {'SS14': (623741, 625045, '-', 434),
             'Nichols': (623858, 625162, '-', 434)},
}
ACCESSIONS = {'SS14': 'NC_021508.1', 'Nichols': 'NC_021490.2'}


def alignment_audit(left, right):
    """Global unit edit costs; retain whether ANY optimal path has a gap.

    End gaps are penalized. A gap in an optimal path blocks direct numbering;
    it is not itself a claim that a biological indel has occurred.
    O(len(left)*len(right)) time, O(len(right)) memory.
    """
    if not left or not right or set(left + right) - set('ACDEFGHIKLMNPQRSTVWY'):
        raise ValueError('Protein audit requires nonempty unambiguous proteins without stops')
    previous = list(range(len(right) + 1))
    previous_gap = [False] + [True] * len(right)
    for i, a in enumerate(left, 1):
        current, current_gap = [i], [True]
        for j, b in enumerate(right, 1):
            diagonal = previous[j-1] + (a != b)
            deletion, insertion = previous[j] + 1, current[j-1] + 1
            best = min(diagonal, deletion, insertion)
            current.append(best)
            current_gap.append((diagonal == best and previous_gap[j-1]) or
                               deletion == best or insertion == best)
        previous, previous_gap = current, current_gap
    nogap = len(left) == len(right) and not previous_gap[-1]
    return {'global_edit_distance': previous[-1],
            'any_optimal_alignment_has_gap': previous_gap[-1],
            'same_numbering_verified': nogap,
            'same_index_identity': sum(a == b for a, b in zip(left, right)) / len(left)
                if len(left) == len(right) else None,
            'alignment_method': 'global_unit_edit_cost_match0_substitution1_insertion1_deletion1; all_optimal_paths_gap_audit'}


def reference_audit(ss14, nichols, outdir, coordinates=None):
    """Write a fresh audit and fail closed on uncertain direct numbering."""
    from prepare_existing_sequences import records, revcomp, protein
    coordinates = COORDINATES if coordinates is None else coordinates
    if not isinstance(coordinates, dict) or not coordinates:
        raise ValueError('Reference audit requires a nonempty gene-coordinate scope')
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=False)
    paths = {'SS14': Path(ss14).resolve(), 'Nichols': Path(nichols).resolve()}
    report = {'status': 'running', 'scope': sorted(coordinates), 'references': {}, 'genes': {},
              'interpretation': 'Checks local reference CDS and direct protein numbering only; no sample indel or annotation-origin validation'}
    proteins, rows = {}, []
    try:
        for label, path in paths.items():
            before = sha256(path)
            items = list(records(path))
            if len(items) != 1:
                raise ValueError(label + ': expected exactly one chromosome reference')
            identifier, genome = items[0]
            if ACCESSIONS[label] not in identifier.replace('|', ' ').split():
                raise ValueError(label + ': expected accession ' + ACCESSIONS[label] + '; found ' + identifier)
            if set(genome) - set('ACGTRYSWKMBDHVN'):
                raise ValueError(label + ': reference must be gap-free nucleotide chromosome')
            report['references'][label] = {'path': str(path), 'id': identifier,
                                           'length': len(genome), 'sha256': before}
            proteins[label] = {}
            for gene, entries in coordinates.items():
                start, end, strand, expected = entries[label]
                if not 1 <= start <= end <= len(genome) or strand not in ('+', '-'):
                    raise ValueError(label + '/' + gene + ': invalid chromosome interval or strand')
                cds = genome[start-1:end]
                if strand == '-':
                    cds = revcomp(cds)
                if set(cds) - set('ACGT'):
                    raise ValueError(label + '/' + gene + ': ambiguous reference CDS')
                aa = protein(cds)
                if (len(cds) != (expected+1)*3 or len(aa) != expected+1 or
                        not aa.endswith('*') or '*' in aa[:-1]):
                    raise ValueError(label + '/' + gene + ': CDS length/frame/terminal-stop check failed')
                proteins[label][gene] = aa[:-1]
                report['genes'].setdefault(gene, {})[label] = {
                    'coordinates': list(entries[label]), 'cds_length': len(cds), 'protein_length': len(aa)-1,
                    'cds_sha256': hashlib.sha256(cds.encode('ascii')).hexdigest(),
                    'protein_sha256': hashlib.sha256(aa[:-1].encode('ascii')).hexdigest()}
            if sha256(path) != before:
                raise ValueError('Reference changed while reading: ' + str(path))
        for gene in coordinates:
            left, right = proteins['Nichols'][gene], proteins['SS14'][gene]
            result = alignment_audit(left, right)
            # A high concordance guard catches wrong gene extraction. This is
            # a coordinate QC threshold, not an evolutionary hypothesis test.
            result['minimum_same_index_identity'] = 0.95
            result['pass'] = result['same_numbering_verified'] and result['same_index_identity'] >= 0.95
            report['genes'][gene]['comparison'] = result
            if result['same_numbering_verified']:
                rows.extend({'gene': gene, 'position': i, 'Nichols_aa': a, 'SS14_aa': b}
                            for i, (a, b) in enumerate(zip(left, right), 1) if a != b)
        table(outdir/'reference_residue_differences.tsv', rows, ['gene','position','Nichols_aa','SS14_aa'])
        if any(not entry['comparison']['pass'] for entry in report['genes'].values()):
            raise ValueError('Direct Nichols/SS14 protein numbering not verified; inspect reference audit; no automatic remapping')
        report['status'] = 'complete'
        return report, proteins
    except Exception as exc:
        report.update(status='failed', error=str(exc))
        raise
    finally:
        write_json(outdir/'REFERENCE_NUMBERING.json', report)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--ss14', default='refs/TPA_SS14.fasta')
    parser.add_argument('--nichols', default='refs/TPA_Nichols.fasta')
    parser.add_argument('--outdir', required=True)
    args = parser.parse_args()
    report, _ = reference_audit(args.ss14, args.nichols, args.outdir)
    print('gene\tNichols_aa\tSS14_aa\tedit_distance\tany_optimal_gap\tstatus')
    for gene, entry in report['genes'].items():
        check = entry['comparison']
        print('\t'.join(map(str, (gene, entry['Nichols']['protein_length'],
              entry['SS14']['protein_length'], check['global_edit_distance'],
              check['any_optimal_alignment_has_gap'], 'PASS' if check['pass'] else 'REVIEW'))))
    print('PASS: local reference numbering audit; report=' + str(Path(args.outdir)/'REFERENCE_NUMBERING.json'))


if __name__ == '__main__':
    main()
