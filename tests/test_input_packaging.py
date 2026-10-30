import importlib.util
import json
from pathlib import Path
import sys
import tarfile
import tempfile
import unittest

CODE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(CODE / 'scripts'))
SPEC = importlib.util.spec_from_file_location('study_inputs', CODE / 'scripts/study_inputs.py')
inputs = importlib.util.module_from_spec(SPEC);SPEC.loader.exec_module(inputs)
from integrity import file_manifest, sha, verify_data

class ExistingInputPackagingTests(unittest.TestCase):
    def test_existing_upstream_inputs_are_repackaged_once(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            project, data, out = (root / p for p in ('project', 'data', 'package'))
            project.mkdir(); data.mkdir()
            excluded = set((inputs.UPSTREAM / 'remove_samples.txt').read_text().split())
            records = []
            for group in ('ss14', 'nichols'):
                for sample in (inputs.UPSTREAM / 'list' / ('tpa.'+group+'.sampleList.txt')).read_text().split():
                    if sample not in excluded:
                        records.append(dict(sample=sample, year='2000', lineage='SS14' if group=='ss14' else 'Nichols', country='Country'))
            assemblies = inputs.rows(inputs.UPSTREAM / 'config/complete_genome_manifest.tsv')
            for row in assemblies:
                records.append(dict(sample=row['sample_id'], year='2000', lineage='Nichols', country='Country'))
                dest = project / row['fasta'];dest.parent.mkdir(parents=True, exist_ok=True);dest.write_text('>reference\nACGT\n')
            metadata = data / 'Supplementary_Tables/TableS1_sample_metadata.tsv'
            inputs.write_tsv(metadata, records, ['sample', 'year', 'lineage', 'country'])
            (project / 'sample_dates.tsv').write_text(''.join(r['sample']+'\t2000\n' for r in records))
            for name in ('TPA_SS14.fasta', 'TPA_Nichols.fasta'):
                dest = project / 'refs' / name;dest.parent.mkdir(exist_ok=True);dest.write_text('>reference\nACGT\n')
            old_inputs = data / 'Upstream_Inputs';old_inputs.mkdir()
            (old_inputs / 'reference_proteins.fasta').write_text('>mrcA\nMA\n')
            (old_inputs / 'esm2_model.json').write_text(json.dumps(dict(model_files=[dict(path='weights.bin', sha256='a'*64)])))
            (old_inputs / 'README.md').write_text('Previous study input instructions\n')
            frozen = data / 'frozen_run'
            prepared = frozen / 'work/figures/data';prepared.mkdir(parents=True)
            (prepared / 'values.tsv').write_text('value\n1\n')
            (frozen / 'work/INPUT_MANIFEST.json').write_text(json.dumps(file_manifest(prepared)))
            (frozen / 'source.tsv').write_text('source\n1\n')
            (frozen / 'SOURCE_MANIFEST.json').write_text(json.dumps(dict(source=dict(path='source.tsv',sha256=sha(frozen/'source.tsv')))))
            original = file_manifest(data)
            (data / 'DATA_MANIFEST.json').write_text(json.dumps(original))
            inputs.package(project, data, out)
            status = json.loads((out/'PACKAGE_STATUS.json').read_text())
            self.assertFalse(status['upstream_database_versions_complete'])
            for name, record in original.items(): self.assertEqual(sha(data/name), record['sha256'])
            extract = root / 'extract';extract.mkdir()
            with tarfile.open(out/'TPA_pbp_evolution_data.tar.gz') as archive:
                names = archive.getnames();self.assertEqual(len(names), len(set(names)))
                archive.extractall(extract, filter='data')
            self.assertEqual(verify_data(extract/'Data_for_Zenodo')['state'], 'PASS')
            self.assertEqual(status['genomes'], 1861)

if __name__ == '__main__': unittest.main()
