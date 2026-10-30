import importlib
import io
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

from Bio import Phylo

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))
common = importlib.import_module('analysis.common')
sites = importlib.import_module('analysis.selection_topologies.fubar_sites')


class SharedAnalysisTests(unittest.TestCase):
    def test_modules_use_the_same_shared_implementation(self):
        temporal = importlib.import_module('analysis.clustered_dates.temporal')
        core = importlib.import_module('analysis.clustered_dates.cluster_core')
        ensemble = importlib.import_module('analysis.selection_topologies.ensemble')
        self.assertIs(temporal.prune, common.prune)
        self.assertIs(ensemble.prune, common.prune)
        self.assertIs(core.sha, common.sha)
        self.assertIs(ensemble.json_sites, sites.json_sites)

    def test_pruning_preserves_retained_distances_and_original_tree(self):
        tree = Phylo.read(io.StringIO('((a:1,b:2):3,c:4);'), 'newick')
        distance = tree.distance('a', 'c')
        retained = common.prune(tree, {'a', 'c'})
        self.assertEqual(common.validate_tree(retained), {'a', 'c'})
        self.assertEqual(retained.distance('a', 'c'), distance)
        self.assertEqual(common.validate_tree(tree), {'a', 'b', 'c'})
        with self.assertRaisesRegex(ValueError, 'Samples absent'):
            common.prune(tree, {'absent'})

    def test_fubar_partition_positions_map_to_exact_sites(self):
        record = {
            'input': {'number of sites': 2, 'number of sequences': 3},
            'MLE': {'headers': [['Prob[alpha<beta]', 'probability']],
                    'content': {'0': [[0.2]], '1': [[0.95]]}},
            'data partitions': {'0': {'coverage': [[1]]},
                                '1': {'coverage': [[0]]}},
        }
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / 'sites.json'
            path.write_text(json.dumps(record))
            self.assertEqual(sites.json_sites(path, 'fubar', 2, 3), {1: 0.2, 0: 0.95})
            record['MLE']['content']['1'] = [[1.01]]
            path.write_text(json.dumps(record))
            with self.assertRaisesRegex(ValueError, 'Invalid HyPhy probability'):
                sites.json_sites(path, 'fubar', 2, 3)

    def test_resumable_context_rejects_changed_settings(self):
        with tempfile.TemporaryDirectory() as directory:
            out = Path(directory)
            common.context(out, {'setting': 1})
            common.context(out, {'setting': 1})
            with self.assertRaisesRegex(ValueError, 'settings'):
                common.context(out, {'setting': 2})
            self.assertEqual(json.loads((out / 'CONTEXT.json').read_text()), {'setting': 1})

    def test_direct_and_unified_entry_points_load(self):
        env = dict(os.environ, PYTHONDONTWRITEBYTECODE='1')
        for module, script in [('clustered-dates', 'clustered_dates/temporal.py'),
                               ('selection-topologies', 'selection_topologies/ensemble.py')]:
            commands = [
                [sys.executable, str(ROOT / 'scripts/analysis' / script), '--help'],
                ['bash', str(ROOT / 'master.sh'), 'analysis', module, '--', '--help'],
            ]
            for command in commands:
                result = subprocess.run(command, cwd=ROOT, env=env, text=True,
                                        stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                self.assertEqual(result.returncode, 0, result.stderr)
                self.assertIn('usage:', result.stdout)


if __name__ == '__main__':
    unittest.main()
