"""Check the handoff from reference-context recomputation to figure inputs."""
import importlib.util
import json
from pathlib import Path
import sys
import tempfile
import unittest

SCRIPTS = Path(__file__).resolve().parents[1]/'scripts'
sys.path.insert(0, str(SCRIPTS))
SPEC = importlib.util.spec_from_file_location('esm2_input_contract', SCRIPTS/'esm2_inputs.py')
binding = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(binding)


class Esm2ContractTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.root = Path(self.temporary.name)
        self.inference = self.root/'inference'
        self.inference.mkdir()
        self.reference = self.root/'reference_proteins.fasta'
        self.reference.write_text('reference-context fixture\n')
        self.scores = self.inference/'esm2_all_substitutions.tsv'
        self.scores.write_text('score-table fixture\n')
        self.parent = dict(state='complete', exit_code=0, score_comparison='not_requested',
                           score_comparison_requested=False,
                           reference_fasta_sha256=binding.sha(self.reference))
        self.write(self.root/'RUN_STATUS.json', self.parent)
        self.write(self.root/'REFERENCE_AUDIT.json',
                   dict(state='PASS', protein_fasta_sha256=binding.sha(self.reference)))
        self.write(self.inference/'RUN_STATUS.json', dict(state='complete', exit_code=0))
        self.write(self.inference/'INPUT_CONTEXT.json',
                   dict(protein_fasta_sha256=binding.sha(self.reference)))
        self.write(self.inference/'MODEL_CONTEXT.json', {})
        self.write(self.inference/'ESM2_COVERAGE.json', {})

    def tearDown(self):
        self.temporary.cleanup()

    @staticmethod
    def write(path, value):
        path.write_text(json.dumps(value)+'\n')

    def compared(self):
        self.parent.update(score_comparison='PASS', score_comparison_requested=True)
        self.write(self.root/'RUN_STATUS.json', self.parent)
        self.write(self.root/'SCORE_COMPARISON.json',
                   dict(state='PASS', new_scores_sha256=binding.sha(self.scores)))

    def test_independent_recomputation_needs_no_previous_scores(self):
        reference, audit, comparison = binding.completed_scan(self.inference)
        self.assertEqual(reference, self.reference)
        self.assertEqual(comparison, 'not_requested')
        self.assertTrue(audit.is_file())
        self.assertFalse((self.root/'SCORE_COMPARISON.json').exists())

    def test_requested_comparison_must_pass_and_match_scores(self):
        self.compared()
        self.assertEqual(binding.completed_scan(self.inference)[2], 'PASS')
        self.scores.write_text('changed score-table fixture\n')
        with self.assertRaisesRegex(ValueError, 'does not match'):
            binding.completed_scan(self.inference)

    def test_failed_comparison_cannot_be_adopted(self):
        self.parent.update(score_comparison='REVIEW_REQUIRED', score_comparison_requested=True)
        self.write(self.root/'RUN_STATUS.json', self.parent)
        with self.assertRaisesRegex(ValueError, 'did not pass'):
            binding.completed_scan(self.inference)

    def test_requested_comparison_cannot_be_recorded_as_not_requested(self):
        self.parent['score_comparison_requested'] = True
        self.write(self.root/'RUN_STATUS.json', self.parent)
        with self.assertRaisesRegex(ValueError, 'inconsistent'):
            binding.completed_scan(self.inference)

    def test_changed_reference_is_rejected(self):
        self.reference.write_text('changed reference-context fixture\n')
        with self.assertRaisesRegex(ValueError, 'translation audit'):
            binding.completed_scan(self.inference)

    def test_existing_compared_outputs_remain_compatible(self):
        self.compared()
        self.parent.pop('score_comparison_requested')
        self.write(self.root/'RUN_STATUS.json', self.parent)
        self.assertEqual(binding.completed_scan(self.inference)[2], 'PASS')

    def test_low_level_output_requires_annotated_reference_provenance(self):
        (self.root/'REFERENCE_AUDIT.json').unlink()
        with self.assertRaisesRegex(ValueError, 'annotated-reference recomputation'):
            binding.completed_scan(self.inference)


if __name__ == '__main__':
    unittest.main()
