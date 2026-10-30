import importlib.util
from pathlib import Path
import unittest


MODULE_PATH = Path(__file__).resolve().parents[1] / 'scripts/supplementary_metadata.py'
SPEC = importlib.util.spec_from_file_location('supplementary_metadata', MODULE_PATH)
module = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(module)


class AminoAcidStateTests(unittest.TestCase):
    def setUp(self):
        self.indicators = [('A506T', 'T'), ('A506V', 'V'), ('A506M', 'M')]

    def state(self, values):
        row = dict(zip((field for field, _ in self.indicators), values), sample='one')
        return module.amino_acid_state(row, self.indicators, 'A')

    def test_joint_alternatives_resolve_actual_residue(self):
        self.assertEqual(self.state(['0', '1', '0']), 'V')
        self.assertEqual(self.state(['1', '0', '0']), 'T')
        self.assertEqual(self.state(['0', '0', '0']), 'A')

    def test_missing_or_partial_calls_are_uncallable(self):
        self.assertEqual(self.state(['NA', 'NA', 'NA']), 'uncallable')
        self.assertEqual(self.state(['NA', '0', '0']), 'uncallable')
        self.assertEqual(self.state(['1', 'NA', '0']), 'uncallable')

    def test_mutually_exclusive_alternatives_rejected(self):
        with self.assertRaisesRegex(ValueError, 'Mutually exclusive'):
            self.state(['1', '1', '0'])

    def test_non_ternary_state_rejected(self):
        with self.assertRaisesRegex(ValueError, 'Invalid ternary'):
            self.state(['0', '0.9', '0'])

    def test_stop_codon_preserved(self):
        row = {'sample': 'one', 'E184Ter': '1'}
        self.assertEqual(module.amino_acid_state(row, [('E184Ter', '*')], 'E'), '*')


class SchemaTests(unittest.TestCase):
    @staticmethod
    def fixture(protein_type='NA'):
        metadata = {field: '' for field in module.METADATA_COLUMNS}
        metadata.update(sample='one', mrcA_protein_type=protein_type)
        reverse = {one: three for three, one in module.AA_ONE.items()}
        matrix = {'sample': 'one'}
        for field, alternatives in module.EXPECTED_ALTERNATIVES.items():
            gene, site = field.split(':')
            for alternative in alternatives:
                indicator = '%s:p.%s%s%s' % (gene, reverse[site[0]], site[1:], reverse[alternative])
                matrix[indicator] = '0'
        return metadata, matrix

    def test_duplicate_samples_rejected(self):
        with self.assertRaisesRegex(ValueError, 'duplicate sample'):
            module.keyed_rows([{'sample': 'same'}, {'sample': 'same'}], 'test')

    def test_sample_mismatch_rejected(self):
        with self.assertRaisesRegex(ValueError, 'sample sets differ'):
            module.build_table([{'sample': 'one'}], [{'sample': 'other'}], ['sample'])

    def test_missing_audited_site_rejected(self):
        with self.assertRaisesRegex(ValueError, 'Audited codon set differs'):
            module.indicator_groups(['sample', 'mrcA:p.Ala506Thr'])

    def test_missing_alternative_at_present_site_rejected(self):
        _, matrix = self.fixture()
        del matrix['mrcA:p.Ala506Met']
        with self.assertRaisesRegex(ValueError, 'Audited alternatives differ at mrcA:A506'):
            module.indicator_groups(list(matrix))

    def test_untyped_is_distinct_from_site_uncallable(self):
        metadata, matrix = self.fixture()
        columns, rows, spec = module.build_table([metadata], [matrix], list(matrix))
        self.assertEqual(len(columns), 54)
        self.assertEqual(rows[0]['mrcA_protein_type'], 'untyped')
        self.assertEqual(rows[0]['mrcA:A506'], 'A')
        self.assertEqual(spec['protein_types_untyped'], 1)
        self.assertFalse(set(columns) & set(module.EXCLUDED_COLUMNS))

    def test_one_uncallable_site_does_not_remove_other_states(self):
        metadata, matrix = self.fixture('WT')
        matrix['mrcA:p.Gly708Ser'] = 'NA'
        matrix['mrcA:p.Ala506Thr'] = '1'
        _, rows, spec = module.build_table([metadata], [matrix], list(matrix))
        self.assertEqual(rows[0]['mrcA:G708'], 'uncallable')
        self.assertEqual(rows[0]['mrcA:A506'], 'T')
        self.assertEqual(rows[0]['pbp1:P564'], 'P')
        self.assertEqual(spec['protein_types_assigned'], 1)


if __name__ == '__main__':
    unittest.main()
