import importlib.util
from pathlib import Path
import unittest

SCRIPT = Path(__file__).resolve().parents[1] / 'scripts/external_resources.py'
SPEC = importlib.util.spec_from_file_location('external_resources', SCRIPT)
resources = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(resources)

class ExternalResourceTests(unittest.TestCase):
    def record(self, resource_id):
        return dict(resource_id=resource_id, name=resource_id, identity_complete=True,
                    retrieval_complete=True, source_path='/private/database',
                    source_description='Reference sequences listed in the database recipe',
                    recipe='Build from the pinned reference accessions',
                    files=[dict(path='database.bin', size=12, sha256='a'*64)])

    def test_exact_identities_and_access_routes_are_required(self):
        value = dict(schema_version=1, resources=[self.record(k) for k in resources.REQUIRED])
        result = resources.public_provenance(value)
        self.assertTrue(result['upstream_database_versions_complete'])
        self.assertNotIn('source_path', result['resources'][0])
        value['resources'][0].pop('recipe')
        self.assertFalse(resources.public_provenance(value)['upstream_database_versions_complete'])

    def test_manual_boolean_cannot_replace_file_checksums(self):
        value = self.record('strainseeker');value['files'] = []
        result = resources.public_provenance(dict(schema_version=1, resources=[value]))
        self.assertFalse(result['upstream_database_versions_complete'])
        self.assertFalse(result['resources'][0]['identity_complete'])

    def test_private_or_escaping_file_names_are_rejected(self):
        for name in ('/private/file', '../file', 'folder\\file'):
            value = self.record('strainseeker');value['files'][0]['path'] = name
            with self.assertRaises(ValueError):
                resources.public_provenance(dict(schema_version=1, resources=[value]))

    def test_missing_databases_remain_explicit(self):
        result = resources.public_provenance(dict(schema_version=1, resources=[]))
        self.assertEqual(set(result['unresolved_resources']), set(resources.REQUIRED))

if __name__ == '__main__': unittest.main()
