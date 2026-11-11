"""核对领域原语的业务数值。"""
import unittest
from analysis_data.io import load
from analysis_data.domain import summarize
class DomainTests(unittest.TestCase):
    def test_balanced_and_post(self):
        rows=summarize(load('data/normal.json')[0]);r={x['site']:x['mean'] for x in rows}
        self.assertEqual(r['donor-a'],r['donor-c']);self.assertEqual(r['treated'],10.5)
        self.assertNotEqual(rows,summarize(load('data/changed-post-treatment.json')[0]))
