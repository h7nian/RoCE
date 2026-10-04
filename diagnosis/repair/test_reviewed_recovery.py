#!/usr/bin/env python3
import copy
from pathlib import Path
import sys
import tempfile
import unittest
sys.dont_write_bytecode = True
import recover_reviewed_repeats as recovery


class ReviewedRecoveryTests(unittest.TestCase):
    def test_combined_cache_io_errors_require_every_component_to_be_io(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            failure = root / 'score_derivative_FAILED.txt'
            evidence = root / 'peer.err'
            evidence.write_text('cannot rename cache: Input/output error\n')
            prefix = 'parallel_lapply: worker(s) for item(s) mu1, mu0 failed: '
            connection = 'Error in file(file, mode) : cannot open the connection'
            atomic = 'Error : Could not atomically store the shared nuisance fit.'
            for errors in ((connection, atomic), (atomic, connection), (atomic,)):
                failure.write_text(prefix + ' | '.join(errors) + '\n')
                item = dict(failure_kind='reviewed_R_file_io',
                    reviewed_error_files={failure.name: recovery.pilot.digest(failure)},
                    io_evidence=[dict(path=str(evidence), sha256=recovery.pilot.digest(evidence))])
                recovery.validate_r_failure_review(root, item)
            for error in ('Optimization failed: KKT criterion not satisfied',
                          'Error : Could not atomically store an unknown result.', ''):
                failure.write_text(prefix + connection + ' | ' + error + '\n')
                item['reviewed_error_files'][failure.name] = recovery.pilot.digest(failure)
                with self.assertRaisesRegex(ValueError, 'not a reviewed'):
                    recovery.validate_r_failure_review(root, item)

    def test_only_hash_pinned_file_errors_with_io_evidence_are_eligible(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            failure = root / 'score_derivative_FAILED.txt'
            failure.write_text('parallel_lapply: worker(s) for item(s) mu1, mu0 failed: Error in file(file, mode) : cannot open the connection\n')
            evidence = root / 'peer.err'
            evidence.write_text('cannot open file: Input/output error\n')
            item = dict(failure_kind='reviewed_R_file_io',
                        reviewed_error_files={failure.name: recovery.pilot.digest(failure)},
                        io_evidence=[dict(path=str(evidence), sha256=recovery.pilot.digest(evidence))])
            recovery.validate_r_failure_review(root, item)
            for changed in [dict(item, failure_kind='node_startup_or_io'),
                            dict(item, reviewed_error_files={}), dict(item, io_evidence=[])]:
                with self.assertRaises(ValueError):
                    recovery.validate_r_failure_review(root, changed)
            failure.write_text('Optimization failed: KKT criterion not satisfied\n')
            changed = copy.deepcopy(item)
            changed['reviewed_error_files'][failure.name] = recovery.pilot.digest(failure)
            with self.assertRaisesRegex(ValueError, 'not a reviewed'):
                recovery.validate_r_failure_review(root, changed)
            failure.write_text('Error in file(file, mode) : cannot open the connection\n')
            with self.assertRaisesRegex(ValueError, 'changed'):
                recovery.validate_r_failure_review(root, item)
            changed['reviewed_error_files'][failure.name] = recovery.pilot.digest(failure)
            evidence.write_text('ordinary diagnostic without an I/O failure\n')
            with self.assertRaisesRegex(ValueError, 'evidence changed'):
                recovery.validate_r_failure_review(root, changed)
            changed['io_evidence'][0]['sha256'] = recovery.pilot.digest(evidence)
            with self.assertRaises(ValueError):
                recovery.validate_r_failure_review(root, changed)

    def test_pre_r_recovery_still_refuses_undeclared_r_errors(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            recovery.validate_r_failure_review(root, dict(failure_kind='pre_R_startup_reviewed'))
            with self.assertRaisesRegex(ValueError, 'missing'):
                recovery.validate_r_failure_review(root, dict(failure_kind='reviewed_R_file_io'))
            (root/'recipe_FAILED.txt').write_text('convergence failure\n')
            with self.assertRaisesRegex(ValueError, 'separate review'):
                recovery.validate_r_failure_review(root, dict(failure_kind='pre_R_startup_reviewed'))


if __name__ == '__main__':
    unittest.main()
