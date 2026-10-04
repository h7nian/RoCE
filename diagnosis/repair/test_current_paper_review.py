#!/usr/bin/env python3
import sys
import unittest
sys.dont_write_bytecode = True
import review_current_paper_mc200 as review


class CurrentPaperReviewTests(unittest.TestCase):
    def rows(self):
        return [dict(panel="main", role=role, config="C2", K="8", p="100", rho="1.5",
                     sim_id=str(seed), deviation_mechanism="treated_arm", n_deviated_sites="1")
                for seed in range(1,201) for role in ("method","baseline")]

    def test_pair_index_rejects_duplicate_missing_or_wrong_seeds(self):
        rows = self.rows()
        self.assertEqual(len(review.index_pairs(rows,"main",100)),200)
        with self.assertRaisesRegex(ValueError,"Duplicate"):
            review.index_pairs(rows+[rows[0]],"main",100)
        with self.assertRaisesRegex(ValueError,"one method and baseline"):
            review.index_pairs(rows[:-1],"main",100)
        with self.assertRaisesRegex(ValueError,"exactly200"):
            review.index_pairs(rows[:-2],"main",100)
        for row in rows[-2:]:row["sim_id"]="401"
        with self.assertRaisesRegex(ValueError,"seeds1"):
            review.index_pairs(rows,"main",100)

    def test_pairing_checks_data_references_and_all_requested_baselines(self):
        value = dict(estimate=".25", se=".02", truth=".24")
        method = {name:dict(value) for name in ["one_round_crossfit_ate", "one_round_crossfit_ate_separate_arms",
            "one_round_crossfit_ate_joint_tate", "target_only_ate", "target_anchor_ate"]}
        baseline = {name+"_ate":dict(value) for name in review.pilot.BASELINE_METHODS}
        baseline["target_only_ate"] = dict(value)
        self.assertEqual(len(review.combine_pair({},method,baseline,"abc","abc")),9)
        with self.assertRaisesRegex(ValueError,"hashes"):
            review.combine_pair({},method,baseline,"abc","def")
        missing=dict(baseline);del missing["pooled_dr_ate"]
        with self.assertRaisesRegex(ValueError,"missing"):
            review.combine_pair({},method,missing,"abc","abc")
        baseline["target_only_ate"]["se"]=".03"
        with self.assertRaisesRegex(ValueError,"reference"):
            review.combine_pair({},method,baseline,"abc","abc")

    def test_primary_only_records_do_not_invent_external_baselines(self):
        value = dict(estimate=".25", se=".02", truth=".24")
        methods = {name:dict(value) for name in ["one_round_crossfit_ate", "one_round_crossfit_ate_separate_arms",
            "one_round_crossfit_ate_joint_tate", "target_only_ate", "target_anchor_ate"]}
        rows=review.primary_records({},methods,"abc")
        self.assertEqual(len(rows),5)
        self.assertFalse(any(row["method"]=="pooled_dr_ate" for row in rows))
        methods["target_anchor_ate"]["truth"]=".3"
        with self.assertRaisesRegex(ValueError,"truths differ"):
            review.primary_records({},methods,"abc")


if __name__ == "__main__":
    unittest.main()
