import sys
import unittest
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "04_consensus_screening"))
from consensus_prediction import consensus
from candidate_selection import select


class ConsensusTests(unittest.TestCase):
    def setUp(self):
        self.rows = [dict(peptide_id="x", sequence="ACDEFGH", model=m, prediction="1") for m in ["A", "B"]]
        self.mapping = {"A": "Family_A", "B": "Family_B"}

    def test_votes(self):
        self.assertEqual(consensus(self.rows, ["A", "B"], 2, self.mapping)[0]["votes"], 2)

    def test_duplicates_do_not_inflate(self):
        with self.assertRaises(ValueError):
            consensus(self.rows + self.rows[:1], ["A", "B"], 2, self.mapping)

    def test_missing_is_not_negative(self):
        with self.assertRaises(ValueError):
            consensus(self.rows[:1], ["A", "B"], 2, self.mapping)

    def test_review_required(self):
        rows = [{k: str(v) for k, v in r.items()} for r in consensus(self.rows, ["A", "B"], 2, self.mapping)]
        reviewed = select(rows, [])
        self.assertEqual(reviewed[0]["selected"], 0)

    def test_conflicting_sequence(self):
        self.rows[1]["sequence"] = "AAAAAAA"
        with self.assertRaises(ValueError):
            consensus(self.rows, ["A", "B"], 2, self.mapping)

    def test_same_family_counts_once(self):
        mapping = {"A": "CAMP", "B": "CAMP"}
        result = consensus(self.rows, ["A", "B"], 1, mapping, {"CAMP": "all"})[0]
        self.assertEqual((result["n_models"], result["n_families"], result["votes"]), (2, 1, 1))

    def test_family_rule_not_invented(self):
        with self.assertRaises(ValueError):
            consensus(self.rows, ["A", "B"], 1, {"A": "CAMP", "B": "CAMP"})

    def test_family_threshold_and_majority_tie(self):
        mapping = {"A": "CAMP", "B": "CAMP"}
        with self.assertRaises(ValueError):
            consensus(self.rows, ["A", "B"], 2, mapping, {"CAMP": "all"})
        self.rows[1]["prediction"] = "0"
        result = consensus(self.rows, ["A", "B"], 1, mapping, {"CAMP": "majority"})[0]
        self.assertEqual(result["votes"], 0)

    def test_seven_outputs_five_families(self):
        models = ["AI4AMP", "AMPScanner", "macrel", "AMPlify", "CAMP_Bert", "CAMP_LSTM", "CAMP_Attention"]
        mapping = {m: "CAMP" if m.startswith("CAMP_") else m for m in models}
        rows = [dict(peptide_id="x", sequence="ACDEFGH", model=m, prediction="1") for m in models]
        result = consensus(rows, models, 5, mapping, {"CAMP": "all"})[0]
        self.assertEqual((result["n_models"], result["n_families"], result["votes"]), (7, 5, 5))


if __name__ == "__main__":
    unittest.main()
