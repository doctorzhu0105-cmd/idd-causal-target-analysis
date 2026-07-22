import sys
import unittest
from pathlib import Path


sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "py"))

from idd_rank_targets import (  # noqa: E402
    choose_localization,
    delivery_fraction,
    directional_delivery_fraction,
    is_credible_modulator,
    modulator_direction_aligned,
    pass_smr,
    recommended_local_payload,
    select_distinct_loci,
)


class RankingRulesTests(unittest.TestCase):
    def test_smr_requires_fdr_heidi_and_ten_variants(self):
        base = {
            "status": "estimable",
            "p_smr_fdr": "0.01",
            "p_heidi": "0.2",
            "n_snp_heidi": "10",
        }
        self.assertTrue(pass_smr(base))
        self.assertFalse(pass_smr({**base, "p_smr_fdr": "0.06"}))
        self.assertFalse(pass_smr({**base, "p_heidi": "0.005"}))
        self.assertFalse(pass_smr({**base, "n_snp_heidi": "9"}))

    def test_not_estimable_localization_does_not_pass(self):
        best, passed = choose_localization(
            [{"gene_symbol": "GENE", "status": "gene_not_detected_in_feature_space"}],
            0.05,
        )
        self.assertIsNone(best)
        self.assertFalse(passed)

    def test_localization_requires_donor_and_cell_support(self):
        rows = [
            {
                "status": "localized",
                "broad_cell_type": "nucleus_pulposus",
                "cell_fraction_expressing": "0.20",
                "mean_donor_expression": "0.1",
                "donor_detection_fraction": "1",
                "donors": "1",
                "cells": "100",
            },
            {
                "status": "localized",
                "broad_cell_type": "annulus_fibrosus",
                "cell_fraction_expressing": "0.10",
                "mean_donor_expression": "0.2",
                "donor_detection_fraction": "0.75",
                "donors": "4",
                "cells": "100",
            },
        ]
        best, passed = choose_localization(rows, 0.05)
        self.assertTrue(passed)
        self.assertEqual(best["broad_cell_type"], "annulus_fibrosus")

    def test_delivery_uses_half_weight_for_moderate_feasibility(self):
        self.assertEqual(delivery_fraction("high_for_local_small_molecule"), 1.0)
        self.assertEqual(delivery_fraction("moderate_for_local_siRNA"), 0.5)
        self.assertEqual(delivery_fraction("not_assessed"), 0.0)

    def test_directionless_dgidb_record_is_not_a_credible_modulator(self):
        self.assertFalse(
            is_credible_modulator(
                {"source": "DGIdb", "interaction_type": "", "directionality": ""}
            )
        )
        self.assertTrue(
            is_credible_modulator(
                {"source": "DGIdb", "interaction_type": "inhibitor", "directionality": "INHIBITORY"}
            )
        )

    def test_directional_delivery_does_not_equate_inhibitor_with_expression_increase(self):
        surface_row = {
            "preferred_local_payload_options": "small_molecule;siRNA",
            "extracellular_or_surface_location": "True",
        }
        intracellular_row = {
            "preferred_local_payload_options": "small_molecule;siRNA",
            "extracellular_or_surface_location": "False",
        }
        self.assertEqual(
            directional_delivery_fraction(surface_row, "increase_or_replace_expression", False),
            1.0,
        )
        self.assertEqual(
            directional_delivery_fraction(intracellular_row, "increase_or_replace_expression", False),
            0.5,
        )
        self.assertEqual(
            directional_delivery_fraction(intracellular_row, "inhibit_or_reduce_expression", True),
            1.0,
        )

    def test_modulator_alignment_follows_mr_intervention_direction(self):
        inhibitor = {"interaction_type": "inhibitor", "directionality": "INHIBITORY"}
        agonist = {"interaction_type": "agonist", "directionality": "ACTIVATING"}
        self.assertTrue(modulator_direction_aligned(inhibitor, "inhibit_or_reduce_expression"))
        self.assertFalse(modulator_direction_aligned(inhibitor, "increase_or_replace_expression"))
        self.assertTrue(modulator_direction_aligned(agonist, "increase_or_replace_expression"))

    def test_payload_recommendation_respects_intervention_direction(self):
        surface_row = {"extracellular_or_surface_location": "True"}
        intracellular_row = {"extracellular_or_surface_location": "False"}
        inhibitor = {"canonical_smiles": "CCO"}
        self.assertIn(
            "recombinant_protein",
            recommended_local_payload(
                surface_row, "increase_or_replace_expression", inhibitor, False
            ),
        )
        self.assertEqual(
            recommended_local_payload(
                intracellular_row, "increase_or_replace_expression", inhibitor, False
            ),
            "mrna_nanoparticle_in_hydrogel",
        )
        self.assertIn(
            "small_molecule",
            recommended_local_payload(
                surface_row, "inhibit_or_reduce_expression", inhibitor, True
            ),
        )

    def test_priority_selection_keeps_one_gene_per_locus(self):
        rows = [
            {"gene_symbol": "A", "chr_b37": "1", "start_b37": 100, "end_b37": 200},
            {"gene_symbol": "B", "chr_b37": "1", "start_b37": 300, "end_b37": 400},
            {"gene_symbol": "C", "chr_b37": "2", "start_b37": 100, "end_b37": 200},
        ]
        selected, alternatives = select_distinct_loci(rows, 2, 1000)
        self.assertEqual([row["gene_symbol"] for row in selected], ["A", "C"])
        self.assertEqual(alternatives[0]["alternative_to"], "A")


if __name__ == "__main__":
    unittest.main()
