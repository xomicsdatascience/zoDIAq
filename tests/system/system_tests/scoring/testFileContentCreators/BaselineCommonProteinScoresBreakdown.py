import pandas as pd
from .MaccScoresBreakdown import MaccScoresBreakdown
from abc import abstractmethod


class BaselineCommonProteinScoresBreakdown(MaccScoresBreakdown):
    def _create_input_template_for_scoring_module(self):
        """
        FDR not being calculated, so there's only one rank
        """
        self.inputDf = self._create_input_for_quantifiable_proteins()
        self._add_macc_rank_values_to_input_template_and_remove_rank_label()

    def _create_generic_peptide_line(self):
        genericMarker = "X"
        return pd.DataFrame(
            [[f"peptide{genericMarker}", f"1/protein{genericMarker}", 1, 0, 0, 100.0]],
            columns=["peptide", "protein", "zLIB", "isDecoy", "rank", "ionCount"],
        )

    @abstractmethod
    def _create_input_for_quantifiable_proteins(self):
        pass

    def _make_expected_spectral_output(self, scoredDf):
        return scoredDf.copy()

    def _make_expected_peptide_output(self, scoredDf):
        expectedPeptideOutput = scoredDf.copy()
        expectedPeptideOutput["peptideFDR"] = [0] * len(expectedPeptideOutput.index)
        return expectedPeptideOutput

    def _make_expected_protein_output(self, peptideDf):
        expectedProteinOutput = peptideDf.copy()
        expectedProteinOutput["leadingProtein"] = expectedProteinOutput["protein"]
        expectedProteinOutput["proteinCosine"] = expectedProteinOutput["cosine"]
        expectedProteinOutput["leadingProteinFDR"] = expectedProteinOutput["peptideFDR"]
        return expectedProteinOutput
