import pandas as pd
from .BaselineCommonProteinScoresBreakdown import BaselineCommonProteinScoresBreakdown

"""
NOTE: this is specific to the maxlfq method. Samples with 1 or 0 peptides 
    shared with other samples are considered 0.
 
protein quantification profile:
        peptide1    peptide2    peptide3
sample1 X           X           -
sample2 X           X           -
sample3 X           -           X

protein quantification:
        protein
sample1 X
sample2 X
sample3 0
"""


class OneMatchSample1And2Breakdown(BaselineCommonProteinScoresBreakdown):
    def _create_input_for_quantifiable_proteins(self):
        df = pd.DataFrame(
            [
                ["peptide1", "1/protein", 1, 0, 0, 100.0],
                ["peptide2", "1/protein", 1, 0, 0, 100.0],
            ],
            columns=["peptide", "protein", "zLIB", "isDecoy", "rank", "ionCount"],
        )
        return pd.concat([df, self._create_generic_peptide_line()])

    def _make_expected_protein_output(self, peptideDf):
        proteinDf = super()._make_expected_protein_output(peptideDf)
        proteinDf["uniquePeptide"] = [0, 0, 1]
        return proteinDf


class OneMatchSample3Breakdown(BaselineCommonProteinScoresBreakdown):
    def _create_input_for_quantifiable_proteins(self):
        df = pd.DataFrame(
            [
                ["peptide1", "1/protein", 1, 0, 0, 100.0],
                ["peptide3", "1/protein", 1, 0, 0, 100.0],
            ],
            columns=["peptide", "protein", "zLIB", "isDecoy", "rank", "ionCount"],
        )
        return pd.concat([df, self._create_generic_peptide_line()])

    def _make_expected_protein_output(self, peptideDf):
        proteinDf = super()._make_expected_protein_output(peptideDf)
        proteinDf["uniquePeptide"] = [0, 0, 1]
        return proteinDf
