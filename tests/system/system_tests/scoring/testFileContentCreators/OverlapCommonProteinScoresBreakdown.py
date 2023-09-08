import pandas as pd
from .BaselineCommonProteinScoresBreakdown import BaselineCommonProteinScoresBreakdown

"""
NOTE: this is specific to the maxlfq method. Samples with 1 or 0 peptides 
    shared with other samples are considered 0.
 
protein quantification profile for protein 1:
        peptide1    peptide2    peptide3
sample1 X           X           X
sample2 X           X           X
sample3 -           X           X

protein quantification profile for protein 2:
        peptide2    peptide3    peptide4
sample1 X           X           -
sample2 X           X           -
sample3 X           X           X

protein 1 quantification:
        protein
sample1 X
sample2 X
sample3 0

protein 2 quantification for average method:
        protein
sample1 0
sample2 0
sample3 X

protein 2 quantification for maxlfq method:
        protein
sample1 0
sample2 0
sample3 0
"""


class OverlapSample1And2Breakdown(BaselineCommonProteinScoresBreakdown):
    def _create_input_for_quantifiable_proteins(self):
        df = pd.DataFrame(
            [
                ["peptide1", "1/protein1", 1, 0, 0, 100.0],
                ["peptide2", "2/protein1/protein2", 1, 0, 0, 100.0],
                ["peptide3", "2/protein1/protein2", 1, 0, 0, 100.0],
            ],
            columns=["peptide", "protein", "zLIB", "isDecoy", "rank", "ionCount"],
        )
        return pd.concat([df, self._create_generic_peptide_line()])

    def _make_expected_protein_output(self, peptideDf):
        proteinDf = super()._make_expected_protein_output(peptideDf)
        proteinDf["uniquePeptide"] = [1, 0, 0, 1]
        return proteinDf


class OvelapSample3Breakdown(BaselineCommonProteinScoresBreakdown):
    def _create_input_for_quantifiable_proteins(self):
        df = pd.DataFrame(
            [
                ["peptide2", "2/protein1/protein2", 1, 0, 0, 100.0],
                ["peptide3", "2/protein1/protein2", 1, 0, 0, 100.0],
                ["peptide4", "1/protein2", 1, 0, 0, 100.0],
            ],
            columns=["peptide", "protein", "zLIB", "isDecoy", "rank", "ionCount"],
        )
        return pd.concat([df, self._create_generic_peptide_line()])

    def _make_expected_protein_output(self, peptideDf):
        proteinDf = super()._make_expected_protein_output(peptideDf)
        proteinDf["uniquePeptide"] = [0, 0, 1, 1]
        return proteinDf
