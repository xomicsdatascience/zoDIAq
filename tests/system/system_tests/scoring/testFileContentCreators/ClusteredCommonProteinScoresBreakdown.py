import pandas as pd
from .BaselineCommonProteinScoresBreakdown import BaselineCommonProteinScoresBreakdown

"""
NOTE: this is specific to the maxlfq method. This method tests the viability
        of having two unconnected clusters of peptides that correspond to the
        same protein.
 
protein quantification profile:
        peptide1    peptide2    peptide3    peptide4
sample1 100.0       100.0       -           -
sample2 100.0       100.0       -           -
sample3 -           -           200.0       200.0
sample3 -           -           200.0       200.0

protein quantification:
        protein
sample1 100.0
sample2 100.0
sample3 200.0
sample4 200.0
"""


class ClusterSample1And2Breakdown(BaselineCommonProteinScoresBreakdown):
    def _create_input_for_quantifiable_proteins(self):
        df = pd.DataFrame(
            [
                ["peptide1", "1/protein", 1, 0, 0, 100.0],
                ["peptide2", "1/protein", 1, 0, 0, 100.0],
                ["peptide3", "1/protein", 1, 0, 0, 0.0],
                ["peptide4", "1/protein", 1, 0, 0, 0.0],
            ],
            columns=["peptide", "protein", "zLIB", "isDecoy", "rank", "ionCount"],
        )
        return pd.concat([df, self._create_generic_peptide_line()])

    def _make_expected_protein_output(self, peptideDf):
        proteinDf = super()._make_expected_protein_output(peptideDf)
        proteinDf["uniquePeptide"] = [0, 0, 0, 0, 1]
        return proteinDf


class ClusterSample3And4Breakdown(BaselineCommonProteinScoresBreakdown):
    def _create_input_for_quantifiable_proteins(self):
        df = pd.DataFrame(
            [
                ["peptide1", "1/protein", 1, 0, 0, 0.0],
                ["peptide2", "1/protein", 1, 0, 0, 0.0],
                ["peptide3", "1/protein", 1, 0, 0, 200.0],
                ["peptide4", "1/protein", 1, 0, 0, 200.0],
            ],
            columns=["peptide", "protein", "zLIB", "isDecoy", "rank", "ionCount"],
        )
        return pd.concat([df, self._create_generic_peptide_line()])

    def _make_expected_protein_output(self, peptideDf):
        proteinDf = super()._make_expected_protein_output(peptideDf)
        proteinDf["uniquePeptide"] = [0, 0, 0, 0, 1]
        return proteinDf
