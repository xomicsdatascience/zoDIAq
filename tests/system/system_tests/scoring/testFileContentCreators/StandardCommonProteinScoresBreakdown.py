import pandas as pd
from .BaselineCommonProteinScoresBreakdown import BaselineCommonProteinScoresBreakdown

"""
protein quantification profile:
        peptide1    peptide2 
sample1 X           X
sample2 X           X
sample3 -           -

protein quantification:
        protein
sample1 X
sample2 X
sample3 0
"""


class StandardSample1And2Breakdown(BaselineCommonProteinScoresBreakdown):
    def _create_input_for_quantifiable_proteins(self):
        dfs = [
            pd.DataFrame(
                [[f"peptide{i}", f"1/protein", 1, 0, 0, 100.0]],
                columns=["peptide", "protein", "zLIB", "isDecoy", "rank", "ionCount"],
            )
            for i in range(2)
        ]
        dfs.append(self._create_generic_peptide_line())
        return pd.concat(dfs)

    def _make_expected_protein_output(self, peptideDf):
        proteinDf = super()._make_expected_protein_output(peptideDf)
        proteinDf["uniquePeptide"] = [0, 0, 1]


class StandardSample3Breakdown(BaselineCommonProteinScoresBreakdown):
    def _create_input_for_quantifiable_proteins(self):
        return self._create_generic_peptide_line()

    def _make_expected_protein_output(self, peptideDf):
        proteinDf = super()._make_expected_protein_output(peptideDf)
        proteinDf["uniquePeptide"] = [1]
