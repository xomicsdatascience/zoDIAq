import pandas as pd
from abc import abstractmethod
from .BaselineCommonProteinScoresBreakdown import BaselineCommonProteinScoresBreakdown

"""
protein quantification profile:
        peptide1    peptide2    peptide3
sample1 100.0       200.0       300.0
sample2 100.0       200.0       400.0
sample3 100.0       200.0       500.0

protein quantification for 'average' method:
        protein
sample1 600.0
sample2 700.0
sample3 800.0

protein quantification for 'maxlfq' method:
        protein
sample1 391.479818
sample2 391.487326
sample3 391.493149
"""


class MethodSampleBreakdown(BaselineCommonProteinScoresBreakdown):
    def _create_input_for_quantifiable_proteins(self):
        df = pd.DataFrame(
            [
                [f"peptide{i}", f"1/protein", 1, 0, 0, self._get_ion_counts()[i]]
                for i in range(3)
            ],
            columns=["peptide", "protein", "zLIB", "isDecoy", "rank", "ionCount"],
        )
        return pd.concat([df, self._create_generic_peptide_line()])

    @abstractmethod
    def _get_ion_counts(self):
        pass

    def _make_expected_protein_output(self, peptideDf):
        proteinDf = super()._make_expected_protein_output(peptideDf)
        proteinDf["uniquePeptide"] = [0, 0, 0, 0]


class MethodSample1Breakdown(MethodSampleBreakdown):
    def _get_ion_counts(self):
        return [100.0, 200.0, 300.0]


class MethodSample2Breakdown(MethodSampleBreakdown):
    def _get_ion_counts(self):
        return [100.0, 200.0, 400.0]


class MethodSample3Breakdown(MethodSampleBreakdown):
    def _get_ion_counts(self):
        return [100.0, 200.0, 500.0]
