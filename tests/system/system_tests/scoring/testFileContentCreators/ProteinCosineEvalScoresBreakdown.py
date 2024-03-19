import pandas as pd
from .MaccScoresBreakdown import MaccScoresBreakdown


class ProteinCosineEvalScoresBreakdown(MaccScoresBreakdown):
    def _create_input_template_for_scoring_module(self):
        """
        For testing the ability of the scoring module to calculate FDRs using macc scoring.
        The created dataframe is broken into the following portions:
            - rank 0 (highest scoring): 1 line. has highest scoring peptide.
            - rank 1: 1 line, total 2. has one peptide with the same protein as rank 0
                (should therefore have same proteinCosine value)
            - rank 2: 200 lines, total 202. has peptides and proteins above the FDR cutoff.
            - rank 3: 2 lines, total 204. decoys that are included in the protein and peptide FDR outputs.
            - rank 4: 1 line, total 205. decoy that is removed in every FDR evaluation.
        """
        sameProteinLayerDf = (
            self._create_two_peptides_from_same_protein_with_different_ranks(
                "000", "001", rank1=0, rank2=1
            )
        )
        peptideInputDf = self._create_input_for_generic_peptides(rank=2, startingRow=2)
        allDecoyDf = self._create_input_for_all_decoys(numRows=3)
        firstDecoyLayerInputDf = allDecoyDf.iloc[:2].copy()
        firstDecoyLayerInputDf["rank"] = [3] * len(firstDecoyLayerInputDf.index)
        secondDecoyLayerInputDf = allDecoyDf.iloc[2:].copy()
        secondDecoyLayerInputDf["rank"] = [4] * len(secondDecoyLayerInputDf.index)
        self.inputDf = pd.concat(
            [
                sameProteinLayerDf,
                peptideInputDf,
                firstDecoyLayerInputDf,
                secondDecoyLayerInputDf,
            ]
        )
        self._add_macc_rank_values_to_input_template_and_remove_rank_label()

    def _create_two_peptides_from_same_protein_with_different_ranks(
        self, peptide1AndProteinId, peptide2Id, rank1, rank2
    ):
        peptides = [f"peptide{peptide1AndProteinId}", f"peptide{peptide2Id}"]
        proteins = [
            f"1/protein{peptide1AndProteinId}",
            f"1/protein{peptide1AndProteinId}",
        ]
        charges = [1, 1]
        isDecoy = [0, 0]
        ranks = [rank1, rank2]
        inputDf = pd.DataFrame()
        inputDf["peptide"] = peptides
        inputDf["protein"] = proteins
        inputDf["zLIB"] = charges
        inputDf["isDecoy"] = isDecoy
        inputDf["rank"] = ranks
        inputDf["ionCount"] = [0] * len(inputDf.index)
        return inputDf

    def _make_expected_spectral_output(self, scoredDf):
        expectedSpectralOutput = scoredDf.copy()
        expectedSpectralOutput = expectedSpectralOutput.iloc[:204, :].reset_index(
            drop=True
        )
        spectralFdr = [0] * len(expectedSpectralOutput.index)
        lastTargetIdx = 202
        decoyFdrs = [i / (lastTargetIdx + i) for i in range(1, 3)]
        spectralFdr[-2:] = decoyFdrs
        expectedSpectralOutput["spectralFDR"] = spectralFdr
        expectedSpectralOutput = expectedSpectralOutput.reset_index(drop=True)
        return expectedSpectralOutput

    def _make_expected_peptide_output(self, scoredDf):
        expectedPeptideOutput = scoredDf.copy()
        expectedPeptideOutput = expectedPeptideOutput.iloc[:204, :].reset_index(
            drop=True
        )
        peptideFdr = [0] * len(expectedPeptideOutput.index)
        lastTargetIdx = 202
        decoyFdrs = [i / (lastTargetIdx + i) for i in range(1, 3)]
        peptideFdr[-2:] = decoyFdrs
        expectedPeptideOutput["peptideFDR"] = peptideFdr
        expectedPeptideOutput = expectedPeptideOutput.reset_index(drop=True)
        return expectedPeptideOutput

    def _make_expected_protein_output(self, peptideDf):
        proteinDf = peptideDf.copy()
        proteinDf["leadingProtein"] = list(proteinDf["protein"])
        proteinCosine = list(proteinDf["cosine"])
        proteinCosine[0] = proteinCosine[1]
        proteinDf["proteinCosine"] = proteinCosine
        numUniqueLeadingProteins = len(set(self.idPickerLeadingProteins))
        lastTargetIdx = len(set(proteinDf["leadingProtein"])) - 2
        proteinFdrs = list(proteinDf["peptideFDR"])
        decoyFdrs = [i / (lastTargetIdx + i) for i in range(1, 3)]
        proteinFdrs[-2:] = decoyFdrs
        proteinDf["leadingProteinFDR"] = proteinFdrs
        uniquePeptides = [0, 0] + [1 for i in range(len(proteinDf.index) - 2)]
        proteinDf["uniquePeptide"] = uniquePeptides
        proteinDf = proteinDf.reset_index(drop=True)
        return proteinDf
