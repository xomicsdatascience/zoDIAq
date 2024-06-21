import pandas as pd
from .BaselineScoresBreakdown import BaselineScoresBreakdown


class MaccScoresBreakdown(BaselineScoresBreakdown):
    def _create_input_template_for_scoring_module(self):
        """
        For testing the ability of the scoring module to calculate FDRs using macc scoring.
        The created dataframe is broken into the following portions:
            - rank 0 (highest scoring): 10 lines. has peptides/proteins that are impacted by the id picker algorithm
            - rank 1: 200 lines, total 210. has peptides that are above the FDR cutoff.
            - rank 2: 200 lines, total 410. a duplicate of rank 1, but with different charges. These are cut out to generate the peptide FDR.
            - rank 3: 2 lines, total 412. decoys that are included in the protein and peptide FDR outputs.
            - rank 4: 2 lines, total 414. decoys that are included in the spectral FDR output.
            - rank 5: 1 line, total 415. decoy that is removed in the spectral FDR output.
            - rank 6: 200 lines, total 615. extra peptides that are removed in the spectral FDR output.
        """
        idPickerLayerDf = self._create_id_picker_influenced_peptide_layer(rank=0)
        peptideInputDf = self._create_input_for_generic_peptides(rank=1, startingRow=0)
        duplicatePeptideInputDf = self._create_input_for_generic_peptides(
            rank=2, startingRow=0
        )
        duplicatePeptideInputDf["zLIB"] = [2] * len(duplicatePeptideInputDf.index)
        allDecoyDf = self._create_input_for_all_decoys()
        firstDecoyLayerInputDf = allDecoyDf.iloc[:2].copy()
        firstDecoyLayerInputDf["rank"] = [3] * len(firstDecoyLayerInputDf.index)
        secondDecoyLayerInputDf = allDecoyDf.iloc[2:4].copy()
        secondDecoyLayerInputDf["rank"] = [4] * len(secondDecoyLayerInputDf.index)
        thirdDecoyLayerInputDf = allDecoyDf.iloc[4:].copy()
        thirdDecoyLayerInputDf["rank"] = [5] * len(thirdDecoyLayerInputDf.index)
        extraSpectralInputDf = self._create_input_for_generic_peptides(
            rank=6, startingRow=200
        )
        self.inputDf = pd.concat(
            [
                idPickerLayerDf,
                peptideInputDf,
                duplicatePeptideInputDf,
                firstDecoyLayerInputDf,
                secondDecoyLayerInputDf,
                thirdDecoyLayerInputDf,
                extraSpectralInputDf,
            ]
        )
        self._add_macc_rank_values_to_input_template_and_remove_rank_label()

    def _add_macc_rank_values_to_input_template_and_remove_rank_label(self):
        maccRanking = [
            (14, 0.88),
            (13, 0.89),
            (12, 0.90),
            (11, 0.91),
            (10, 0.92),
            (9, 0.93),
            (8, 0.94),
            (7, 0.95),
            (6, 0.96),
            (5, 0.97),
            (4, 0.98),
            (3, 0.99),
        ]
        maccRankValues = self.inputDf.groupby("rank")["rank"].transform(
            lambda x: [maccRanking[x.iloc[0]]] * len(x)
        )
        self.inputDf["shared"], self.inputDf["cosine"] = zip(*maccRankValues)
        del self.inputDf["rank"]

    def _make_expected_spectral_output(self, scoredDf):
        expectedSpectralOutput = scoredDf.copy()
        expectedSpectralOutput = expectedSpectralOutput.iloc[:414, :].reset_index(
            drop=True
        )
        spectralFdr = [0] * len(expectedSpectralOutput.index)
        lastTargetIdx = 410
        decoyFdrs = [i / (lastTargetIdx + i) for i in range(1, 5)]
        spectralFdr[-4:] = decoyFdrs
        expectedSpectralOutput["spectralFDR"] = spectralFdr
        expectedSpectralOutput = expectedSpectralOutput.reset_index(drop=True)
        return expectedSpectralOutput

    def _make_expected_peptide_output(self, scoredDf):
        expectedPeptideOutput = scoredDf.copy()
        lastTargetIdx = 210
        expectedPeptideOutput = expectedPeptideOutput.iloc[:lastTargetIdx, :]
        decoysInPeptideOutput = scoredDf.iloc[410:412, :]
        expectedPeptideOutput = pd.concat(
            [expectedPeptideOutput, decoysInPeptideOutput]
        ).reset_index(drop=True)
        peptideFdr = [0] * len(expectedPeptideOutput.index)
        decoyFdrs = [i / (lastTargetIdx + i) for i in range(1, 3)]
        peptideFdr[-2:] = decoyFdrs
        expectedPeptideOutput["peptideFDR"] = peptideFdr
        expectedPeptideOutput = expectedPeptideOutput.reset_index(drop=True)
        return expectedPeptideOutput

    def _make_expected_protein_output(self, peptideDf):
        proteinDf = peptideDf.copy()
        proteinDf = proteinDf.reindex(list(proteinDf.index) + [1]).sort_index()
        leadingProteins = list(proteinDf["protein"])
        leadingProteins[: len(self.idPickerLeadingProteins)] = (
            self.idPickerLeadingProteins
        )
        proteinDf["leadingProtein"] = leadingProteins
        proteinDf["proteinCosine"] = proteinDf["cosine"]
        numUniqueLeadingProteins = len(set(self.idPickerLeadingProteins))
        lastTargetIdx = len(set(proteinDf["leadingProtein"])) - 2
        proteinFdrs = list(proteinDf["peptideFDR"])
        decoyFdrs = [i / (lastTargetIdx + i) for i in range(1, 3)]
        proteinFdrs[-2:] = decoyFdrs
        proteinDf["leadingProteinFDR"] = proteinFdrs
        proteinDf["uniquePeptide"] = [0] * len(self.idPickerLeadingProteins) + [1] * (
            len(proteinDf.index) - len(self.idPickerLeadingProteins)
        )
        proteinDf = proteinDf.reset_index(drop=True)
        return proteinDf
