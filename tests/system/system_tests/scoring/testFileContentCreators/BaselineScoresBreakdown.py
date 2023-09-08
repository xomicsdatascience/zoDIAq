from abc import ABC, abstractmethod
import pandas as pd
import os

class BaselineScoresBreakdown(ABC):
    def __init__(self, expectedOutputDirectory):
        self.idPickerLeadingProteins = [
            "1/protein7",
            "2/protein4/protein9",
            "1/protein6",
            "1/protein1",
            "1/protein1",
            "1/protein7",
            "1/protein6",
            "1/protein1",
            "1/protein1",
            "1/protein1",
            "2/protein4/protein9",
        ]
        self._create_input_template_for_scoring_module()
        self._create_outputs_of_scoring_module_and_write_to_file(expectedOutputDirectory)

    @abstractmethod
    def _create_input_template_for_scoring_module(self):
        """
        Creates a fabricated version of the output of the identification module ideal for testing the scoring portion.
            Sections are ordered by expected importance by generic rank, so that future scoring mechanisms can attach
            new scoring criteria and still have the output ordered appropriately.
        """
        pass

    def _create_outputs_of_scoring_module_and_write_to_file(self, expectedOutputDirectory):
        expectedSpectralOutput = self._make_expected_spectral_output(self.inputDf)
        expectedPeptideOutput = self._make_expected_peptide_output(self.inputDf)
        expectedProteinOutput = self._make_expected_protein_output(expectedPeptideOutput)
        self.outputDict = {
            'spectralFDR': expectedSpectralOutput,
            'peptideFDR': expectedPeptideOutput,
            'proteinFDR': expectedProteinOutput,
        }
        for key,value in self.outputDict.items():
            if isinstance(value, pd.DataFrame):
                value.to_csv(os.path.join(expectedOutputDirectory, f'{key}.csv'), index=False)

    @abstractmethod
    def _make_expected_spectral_output(self, scoredDf):
        pass

    @abstractmethod
    def _make_expected_peptide_output(self, scoredDf):
        pass

    @abstractmethod
    def _make_expected_protein_output(self, peptideDf):
        pass

    def _create_id_picker_influenced_peptide_layer(self, rank):
        inputData = [
            ["peptide01", "1/protein7"],
            ["peptide02", "3/protein4/protein6/protein9"],
            ["peptide03", "1/protein1"],
            ["peptide04", "2/protein1/protein5"],
            ["peptide05", "1/protein7"],
            ["peptide06", "2/protein3/protein6"],
            ["peptide07", "1/protein1"],
            ["peptide08", "4/protein1/protein2/protein5/protein8"],
            ["peptide09", "1/protein1"],
            ["peptide10", "2/protein4/protein9"],
        ]
        inputDf = pd.DataFrame(inputData, columns=["peptide", "protein"])
        inputDf["zLIB"] = [1] * len(inputDf.index)
        inputDf["isDecoy"] = [0] * len(inputDf.index)
        inputDf["rank"] = [rank] * len(inputDf.index)
        inputDf['ionCount'] = [100.0] * len(inputDf.index)

        return inputDf

    def _create_input_for_generic_peptides(self, rank, startingRow, numRows=200):
        numbers = [str(i).rjust(3, '0') for i in range(startingRow, numRows + startingRow)]
        peptides = [f'peptide{number}' for number in numbers]
        proteins = [f'1/protein{number}' for number in numbers]
        charges = [1] * numRows
        isDecoy = [0] * numRows
        ranks = [rank] * numRows
        inputDf = pd.DataFrame()
        inputDf['peptide'] = peptides
        inputDf['protein'] = proteins
        inputDf['zLIB'] = charges
        inputDf['isDecoy'] = isDecoy
        inputDf['rank'] = ranks
        inputDf['ionCount'] = [100.0] * len(inputDf.index)
        return inputDf

    def _create_input_for_all_decoys(self, numRows = 5):
        numbers = [str(i).rjust(len(str(numRows)),'0') for i in range(numRows)]
        peptides = [f'DECOY_peptide{number}' for number in numbers]
        proteins = [f'1/DECOY_protein{number}' for number in numbers]
        inputDf = pd.DataFrame()
        inputDf['peptide'] = peptides
        inputDf['protein'] = proteins
        inputDf['zLIB'] = [1]*len(inputDf.index)
        inputDf['isDecoy'] = [1]*len(inputDf.index)
        inputDf['ionCount'] = [100.0]*len(inputDf.index)
        return inputDf
