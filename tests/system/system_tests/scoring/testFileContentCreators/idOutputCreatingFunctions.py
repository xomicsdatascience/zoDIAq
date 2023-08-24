import pandas as pd

def create_input_template_for_scoring_module():
    """
    Creates a fabricated version of the output of the identification module ideal for testing the scoring portion.
        Sections are ordered by expected importance by generic rank, so that future scoring mechanisms can attach
        new scoring criteria and still have the output ordered appropriately.
    The created dataframe is broken into the following portions:
        - rank 0 (highest scoring): has peptides/proteins that are impacted by the id picker algorithm
        - rank 1: has peptides that are above the FDR cutoff.
        - rank 2: a duplicate of rank 1, but with different charges. These are cut out to generate the peptide FDR.
        - rank 3: decoys that are included in the protein and peptide FDR outputs.
        - rank 4: decoys that are included in the spectral FDR output.
        - rank 5: extra peptides that are removed in the spectral FDR output.
    """
    idPickerLayerDf = create_id_picker_influenced_peptide_layer(rank=0)
    peptideInputDf = create_input_for_generic_peptides(rank=1, startingRow=0)
    duplicatePeptideInputDf = create_input_for_generic_peptides(rank=2, startingRow=0)
    duplicatePeptideInputDf['zLIB'] = [2]*len(duplicatePeptideInputDf.index)
    allDecoyDf = create_input_for_all_decoys()
    firstDecoyLayerInputDf = allDecoyDf.iloc[:3]
    firstDecoyLayerInputDf['rank'] = [3]*len(firstDecoyLayerInputDf.index)
    secondDecoyLayerInputDf = allDecoyDf.iloc[3:]
    secondDecoyLayerInputDf['rank'] = [4]*len(secondDecoyLayerInputDf.index)
    extraSpectralInputDf = create_input_for_generic_peptides(rank=5, startingRow=200)
    return pd.concat([
        idPickerLayerDf,
        peptideInputDf,
        duplicatePeptideInputDf,
        firstDecoyLayerInputDf,
        secondDecoyLayerInputDf,
        extraSpectralInputDf,
    ])

def create_id_picker_influenced_peptide_layer(rank):
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
    inputDf = pd.DataFrame(inputData, columns = ["peptide","protein"])
    inputDf["zLIB"] = [1] * len(inputDf.index)
    inputDf["isDecoy"] = [0] * len(inputDf.index)
    inputDf["rank"] = [rank] * len(inputDf.index)
    return inputDf

def create_input_for_generic_peptides(rank, startingRow):
    numRows = 200
    numbers = [str(i).rjust(3,'0') for i in range(startingRow, 200+startingRow)]
    peptides = [f'peptide{number}' for number in numbers]
    proteins = [f'protein{number}' for number in numbers]
    charges = [1]*numRows
    isDecoy = [0]*numRows
    ranks = [rank]*numRows
    inputDf = pd.DataFrame()
    inputDf['peptide'] = peptides
    inputDf['protein'] = proteins
    inputDf['zLIB'] = charges
    inputDf['isDecoy'] = isDecoy
    inputDf['rank'] = ranks
    return inputDf

def create_input_for_all_decoys():
    numRows = 5
    numbers = [str(i).rjust(2,'0') for i in range(5)]
    peptides = [f'DECOY_peptide{number}' for number in numbers]
    proteins = [f'DECOY_protein{number}' for number in numbers]
    inputDf = pd.DataFrame()
    inputDf['peptide'] = peptides
    inputDf['protein'] = proteins
    inputDf['zLIB'] = [1]*len(inputDf.index)
    inputDf['isDecoy'] = [1]*len(inputDf.index)
    return inputDf

def create_expected_protein_output():
    outputData = [
        ["1/protein7", 0],
        ["2/protein4/protein9", 1], #peptide02
        ["1/protein6", 0], #peptide02 duplicate row
        ["1/protein1", 0],
        ["1/protein1", 0],
        ["1/protein7", 0],
        ["1/protein6", 0],
        ["1/protein7", 0],
        ["1/protein7", 0],
        ["1/protein7", 0],
    ]

