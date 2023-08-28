import pandas as pd

def create_input_template_for_scoring_module():
    """
    Creates a fabricated version of the output of the identification module ideal for testing the scoring portion.
        Sections are ordered by expected importance by generic rank, so that future scoring mechanisms can attach
        new scoring criteria and still have the output ordered appropriately.
    The created dataframe is broken into the following portions:
        - rank 0 (highest scoring): 10 lines. has peptides/proteins that are impacted by the id picker algorithm
        - rank 1: 200 lines, total 210. has peptides that are above the FDR cutoff.
        - rank 2: 200 lines, total 410. a duplicate of rank 1, but with different charges. These are cut out to generate the peptide FDR.
        - rank 3: 2 lines, total 412. decoys that are included in the protein and peptide FDR outputs.
        - rank 4: 2 lines, total 414. decoys that are included in the spectral FDR output.
        - rank 5: 1 line, total 415. decoy that is removed in the spectral FDR output.
        - rank 6: 200 lines, total 615. extra peptides that are removed in the spectral FDR output.
    """
    idPickerLayerDf = create_id_picker_influenced_peptide_layer(rank=0)
    peptideInputDf = create_input_for_generic_peptides(rank=1, startingRow=0)
    duplicatePeptideInputDf = create_input_for_generic_peptides(rank=2, startingRow=0)
    duplicatePeptideInputDf['zLIB'] = [2]*len(duplicatePeptideInputDf.index)
    allDecoyDf = create_input_for_all_decoys()
    firstDecoyLayerInputDf = allDecoyDf.iloc[:2]
    firstDecoyLayerInputDf['rank'] = [3] * len(firstDecoyLayerInputDf.index)
    secondDecoyLayerInputDf = allDecoyDf.iloc[2:4]
    secondDecoyLayerInputDf['rank'] = [4] * len(secondDecoyLayerInputDf.index)
    thirdDecoyLayerInputDf = allDecoyDf.iloc[4:]
    thirdDecoyLayerInputDf['rank'] = [5] * len(thirdDecoyLayerInputDf.index)
    extraSpectralInputDf = create_input_for_generic_peptides(rank=6, startingRow=200)
    return pd.concat([
        idPickerLayerDf,
        peptideInputDf,
        duplicatePeptideInputDf,
        firstDecoyLayerInputDf,
        secondDecoyLayerInputDf,
        thirdDecoyLayerInputDf,
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
    inputDf['ionCount'] = [0]*len(inputDf.index)

    return inputDf

def create_input_for_generic_peptides(rank, startingRow):
    numRows = 200
    numbers = [str(i).rjust(3,'0') for i in range(startingRow, 200+startingRow)]
    peptides = [f'peptide{number}' for number in numbers]
    proteins = [f'1/protein{number}' for number in numbers]
    charges = [1]*numRows
    isDecoy = [0]*numRows
    ranks = [rank]*numRows
    inputDf = pd.DataFrame()
    inputDf['peptide'] = peptides
    inputDf['protein'] = proteins
    inputDf['zLIB'] = charges
    inputDf['isDecoy'] = isDecoy
    inputDf['rank'] = ranks
    inputDf['ionCount'] = [0]*len(inputDf.index)

    return inputDf

def create_input_for_all_decoys():
    numRows = 5
    numbers = [str(i).rjust(2,'0') for i in range(5)]
    peptides = [f'DECOY_peptide{number}' for number in numbers]
    proteins = [f'1/DECOY_protein{number}' for number in numbers]
    inputDf = pd.DataFrame()
    inputDf['peptide'] = peptides
    inputDf['protein'] = proteins
    inputDf['zLIB'] = [1]*len(inputDf.index)
    inputDf['isDecoy'] = [1]*len(inputDf.index)
    inputDf['ionCount'] = [0]*len(inputDf.index)
    return inputDf
