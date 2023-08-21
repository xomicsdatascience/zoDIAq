import os
import warnings

import pandas as pd
from csodiaq import set_args_from_command_line_input, check_for_conflicting_args
from csodiaq.identification import Identifier
from csodiaq.utils import (
    create_outfile_header,
    confirm_proteins_in_list_are_in_appropriate_format,
    Printer,
)
from csodiaq.scoring import (
    create_spectral_fdr_output_from_full_output,
    create_peptide_fdr_output_from_full_output,
    create_protein_fdr_output_from_peptide_fdr_output,
    calculate_ion_count_for_each_protein_in_protein_fdr_df,
    compile_ion_count_comparison_across_runs_df,
)
from csodiaq.targetedReanalysis import (
    create_mass_spec_input_dataframes_for_targeted_reanalysis_of_identified_peptides,
)
from csodiaq.gui import run_gui


def main():
    parser = set_args_from_command_line_input()
    args = vars(parser.parse_args())
    check_for_conflicting_args(args)
    if args["command"] == "gui" or args["command"] is None:
        run_gui()
    elif args["command"] == "id":
        run_identification(args)
    elif args["command"] == "score":
        run_scoring(args)
    elif args["command"] == "targetedReanalysis":
        run_targeted_reanalysis(args)


def run_identification(args):
    printer = Printer()
    printer("Begin Peptide Identification Process")
    identifier = Identifier(args)
    os.mkdir(args["output"])
    for queryFile in args["input"]:
        printer(f"Beginning Identification for {queryFile} input file")
        identificationFullOutputDf = identifier.identify_library_spectra_in_query_file(
            queryFile
        )
        if isinstance(identificationFullOutputDf, str):
            warnings.warn(f'{identificationFullOutputDf} Skipping {queryFile} file.', UserWarning)
            continue

        outFileHeader = create_outfile_header(
            args["output"], queryFile, args["correctionDegree"]
        )
        identificationFullOutputDf.to_csv(
            f"{outFileHeader}_fullOutput.csv", index=False
        )
    printer("End Peptide Identification Process")


def run_scoring(args):
    printer = Printer()
    printer("Begin Scoring")
    outputDir = os.path.join(
        args["input"]["csodiaqDirectory"], f'fdrScores-{args["score"]}'
    )
    if not os.path.exists(outputDir):
        os.mkdir(outputDir)
    peptideDfs = {}
    proteinDfs = {}
    for idDfFile in args["input"]["idFiles"]:
        fileHeader = extract_file_name_without_file_type(idDfFile)
        idDf = pd.read_csv(os.path.join(args["input"]["csodiaqDirectory"], idDfFile))
        spectralDf = create_spectral_fdr_output_from_full_output(idDf)
        spectralDf.to_csv(
            os.path.join(outputDir, f"{fileHeader}_spectralFDR.csv"), index=False
        )
        peptideDf = create_peptide_fdr_output_from_full_output(idDf)
        peptideDf.to_csv(
            os.path.join(outputDir, f"{fileHeader}_peptideFDR.csv"), index=False
        )
        peptideDfs[fileHeader] = peptideDf[["peptide", "ionCount"]]
        if confirm_proteins_in_list_are_in_appropriate_format(peptideDf["protein"]):
            proteinDf = create_protein_fdr_output_from_peptide_fdr_output(peptideDf)
            proteinDf.to_csv(
                os.path.join(outputDir, f"{fileHeader}_proteinFDR.csv"), index=False
            )
            proteinDfs[
                fileHeader
            ] = calculate_ion_count_for_each_protein_in_protein_fdr_df(proteinDf)
    printer("Begin Quantifying Common Peptides")
    commonPeptideDf = compile_ion_count_comparison_across_runs_df(peptideDfs, "peptide")
    commonPeptideDf.to_csv(os.path.join(outputDir, "commonPeptides.csv"))
    if len(proteinDfs) > 0:
        printer("Begin Quantifying Common Proteins")
        commonProteinDf = compile_ion_count_comparison_across_runs_df(
            proteinDfs, "protein"
        )
        commonProteinDf.to_csv(os.path.join(outputDir, "commonProteins.csv"))
    printer("Finish Scoring")


def run_targeted_reanalysis(args):
    printer = Printer()
    printer("Begin Targeted Reanalysis File Generation")
    outputDir = make_targeted_reanalysis_output_directory_name(args)
    if not os.path.exists(outputDir):
        os.mkdir(outputDir)
    if args["protein"]:
        scoreType = "protein"
    else:
        scoreType = "peptide"
    for scoreFdrFile in args["input"][scoreType]:
        fileHeader = extract_file_name_without_file_type(scoreFdrFile)
        scoreDf = pd.read_csv(
            os.path.join(args["input"]["csodiaqDirectory"], scoreFdrFile)
        )
        targetedOutputDict = create_mass_spec_input_dataframes_for_targeted_reanalysis_of_identified_peptides(
            scoreDf,
            isIncludeHeavyIsotopes=args["heavyIsotope"],
            maximumPeptidesPerProtein=args["protein"],
            binValueProximity=args["binValueProximity"],
        )
        for name, df in targetedOutputDict.items():
            if name == "fullDf":
                df.to_csv(
                    os.path.join(outputDir, f"{name}_{fileHeader}.csv"), index=False
                )
            else:
                df.to_csv(
                    os.path.join(outputDir, f"{name}_{fileHeader}.txt"),
                    sep="\t",
                    index=False,
                )
    printer("End Targeted Reanalysis File Generation")


def make_targeted_reanalysis_output_directory_name(args):
    if args["protein"]:
        proteinHeader = f"maxPeptidesPerProtein{args['protein']}"
    else:
        proteinHeader = "peptidesNoProteins"
    if args["heavyIsotope"]:
        heavyHeader = "includesHeavyIsotopes"
    else:
        heavyHeader = "noHeavyIsotopes"
    binWidthHeader = (
        f'binValueProximity{str(args["binValueProximity"]).replace(".","p")}'
    )
    newDirectoryName = (
        f"targetedReanalysis_{proteinHeader}_{heavyHeader}_{binWidthHeader}"
    )
    return os.path.join(args["input"]["csodiaqDirectory"], newDirectoryName)


def extract_file_name_without_file_type(file):
    return ".".join(file.split(".")[:-1])


if __name__ == "__main__":
    main()
