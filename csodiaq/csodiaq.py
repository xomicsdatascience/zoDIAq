import os
import pandas as pd
from csodiaq import set_command_line_settings
from csodiaq.identification import Identifier
from csodiaq.utils import create_outfile_header
from csodiaq.scoring import create_spectral_fdr_output_from_full_output, create_peptide_fdr_output_from_full_output, create_protein_fdr_output_from_peptide_fdr_output

def main():
    parser = set_command_line_settings()
    args = vars(parser.parse_args())
    if args['command'] == 'gui' or args['command'] is None: pass
    elif args['command'] == 'id':
        run_identification(args)
    elif args["command"] == 'score':
        run_scoring(args)
    elif args["command"] == 'targetedReanalysis': pass

def run_identification(args):
    identifier = Identifier(args)
    for queryFile in args["input"]:
        identificationFullOutputDf = identifier.identify_library_spectra_in_query_file(queryFile)
        outFileHeader = create_outfile_header(args['output'], queryFile, args['correctionDegree'])
        identificationFullOutputDf.to_csv(f'{outFileHeader}_fullOutput.csv', index=False)

def run_scoring(args):
    outputDir = os.path.join(args["input"]["identificationDirectory"], f'fdrScores-{args["score"]}')
    if not os.path.exists(outputDir):
        os.mkdir(outputDir)
    for idDfFile in args["input"]["idFiles"]:
        fileHeader = extract_file_name_without_file_type(idDfFile)
        idDf = pd.read_csv(os.path.join(args["input"]["identificationDirectory"],idDfFile))
        spectralDf = create_spectral_fdr_output_from_full_output(idDf)
        spectralDf.to_csv(os.path.join(outputDir, f'{fileHeader}_spectralFDR.csv'), index=False)
        peptideDf = create_peptide_fdr_output_from_full_output(idDf)
        peptideDf.to_csv(os.path.join(outputDir, f'{fileHeader}_peptideFDR.csv'), index=False)
        proteinDf = create_protein_fdr_output_from_peptide_fdr_output(peptideDf) #TODO: This should be conditional on the presence of an appropriately-formatted protein column
        proteinDf.to_csv(os.path.join(outputDir, f'{fileHeader}_proteinFDR.csv'), index=False)

def extract_file_name_without_file_type(file):
    return '.'.join(file.split('.')[:-1])

if __name__ == "__main__":
    main()
