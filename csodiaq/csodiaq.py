#!/usr/bin/env python
import sys
import csv
import argparse
import os
from os import path
from os.path import join
from . import csodiaq_gui as gui
import pickle
import numpy as np
from . import csodiaq_identification_functions as cif
from . import csodiaq_quantification_functions as cqf
from . import csodiaq_mgf_cleaning_functions as cmf
from . import peptide_quantification
from csodiaq.identifier import Identifier, identify_high_confidence_proteins, calculate_fdr_rates_of_decoy_array
from csodiaq.utils import create_outfile_header
import pandas as pd

def main():
    arg_parser = set_command_line_settings()
    args = vars(arg_parser.parse_args())
    check_for_invalid_input(arg_parser, args)

    if args['command'] == 'gui' or args['command'] is None:  # Default to GUI
        gui.main()
        return

    if args['command'] == 'mgf':
        cmf.clean_mgf_file(args['mgf'], args['fasta'], ions=args['ions'])
        return

    if args['command'] == 'id':
        identifier = Identifier(args)
        for queryFile in args["files"]:
            idDf = identifier.identify_library_spectra_in_query_file(queryFile)
            outFileHeader = create_outfile_header(args['outDirectory'], queryFile, args['correction'])
            write_identification_outputs(idDf, outFileHeader, args['proteinTargets'])
        #lib = cif.library_file_to_dict(args['library'])
        #maxQuerySpectraToPool = queryPooling = args['query']
        #if not maxQuerySpectraToPool:
        #    maxQuerySpectraToPool = np.inf
        for i in range(len(args['files'])):
            outFileHeader = args['outDirectory'] + 'CsoDIAq-file' + '_' + '.'.join(args['files'][i].split('/')[-1].split('.')[:-1])
            if args['correction'] != -1:
                outFileHeader += '_corrected'
            outFile = outFileHeader + '.csv'
            '''
            if args['histogram']:
                histFile = outFileHeader + '_histogram.png'
            else:
                histFile = ''
            if args['peaks'] == 0:
                args['peaks'] = np.inf

            cif.perform_spectra_pooling_and_analysis(args['files'][i],
                                                     outFile,
                                                     lib,
                                                     args['fragmentMassTolerance'],
                                                     maxQuerySpectraToPool,
                                                     args['correction'],
                                                     histFile,
                                                     num_peaks=args['peaks'])
            '''
            spectralFile = outFileHeader + '_spectralFDR.csv'
            peptideFile = outFileHeader + '_peptideFDR.csv'
            if args['proteinTargets']:
                proteinFile = outFileHeader + '_proteinFDR.csv'
            else:
                proteinFile = ''
            #cif.write_fdr_outputs(outFile, spectralFile,
            #                      peptideFile, proteinFile)

            reanalysisHeader = outFileHeader + '_mostIntenseTargs'
            if args['proteinTargets']:
                inFile = proteinFile
            else:
                inFile = peptideFile
            fdrDf = cif.filter_fdr_output_for_targeted_reanalysis(
                inFile, args['proteinTargets'], args['heavyMz'])
            cif.write_targeted_reanalysis_outputs(
                reanalysisHeader, fdrDf, args['heavyMz'])

        # Post-processing
        # Check whether we should do anything
        if args['commonpeptide']:
            print(f'Performing protein quantification and identifying common peptpides...', flush=True)
            peptide_quantification.get_all_peptide_quantities(file_list=args['files'],
                                                              library_file=args['library'],
                                                              csodiaq_output_dir=args['outDirectory'],
                                                              num_library_fragments=args['peaks'],
                                                              save_file=os.path.join(args['outDirectory'], 'common_peptides.csv'))
            print('Done.')
        if args['commonprotein']:
            print('Extracting proteins common across output files...', flush=True)
            peptide_quantification.get_protein_quantities(file_list=args['files'],
                                                          library_file=args['library'],
                                                          csodiaq_output_dir=args['outDirectory'],
                                                          num_library_fragments=args['peaks'],
                                                          save_file=os.path.join(args['outDirectory'], 'common_proteins.csv'))
            print('Done')

    if args['command'] == 'quant':
        scanToCsodiaqDict, scanToLibPeaksDict = cqf.connect_mzxml_to_csodiaq_and_library(
            args['idFile'], args['library'], args['files'], args['libraryPeaks'])

        if args['histogram']:
            hist = args['outDirectory'] + 'SILAC_Quantification_histogram.png'
        else:
            hist = ''
        df = cqf.heavy_light_quantification(scanToCsodiaqDict, scanToLibPeaksDict, args['files'], args['outDirectory'], args[
                                            'fragmentMassTolerance'], args['minimumMatches'], args['ratioType'], args['correction'], hist)

        df.to_csv(args['outDirectory'] +
                  'CsoDIAq_output_SILAC_Quantification.csv', index=False)


def set_command_line_settings():
    arg_parser = argparse.ArgumentParser(description='')
    subparsers = arg_parser.add_subparsers(
        dest='command', help='CsoDIAq Functions')
    arg_parser.version = '0.1'
    GUI_parser = subparsers.add_parser(
        'gui', help='Launches the (optional) GUI application for using CsoDIAq.')
    mgf_parser = subparsers.add_parser(
        'mgf', help='Processes an MGF file by adding proteins (from a given FASTA file) and optionally filtering all peaks that don\'t correspond to a y or b ion.')
    mgf_parser.add_argument('-m', '--mgf', type=str,
                            required=True, help='MGF file to process.')
    mgf_parser.add_argument('-f', '--fasta', type=str, required=True,
                            help='Protein FASTA reference file for matching peptides to proteins.')
    mgf_parser.add_argument('-i', '--ions', action='store_true',
                            help='Filter by Y and B ions - This flag indicates if only Y or B ion peaks corresponding to the given peptide should be included in the library.')
    parent_parser = argparse.ArgumentParser(add_help=False)
    parent_parser.add_argument('-f', '--files', type=str, action='append', required=True,
                               help='DISPA Data Files - argument can be used multiple times for multiple files.')
    parent_parser.add_argument('-l', '--library', type=str, required=True,
                               help='Library Spectra Data File - CsoDIAq accepts MGF and TraML formats.')
    parent_parser.add_argument('-o', '--outDirectory', type=str, required=True,
                               help='Outfile Directory - Folder to write output files to.')
    parent_parser.add_argument('-t', '--fragmentMassTolerance', type=restricted_int, default=30,
                               help='Fragment Mass Tolerance - Initial tolerance allowed between mz values of peaks to be considered a match (in PPM). Default value is 20.')
    parent_parser.add_argument('-c', '--correction', type=restricted_float, nargs='?', const=0, default=-1, help='PPM Correction - The presence of this flag indicates that a second, corrected analysis should be performed. During the initial uncorrected run, PPM differences between all matched peaks for peptides of interest are stored. Generally, these PPM differences concentrate at a non-zero location. Thus, a second, corrected analysis using this point of concentration can provide more accurate results. Floats between 0.5 and 2 can be accepted as arguments, indicating the number of standard deviations around the mean that will be used as the new tolerance (see the fragmentMassTolerance flag). Alternatively, if the flag is present but no standard deviation is specified, a customized offset and tolerance based on the histogram of the uncorrected analysis is used instead. In this analysis, offset is the PPM value of the highest peak, and tolerance is the width of the peak before bars are below the average of the bars representing noise.')
    parent_parser.add_argument('-hist', '--histogram', action='store_true',
                               help='Histogram - This flag indicates a histogram of the uncorrected PPM values (with lines for the chosen offset/tolerance) should be generated. Flag can only be used if the -c or --correction flag is used.')

    id_parser = subparsers.add_parser('id', parents=[
                                      parent_parser], help='Identify peptides and/or proteins from DISPA data and prepare files for DISPA re-analysis.')
    quant_parser = subparsers.add_parser('quant', parents=[
                                         parent_parser], help='Quantify peptides and/or proteins from DISPA targetted re-analysis data corresponding to a SILAC experiment.')

    id_parser.add_argument('-p', '--proteinTargets', type=restricted_int, default=0,
                           help='Protein Targets - The number of peptide targets per protein that should be written to the DISPA re-analysis files. Value must be a positive, non-zero integer. When not specified, the program will do an untargeted peptide analysis by default.')
    id_parser.add_argument('-q', '--query', type=restricted_int, default=0,
                           help='Enable Query Pooling - This flag indicates that query spectra should NOT be pooled. This action increases the overall time complexity of the algorithm while reducing memory complexity.')
    id_parser.add_argument('-heavy', '--heavyMz', action='store_true',
                           help='Heavy Mz - Files for targetted re-analysis include heavy fragment isotopes for SILAC quantification.')
    id_parser.add_argument('--peaks', type=int, default=0,
                           help='Number of peaks from the query spectrum to use for quantification of ions.')
    id_parser.add_argument('--commonpeptide', default=False, action='store_true')
    id_parser.add_argument('--commonprotein', default=False, action='store_true')

    quant_parser.add_argument('-i', '--idFile', type=str, required=True,
                              help='Protein/Identification File - Output from the identification portion of QsoDIAq. The file ending will have "all_CVs.csv"')
    quant_parser.add_argument('-p', '--libraryPeaks', type=restricted_int, default=0,
                              help='Maximum Number of Library Peaks - Maximum number of library peaks allowed in library spectra, prioritizing intensity. Program allows all peaks by default.')
    quant_parser.add_argument('-m', '--minimumMatches', type=restricted_int, default=0,
                              help='Minimum Number of Matching Peaks - Minimum number of matched peaks . Program requires at least one match from the top 3 most intense peaks in the library spectrum by default.')
    quant_parser.add_argument('-r', '--ratioType', type=restricted_ratio_type, default='median',
                              help='Ratio Determination Method - method for picking the ratio for a given peptide when several peaks match. Options are median or mean, default is median.')
    return arg_parser


def check_for_invalid_input(arg_parser, args):
    if args['command'] == 'id' or args['command'] == 'quant':
        if (args['histogram'] and args['correction'] == -1):
            arg_parser.error(
                'The -hist or --histogram argument requires the -c or -correction argument')
        for file in args['files']:
            if not restricted_file(file, permittedTypes=['mzxml']):
                arg_parser.error(
                    'The -f or --files argument must be an existing file of type mzxml')
        if not restricted_file(args['library'], permittedTypes=['mgf', 'csv', 'tsv']):
            arg_parser.error(
                'The -l or --library argument must be an existing file of type MGF (.mgf) or TraML (.tsv or .csv) format')
        if not restricted_file(args['outDirectory']):
            arg_parser.error(
                'The -o or --outDirectory argument must be an existing directory')


def restricted_float(x):
    try:
        x = float(x)
    except ValueError:
        raise argparse.ArgumentTypeError(
            "%r not a floating-point literal" % (x,))

    if x < 0.5 or x > 2.0:
        raise argparse.ArgumentTypeError("%r not between 0.5 and 2.0" % (x,))
    return x


def restricted_int(x):
    try:
        x = int(x)
    except ValueError:
        raise argparse.ArgumentTypeError("%r not an integer literal" % (x,))

    if x < 1:
        raise argparse.ArgumentTypeError(
            "%r must be an integer greater than 0" % (x,))
    return x


def restricted_file(x, permittedTypes=[]):
    if not len(permittedTypes):
        if not path.isdir(x):
            return False
        else:
            return True
    if x.split('.')[-1].lower() not in permittedTypes or not path.isfile(x):
        return False
    return True


def restricted_ratio_type(x):
    allowed_types = ['median', 'mean']
    if x not in allowed_types:
        raise argparse.ArgumentTypeError(
            "%r not a valid ratio return type (allowed types: mean, median)" % (x,))
    return x


if __name__ == "__main__":
    main()
