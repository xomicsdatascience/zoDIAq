from csodiaq import set_command_line_settings
from csodiaq.identification import Identifier
from csodiaq.utils import create_outfile_header

def main():
    parser = set_command_line_settings()
    args = vars(parser.parse_args())
    if args['command'] == 'gui' or args['command'] is None: pass
    elif args['command'] == 'id':
        run_identification(args)
    elif args["command"] == 'score': pass
    elif args["command"] == 'targetedReanalysis': pass

def run_identification(args):
    identifier = Identifier(args)
    for queryFile in args["input"]:
        outputDict = identifier.identify_library_spectra_in_query_file(queryFile)
        outFileHeader = create_outfile_header(args['output'], queryFile, args['correctionDegree'])
        for dfType, df in outputDict.items():
            df.to_csv(f'{outFileHeader}_{dfType}.csv', index=False)

if __name__ == "__main__":
    main()
