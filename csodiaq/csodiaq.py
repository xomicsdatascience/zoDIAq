import argparse


def set_command_line_settings():
    parser = argparse.ArgumentParser(description="")
    commandParser = parser.add_subparsers(dest="command", help="CsoDIAq Functions")
    guiParser = commandParser.add_parser(
        'gui', help='Launches the (optional) GUI application for using CsoDIAq.')
    idParser = commandParser.add_parser('id', help='Identify peptides from a designated peptide library in mass spectrometry run (query) files.')
    scoringParser = commandParser.add_parser('score', help='Scores confidence in identified peptides and removes those above a False Discovery Rate (FDR) of 0.01.\nNOTE: The identification output must include decoys for accurate FDR scoring.')
    reanalysisParser = commandParser.add_parser('targetedReanalysis', help='Creates files readable by a mass spectrometer for targeted reanalysis of identified peptides.')
    return parser
