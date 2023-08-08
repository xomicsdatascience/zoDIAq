import pytest
import argparse
from csodiaq import set_command_line_settings
from unittest.mock import Mock


@pytest.fixture
def parser():
    def mock_error(message):
        raise argparse.ArgumentTypeError(message)
    argparse.ArgumentParser.error = Mock(side_effect=mock_error)
    parser = set_command_line_settings()
    return parser

def test__csodiaq__set_command_line_settings__gui__explicitly_flagged_succeeds(parser):
    args = vars(parser.parse_args(["gui"]))
    assert args["command"] == "gui"

def test__csodiaq__set_command_line_settings__gui__succeeds_with_nothing(parser):
    args = vars(parser.parse_args([]))
    assert not args["command"]


@pytest.fixture
def idArgs():
    return ["id"]

@pytest.fixture
def parsedIdArgs(parser, idArgs):
    return vars(parser.parse_args(idArgs))

def test__csodiaq__set_command_line_settings__initialize_identification(parsedIdArgs):
    assert parsedIdArgs["command"] == "id"

@pytest.fixture
def scoreArgs():
    return ["score"]

@pytest.fixture
def parsedScoreArgs(parser, scoreArgs):
    return vars(parser.parse_args(scoreArgs))

def test__csodiaq__set_command_line_settings__initialize_scoring(parsedScoreArgs):
    assert parsedScoreArgs["command"] == "score"

@pytest.fixture
def reanalysisArgs():
    return ["targetedReanalysis"]

@pytest.fixture
def parsedReanalysisArgs(parser, reanalysisArgs):
    return vars(parser.parse_args(reanalysisArgs))

def test__csodiaq__set_command_line_settings__initialize_targeted_reanalysis(parsedReanalysisArgs):
    assert parsedReanalysisArgs["command"] == "targetedReanalysis"
