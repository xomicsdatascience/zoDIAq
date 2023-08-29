import os
from tempfile import TemporaryDirectory
import pytest

@pytest.fixture(scope="module")
def systemTestFileDirectory():
    return TemporaryDirectory(prefix="csodiaq_scoring_system_test_files_")

@pytest.fixture(scope="module")
def inputFileDirectory(systemTestFileDirectory):
    inputDirectory = os.path.join(systemTestFileDirectory.name, "input_files")
    os.mkdir(inputDirectory)
    return inputDirectory

@pytest.fixture(scope="module")
def expectedOutputDirectory(systemTestFileDirectory):
    expectedOutputDirectory = os.path.join(systemTestFileDirectory.name, "expected_output_files")
    os.mkdir(expectedOutputDirectory)
    return expectedOutputDirectory

@pytest.mark.skip(
    'Targeted reanalysis unit tests were very thorough - writing them at the system test level felt redundant. This may change with future updates.'
)
def test__targetedReanalysis__baseline_run(inputFileDirectory, expectedOutputDirectory):
    pass