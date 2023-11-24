from zodiaq.loaders.library.libraryLoaderStrategy import (
    create_peaks_from_mz_intensity_lists_and_zodiaq_key_id,
    remove_low_intensity_peaks_below_max_peak_num,
    format_proteins_into_list_format,
)
import pytest


@pytest.fixture
def testPeaks():
    return [
        (0, 0, 0),
        (2, 1, 0),
        (4, 2, 0),
        (6, 3, 0),
        (8, 4, 0),
        (10, 5, 0),
        (12, 6, 0),
        (14, 7, 0),
        (1, 8, 0),
        (3, 9, 0),
        (5, 10, 0),
        (7, 11, 0),
        (9, 12, 0),
        (11, 13, 0),
        (13, 14, 0),
    ]


def test__library_loader_strategy__create_peaks_from_mz_intensity_lists_and_zodiaq_key_id(
    testPeaks,
):
    numPeaks = 15
    mzList = [i for i in range(0, numPeaks, 2)] + [i for i in range(1, numPeaks - 1, 2)]
    intensityList = [i for i in range(numPeaks)]
    id = 0
    peaks = create_peaks_from_mz_intensity_lists_and_zodiaq_key_id(
        mzList, intensityList, id
    )
    assert peaks == testPeaks


def test__library_loader_strategy__remove_low_intensity_peaks_below_max_peak_num(
    testPeaks,
):
    maxPeakNum = 10
    expectedReducedTestPeaks = [
        (13, 14, 0),
        (11, 13, 0),
        (9, 12, 0),
        (7, 11, 0),
        (5, 10, 0),
        (3, 9, 0),
        (1, 8, 0),
        (14, 7, 0),
        (12, 6, 0),
        (10, 5, 0),
    ]
    reducedTestPeaks = remove_low_intensity_peaks_below_max_peak_num(
        testPeaks, maxPeakNum
    )
    assert reducedTestPeaks == expectedReducedTestPeaks


def test__library_loader_strategy__remove_low_intensity_peaks_below_max_peak_num__all_peaks_returned_when_length_fewer_than_max_peak_num(
    testPeaks,
):
    maxPeakNum = 10
    expectedReducedShortTestPeaks = [
        (1, 8, 0),
        (14, 7, 0),
        (12, 6, 0),
        (10, 5, 0),
        (8, 4, 0),
        (6, 3, 0),
        (4, 2, 0),
        (2, 1, 0),
        (0, 0, 0),
    ]
    shortTestPeaks = testPeaks[: maxPeakNum - 1]
    reducedShortTestPeaks = remove_low_intensity_peaks_below_max_peak_num(
        shortTestPeaks, maxPeakNum
    )
    assert reducedShortTestPeaks == expectedReducedShortTestPeaks


def test__library_loader_strategy__format_proteins_into_list_format():
    zodiaqLibraryDict = {
        (100.0, "PEPTIDE"): {
            "proteinName": "PROTEIN",
        }
    }
    format_proteins_into_list_format(zodiaqLibraryDict)
    assert zodiaqLibraryDict[(100.0, "PEPTIDE")]["proteinName"] == "1/PROTEIN"

    zodiaqLibraryDict = {
        (200.0, "PEPTIDE"): {
            "proteinName": "",
        }
    }
    format_proteins_into_list_format(zodiaqLibraryDict)
    assert zodiaqLibraryDict[(200.0, "PEPTIDE")]["proteinName"] == ""
