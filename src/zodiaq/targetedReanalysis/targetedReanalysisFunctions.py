from zodiaq.utils import format_protein_list_to_string
import pandas as pd
import numpy as np


def calculate_mz_of_heavy_version_of_peptide(peptide, lightMz, z):
    lightAndHeavyLysKMassDiff = 8.014199
    lightAndHeavyArgRMassDiff = 10.00827

    numLysK = peptide.count("K")
    numArgR = peptide.count("R")

    return (
        lightMz
        + (numLysK * lightAndHeavyLysKMassDiff) / z
        + (numArgR * lightAndHeavyArgRMassDiff) / z
    )


def filter_to_only_keep_peptides_with_possibly_heavy_K_or_R_terminal_residue(fdrDf):
    return fdrDf[
        fdrDf["peptide"].str.endswith("R") | fdrDf["peptide"].str.endswith("K")
    ].reset_index(drop=True)


def filter_to_only_keep_top_peptides_unique_to_protein(fdrDf, maxPeptidesPerProtein):
    fdrDf = (
        fdrDf[fdrDf["uniquePeptide"] == 1]
        .sort_values(["ionCount", "leadingProtein"], ascending=[False, True])
        .reset_index(drop=True)
    )
    return (
        fdrDf.groupby(["leadingProtein"])
        .head(maxPeptidesPerProtein)
        .reset_index(drop=True)
    )


def calculate_mz_of_heavy_isotope_of_each_peptide(fdrDf):
    return fdrDf.apply(
        lambda x: calculate_mz_of_heavy_version_of_peptide(
            x["peptide"], x["MzLIB"], x["zLIB"]
        ),
        axis=1,
    )


def make_bin_assignments_for_mz_values(mzValues, maxDistanceFromBinValue):
    fullBinWidth = maxDistanceFromBinValue * 2
    bins = np.arange(
        int(min(mzValues)) - 1, int(max(mzValues) + fullBinWidth) + 1, fullBinWidth
    )
    values = np.digitize(mzValues, bins)
    mzBinValues = bins[values]
    mzBinValues -= maxDistanceFromBinValue
    return mzBinValues


def check_if_is_protein_fdr_analysis(maximumPeptidesPerProtein):
    if maximumPeptidesPerProtein:
        return True
    return False


def filter_out_peptides_based_on_user_settings(
    fdrDf, isIncludeHeavyIsotopes, maximumPeptidesPerProtein
):
    if isIncludeHeavyIsotopes:
        fdrDf = (
            filter_to_only_keep_peptides_with_possibly_heavy_K_or_R_terminal_residue(
                fdrDf
            )
        )
    if check_if_is_protein_fdr_analysis(maximumPeptidesPerProtein):
        fdrDf = filter_to_only_keep_top_peptides_unique_to_protein(
            fdrDf, maximumPeptidesPerProtein
        )
    return fdrDf


def calculate_binning_information_by_compensation_voltage(
    fdrDf,
    isIncludeHeavyIsotopes,
    binValueProximity,
):
    dfs = [
        organize_for_targeted_reanalysis_of_identified_peptides(
            df, isIncludeHeavyIsotopes, binValueProximity
        )
        for _, df in fdrDf.groupby(["CompensationVoltage"])
    ]
    return pd.concat(dfs)


def organize_for_targeted_reanalysis_of_identified_peptides(
    fdrDf, isIncludeHeavyIsotopes, binValueProximity
):
    fdrDf["lightMzBin"] = make_bin_assignments_for_mz_values(
        fdrDf["MzLIB"], maxDistanceFromBinValue=binValueProximity
    )
    if isIncludeHeavyIsotopes:
        fdrDf["heavyMz"] = calculate_mz_of_heavy_isotope_of_each_peptide(fdrDf)
        fdrDf["heavyMzBin"] = make_bin_assignments_for_mz_values(
            fdrDf["heavyMz"], maxDistanceFromBinValue=binValueProximity
        )
    return fdrDf


def make_targeted_reanalysis_line(peptide, mz, id):
    formula = ""
    adduct = "(no adduct)"
    genericCharge = 2
    return [peptide, formula, adduct, mz, genericCharge, id + 1]


def consolidate_peptides_by_bin_values(df, isIncludeHeavyIsotopes):
    bins = ["lightMzBin"]
    if isIncludeHeavyIsotopes:
        bins.append("heavyMzBin")
    return (
        df.groupby(bins)
        .apply(lambda x: format_protein_list_to_string(x["peptide"]))
        .reset_index(name="peptide")
        .sort_values(bins)
    )


def organize_binned_data_for_targeted_reanalysis(condensedDf, isIncludeHeavyIsotopes):
    """
    Compiles data into a format that can be read by a mass spectrometer. Note that, in cases where the SILAC protocol
        is being used, the mass spectrometer can automatically match windows with related peptides provided they have
        same MSXID (the last column of the output). In this case the output will have paired rows with identical
        MSXIDs but differing m.z column values, one for the light mz value and one for the heavy mz value.
    """
    data = [
        make_targeted_reanalysis_line(
            condensedDf.loc[i]["peptide"], condensedDf.loc[i]["lightMzBin"], i
        )
        for i in range(len(condensedDf.index))
    ]
    if isIncludeHeavyIsotopes:
        data.extend(
            [
                make_targeted_reanalysis_line(
                    condensedDf.loc[i]["peptide"], condensedDf.loc[i]["heavyMzBin"], i
                )
                for i in range(len(condensedDf.index))
            ]
        )
    return data


def create_targeted_reanalysis_dataframe(df, isIncludeHeavyIsotopes):
    consolidatedDf = consolidate_peptides_by_bin_values(df, isIncludeHeavyIsotopes)
    targetedReanalysisData = organize_binned_data_for_targeted_reanalysis(
        consolidatedDf, isIncludeHeavyIsotopes
    )
    return (
        pd.DataFrame(
            targetedReanalysisData,
            columns=["Compound", "Formula", "Adduct", "m.z", "z", "MSXID"],
        )
        .sort_values(["MSXID", "m.z"])
        .reset_index(drop=True)
    )


def make_cv_header(cvValue):
    if cvValue:
        return f"CV_{str(abs(int(cvValue)))}"
    else:
        return "noCV"


def create_targeted_reanalysis_dataframes_by_compensation_voltage(
    fdrDf, isIncludeHeavyIsotopes
):
    return {
        make_cv_header(cv): create_targeted_reanalysis_dataframe(
            df, isIncludeHeavyIsotopes
        )
        for cv, df in fdrDf.groupby(["CompensationVoltage"])
    }


def create_mass_spec_input_dataframes_for_targeted_reanalysis_of_identified_peptides(
    fdrDf,
    isIncludeHeavyIsotopes,
    maximumPeptidesPerProtein,
    binValueProximity,
):
    """
    Creates dataframes with data to be fed into a mass spectrometer for targetted reanalysis of identified peptides.

    Extended Summary
    ----------------
    After peptide identification, additional experimentation is often required to verify or expand upon the
        significance of what was identified. For instance, some researchers would want to identify the relative
        quantity of peptides under different experimental conditions. Using Stable Isotope Labeling by Amino acids
        in Cell culture (SILAC), one can develop a cell culture where proteins use a heavier isotope of lysine and
        arginine. When targeting these peptides of interest, one can then compare the intensity of the "lighter"
        peptide isotope in one condition to the "heavier" peptide isotope in another, determining relative quantity
        changes between conditions. Thus, after identifying specific proteins and/or peptides using the identification
        workflow, one would need to prime the next mass spectrometer run to target those specific peptides only, one run
        for each category (light and heavy). This function generates tables that can be fed into a mass spectrometer
        to target specific peptides of interest.

    Note that as an aid to the researchers running the mass spectrometer, all targeted m/z values within a given
        m/z window are binned to reduce redundant window analysis.

    Also note that, for targeted reanalysis of data that utilized compensation voltage, a separate mass spectrometer
        input file must be generated for each. As such, a different dataset is made for each one. When no compensation
        voltage is provided, a single file will be generated. To keep these potentially variable compensation voltage
        tables separate, the return value is a dictionary highlighting the exact compensation voltage (or lack thereof)
        in the key name with the dataframe as the value. A summary file of all tables is also provided.

    Parameters
    ----------
    fdrDf : pandas DataFrame
    A dataframe containing peptides of interest. Note that this could be peptides belonging to specific identified
        proteins or just identified peptides generally, as indicated by the maximumPeptidesPerProtein variable.

    isIncludeHeavyIsotopes : boolean
    A boolean indicating if the experiment is to be primed for a SILAC-based mass spectrometry run.

    maximumPeptidesPerProtein : int
    The maximum number of peptides to be identified per protein of interest. In the case where peptides generally are
        being targetted this parameter is set to 0.

    Returns
    -------
    outputDfDict : dict
        A dictionary containing the targeted reanalysis files in addition to a summary file.
        'fullDf' : summary dataframe
        all other keys: target reanalysis dataframes, broken down by compensation voltage.

    """
    fdrDf = filter_out_peptides_based_on_user_settings(
        fdrDf, isIncludeHeavyIsotopes, maximumPeptidesPerProtein
    )
    fdrDf = organize_for_targeted_reanalysis_of_identified_peptides(
        fdrDf, isIncludeHeavyIsotopes, binValueProximity=binValueProximity
    )
    outputDfDict = create_targeted_reanalysis_dataframes_by_compensation_voltage(
        fdrDf, isIncludeHeavyIsotopes
    )
    outputDfDict["fullDf"] = fdrDf
    return outputDfDict
