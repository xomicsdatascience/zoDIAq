from csodiaq.utils import format_protein_list_to_string
import pandas as pd
import numpy as np

def calculate_mz_of_heavy_version_of_peptide(peptide, lightMz, z):
    lightAndHeavyLysKMassDiff = 8.014199
    lightAndHeavyArgRMassDiff = 10.00827

    numLysK = peptide.count('K')
    numArgR = peptide.count('R')

    return (lightMz +
            (numLysK * lightAndHeavyLysKMassDiff) / z +
            (numArgR * lightAndHeavyArgRMassDiff) /  z)

def filter_to_only_keep_peptides_with_possibly_heavy_K_or_R_terminal_residue(fdrDf):
    return fdrDf[
        fdrDf["peptide"].str.endswith("R") | fdrDf["peptide"].str.endswith("K")
        ].reset_index(drop=True)

def filter_to_only_keep_top_peptides_unique_to_protein(fdrDf, maxPeptidesPerProtein):
    fdrDf = (
        fdrDf[fdrDf["uniquePeptide"] == 1]
        .sort_values("ionCount", ascending=False)
        .reset_index(drop=True)
    )
    return fdrDf.groupby(["leadingProtein"]).head(maxPeptidesPerProtein).reset_index(drop=True)

def calculate_mz_of_heavy_isotope_of_each_peptide(fdrDf):
    return fdrDf.apply(lambda x: calculate_mz_of_heavy_version_of_peptide(x["peptide"],x["MzLIB"],x["zLIB"]), axis=1)

def make_bin_assignments_for_mz_values(mzValues, binWidth=0.75):
    bins = np.arange(int(min(mzValues)) - 1, int(max(mzValues) + binWidth*3)+1, binWidth*2)
    values = np.digitize(mzValues, bins)
    mzBinValues = bins[values]
    mzBinValues -= binWidth
    return mzBinValues

def filter_out_peptides_based_on_user_settings(fdrDf, isIncludeHeavy, maximumPeptidesPerProtein):
    if isIncludeHeavy:
        fdrDf = filter_to_only_keep_peptides_with_possibly_heavy_K_or_R_terminal_residue(fdrDf)
    if maximumPeptidesPerProtein:
        fdrDf = filter_to_only_keep_top_peptides_unique_to_protein(fdrDf, maximumPeptidesPerProtein)
    return fdrDf

def calculate_binning_information_by_compensation_voltage(fdrDf, isIncludeHeavy):
    dfs = [
        organize_for_targeted_reanalysis_of_identified_peptides(df, isIncludeHeavy) for _, df in fdrDf.groupby(["CompensationVoltage"])
    ]
    return pd.concat(dfs)

def organize_for_targeted_reanalysis_of_identified_peptides(fdrDf, isIncludeHeavy):
    fdrDf["lightMzBin"] = make_bin_assignments_for_mz_values(fdrDf["MzLIB"])
    if isIncludeHeavy:
        fdrDf["heavyMz"] = calculate_mz_of_heavy_isotope_of_each_peptide(fdrDf)
        fdrDf["heavyMzBin"] = make_bin_assignments_for_mz_values(fdrDf["heavyMz"])
    return fdrDf

def make_targeted_reanalysis_line(peptide, mz, id):
    formula = ""
    adduct = "(no adduct)"
    genericCharge = 2
    return [peptide, formula, adduct, mz, genericCharge, id+1]

def consolidate_peptides_by_bin_values(df, isIncludeHeavy):
    bins = ["lightMzBin"]
    if isIncludeHeavy:
        bins.append("heavyMzBin")
    return df.groupby(bins).apply(lambda x: format_protein_list_to_string(x["peptide"])).reset_index(name="peptide").sort_values(bins)

def organize_binned_data_for_targeted_reanalysis(condensedDf, isIncludeHeavy):
    data = [
         make_targeted_reanalysis_line(condensedDf.loc[i]["peptide"], condensedDf.loc[i]["lightMzBin"], i) for i in range(len(condensedDf.index))
    ]
    if isIncludeHeavy:
        data.extend(
            [
                make_targeted_reanalysis_line(condensedDf.loc[i]["peptide"], condensedDf.loc[i]["heavyMzBin"], i) for i
                in range(len(condensedDf.index))

            ]
        )
    return data

def create_targeted_reanalysis_dataframe(df, isIncludeHeavy):
    consolidatedDf = consolidate_peptides_by_bin_values(df, isIncludeHeavy)
    targetedReanalysisData = organize_binned_data_for_targeted_reanalysis(consolidatedDf, isIncludeHeavy)
    return pd.DataFrame(targetedReanalysisData, columns=["Compound","Formula","Adduct","m.z","z","MSXID"]).sort_values(["MSXID","m.z"]).reset_index(drop=True)

def make_cv_header(cvValue):
    if cvValue:
        return f'CV_{str(abs(int(cvValue)))}'
    else:
        return 'noCV'

def create_targeted_reanalysis_dataframe_by_compensation_voltage(fdrDf, isIncludeHeavy):
    return {
        make_cv_header(cv[0]): create_targeted_reanalysis_dataframe(df, isIncludeHeavy) for cv, df in fdrDf.groupby(["CompensationVoltage"])
    }

def create_mass_spec_input_dataframes_for_targeted_reanalysis_of_identified_peptides(fdrDf, isIncludeHeavy=False, maximumPeptidesPerProtein=0):
    fdrDf = filter_out_peptides_based_on_user_settings(fdrDf, isIncludeHeavy, maximumPeptidesPerProtein)
    fdrDf = organize_for_targeted_reanalysis_of_identified_peptides(fdrDf, isIncludeHeavy)
    outputDfDict = create_targeted_reanalysis_dataframe_by_compensation_voltage(fdrDf, isIncludeHeavy)
    outputDfDict["fullDf"] = fdrDf
    return outputDfDict
