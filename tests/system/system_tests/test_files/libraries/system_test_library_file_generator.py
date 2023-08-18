import pandas as pd
import os

def get_parent_dir():
    return os.path.dirname(os.path.abspath(__file__))

minMz = 400
maxMz = 1600
precursorMzDiff = 1
peakMzDiff = 0.1
peakIntensityDiff = 1000.0
data = []
for precursorMz in range(minMz, maxMz):
    if precursorMz < (maxMz+minMz) / 2:
        protein = f'1/protein{precursorMz}'
    else:
        protein = f'1/DECOY_protein{precursorMz}'
    variableDict = {
        'precursorMz': precursorMz,
        'peptide': f'peptide{precursorMz}',
        'precursorCharge': 1,
        'id': f'id{precursorMz}',
        'protein': protein
    }
    for i in range(1, 11):
        peakMz = precursorMz - (precursorMzDiff/2) + (i-1)*peakMzDiff
        peakIntensity = i * peakIntensityDiff
        data.append([
            variableDict["precursorMz"],
            variableDict["peptide"],
            variableDict["precursorCharge"],
            variableDict["id"],
            variableDict["protein"],
            peakMz,
            peakIntensity,
        ])

df = pd.DataFrame(data, columns=['precursorMz', 'peptide', 'precursorCharge', 'id', 'protein', 'peakMz', 'peakIntensity'])
df.to_csv(os.path.join(get_parent_dir(), 'template_test_library.csv'), index=False)
df.columns = ["PrecursorMz", "FullUniModPeptideName", "PrecursorCharge", "transition_group_id", "ProteinName", "ProductMz", "LibraryIntensity"]
df.to_csv(os.path.join(get_parent_dir(), 'spectrast_test_library.csv'), index=False)
