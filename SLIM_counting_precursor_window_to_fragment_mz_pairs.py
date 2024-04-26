from pyteomics import mzxml
import sys

filePath = sys.argv[1]

pairs = []
with mzxml.read(filePath) as spectra:
    for spec in spectra:
        if "precursorMz" not in spec: # skip if an MS1 spectra
            continue
        precMz = spec["precursorMz"][0]["precursorMz"]
        windowWidth = spec["precursorMz"][0]["windowWideness"]
        for fragmentMz in spec["m/z array"]:
            pairs.append((precMz, windowWidth, fragmentMz))

print(f'total pairs: {len(pairs)}')
print(f'unique pairs: {len(set(pairs))}')
