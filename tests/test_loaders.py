from csodiaq.loaders import library
import unittest
from io import StringIO
import os
from os.path import join
"""
This file tests the loaders found in csodiaq/loaders. 
"""

pwd = os.path.dirname(__file__)
expected_cols = ['PrecursorMz', 'FullUniModPeptideName', 'PrecursorCharge', 'ProductMz', 'LibraryIntensity',
                 'transition_group_id', 'ProteinName']

class TestMGFLoader(unittest.TestCase):
    def test_default_load(self):
        file = os.path.join(pwd, 'test_files', 'sample_mgf.mgf')
        # Load file
        data = library.mgf_library_upload(file_name=file, max_peaks=10)

        # Assert file properties
        self.assertEqual(len(data), 3)
        k1 = (798.93, 'AGAAGAAAAAAAAAAAAAAGAK')  # expected entry from .mgf
        spec_dat = data[k1]
        self.assertEqual(spec_dat['PrecursorCharge'], 2)
        self.assertEqual(spec_dat['ProteinName'], '1/DECOY_0_sp|P55011|S12A2_HUMAN')
        self.assertEqual(len(spec_dat['Peaks']), 10)
        return

    def test_max_peaks(self):
        file = os.path.join(pwd, 'test_files', 'sample_mgf.mgf')

        # Load file
        for max_peaks in [5,10,20,50]:
            data = library.mgf_library_upload(file_name=file, max_peaks=max_peaks)
            k2 = (532.95, 'AAAAAAAAAAAAAAAGAGAGAK')  # expected entry from .mgf
            spec_dat = data[k2]
            self.assertEqual(len(spec_dat['Peaks']), max_peaks)
        return


class TestFragpipeLoader(unittest.TestCase):
    def test_default_load(self):
        print(join(pwd, 'test_files', 'fragpipe.tsv'))
        fragpipe = library.load_library(join(pwd, 'test_files', 'fragpipe.tsv'))
        # Confirm expected columns
        frag_cols = fragpipe.columns
        for col in expected_cols:
            if(col not in frag_cols):
                print(col)
            self.assertTrue(col in frag_cols)
        return