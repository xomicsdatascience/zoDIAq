import numpy as np
from sklearn.metrics.pairwise import cosine_similarity

class Spectrum():
    def __init__(self,
                 mz: np.array = None,
                 intensity: np.array = None,
                 num_fragments: int = 10):
        """
        Creates a mass spectrometry spectrum, defined by its m/z bins and corresponding intensity values.
        Parameters
        ----------
        mz : np.array
            Array of mass:charge ratios.
        intensity : np.array
            Array of intensity values at corresponding indices of m/z.
        num_fragments : int
            Number of most intense fragments to extract.
        """
        self.mz = mz
        self.intensity = intensity
        self.extracted_intensity = None
        self.extracted_mz = None
        self._extract_most_intense(num_fragments=num_fragments)
        return

    def _extract_most_intense(self,
                             num_fragments: int = 10):
        """
        Extracts the most intense m/z bins, and stores the intensity and m/z sorted values in self.extracted_mz and
        self.extracted_intensity
        Parameters
        ----------
        spectrum_intensity : np.array
            Array containing the spectrum intensities to extract.
        spectrum_mz : np.array
            Array containing the m/z values corresponding to the spectrum intensities. If None, the index corresponding to
            the intensities is returned. Default: None.
        num_fragments : int
            Number of fragments to extract. Default: 10.

        Returns
        -------
        tuple
            Tuple of np.arrays sorted in decreasing order: (mz, intensity) corresponding to the most intense peaks. I.e.,
            the most intense peak is first.
        """
        if len(self.mz) == 0:
            self.mz = np.arange(len(self.intensity))

        # Sort fragments
        arg_idx = np.argsort(self.intensity)
        # Reverse list; take only top num_fragments
        # Check against np.inf
        if num_fragments == np.inf:
            # Take all fragments
            arg_idx = arg_idx[::-1]
        else:
            arg_idx = arg_idx[-1:-(num_fragments + 1):-1]
        self.extracted_mz = self.mz[arg_idx]
        self.extracted_intensity = self.intensity[arg_idx]
        return


    def get_matching_mz_indices(self,
                                spectrum_to_match: 'Spectrum',
                                match_tolerance_ppm: int = 30) -> tuple:
        """
        Checks whether the extracted fragments from this spectrum matches those of another. Matches are defined as two
         spectra having high intensity values at the same m/z bin, within tolerance.
        Parameters
        ----------
        spectrum_to_match : Spectrum
            Other Spectrum object to compare.
        min_fragment : int
            Minimum number of fragments to consider a match
        match_tolerance_ppm : int
            The margin of error given to bins.

        Returns
        -------
        tuple
            Pair of lists corresponding to the indices that are matched across the spectra. (self_idx, spectrum_to_match_idx)
        """

        # Spectra might be of different length; no clean way to compare them in bulk
        # First sort by mz, then crawl through each. Use argsort to restore original indices later
        self_arg_sorted = np.argsort(self.mz)
        sorted_self_mz = self.mz[self_arg_sorted]
        other_arg_sorted = np.argsort(spectrum_to_match.extracted_mz)
        sorted_other_mz = spectrum_to_match.extracted_mz[other_arg_sorted]

        self_start_idx = 0
        self_idx = []
        other_idx = []
        for other_mz_idx, other_mz in enumerate(sorted_other_mz):
            smallest_difference = np.inf
            for self_mz_idx, self_mz in zip(range(self_start_idx, len(sorted_self_mz)),
                                            sorted_self_mz[self_start_idx:]):
                # Since they're sorted, we can skip the ones we've already checked
                # Check whether we're within tolerance
                # (ref_val - val) * (1e6) / ref_val
                if np.isclose(self_mz, other_mz, rtol=match_tolerance_ppm):
                    # Mz values are sorted; once the difference increases we've found the minimum
                    if abs(self_mz - other_mz) < smallest_difference:
                        smallest_difference = abs(self_mz - other_mz)
                        closest_self_idx = self_mz_idx
                        closest_other_idx = other_mz_idx
                        continue
                    self_idx.append(closest_self_idx)
                    other_idx.append(closest_other_idx)
                    self_start_idx = self_mz_idx - 1  # Previous one might also be closest to the next bin
                    break
                # Check whether other_mz is larger; indicates we should increment self_mz
                elif self_mz > other_mz:
                    # Get next self_mz; ignore previously-checked mz from input spectra in next round
                    self_start_idx = self_mz_idx
                    break
        # Need to remap sorted idx to original idx:
        return self_arg_sorted[self_idx], other_arg_sorted[other_idx]

    def extracted_cosine_similarity(self,
                                    spectrum_to_compare: 'Spectrum') -> float:
        """
        Computes the cosine similarity between the extracted fragments of this spectrum and the input.

        Parameters
        ----------
        spectrum_to_compare : Spectrum
            Spectrum against which to compare
        Returns
        -------
        float
            Cosine similarity between this spectrum and the input.
        """

        # Get fragment idx
        self_fragment_idx, other_fragment_idx = self.get_matching_mz_indices(spectrum_to_match=spectrum_to_compare,
                                                                             match_tolerance_ppm=30)
        # if len(self_fragment_idx) == 0:
        #     return cosine_similarity(np.ones((1,1)), np.zeros((1,1)))
        return cosine_similarity(np.expand_dims(self.extracted_intensity[self_fragment_idx], axis=0),
                                 np.expand_dims(spectrum_to_compare.extracted_intensity[other_fragment_idx], axis=0))