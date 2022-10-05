import pyteomics.mzxml
from matplotlib import pyplot as plt
from csodiaq.spectrum import Spectrum


def spectrum_lineplot(spectrum: Spectrum,
                      positive: bool = True,
                      color: str = 'black',
                      label: str = None):
    """
    Plots the input Spectrum as a line plot
    Parameters
    ----------
    spectrum : Spectrum
        Spectrum object
    positive : bool
        Whether plot should have positive or negative intensities
    color : str
        Color for the line in the plot.
    label : str
        Label to use for the legend.

    Returns
    -------
    None
    """
    intensity = spectrum.intensity
    max_intensity = max(intensity)
    intensity = [i/max_intensity for i in intensity]
    if not positive:
        intensity = [-i for i in intensity]
    plt.vlines(spectrum.mz, [0]*len(intensity), intensity, colors=color, label=label)
    plt.hlines(0, min(spectrum.mz), max(spectrum.mz), colors='black')
    return


def spectrum_comparison_lineplot(spectrum_positive: Spectrum,
                                 spectrum_negative: Spectrum,
                                 colors: list = ('dodgerblue', 'purple'),
                                 labels: list = None):
    """
    Produces a plot with input Spectrum in positive/negative y directions.
    Parameters
    ----------
    spectrum_positive : Spectrum
        Spectrum to plot with positive intensity.
    spectrum_negative : Spectrum
        Spectrum to plot with negative intensity.
    colors : list
        List of colors to use for the plots.
    labels : list
        List of labels to use for the legends.
    Returns
    -------
    None
    """
    plt.figure()
    spectrum_lineplot(spectrum_positive, positive=True, color=colors[0], label=labels[0])
    spectrum_lineplot(spectrum_negative, positive=False, color=colors[1], label=labels[1])
    return