import numpy as np
from numba import njit, typeof
import matplotlib.pyplot as pyplot
import csv

class quantSpectraMatcher:
    def __init__(self):
        self.libraryTags = np.array([],dtype=int)
        self.queryTags = np.array([],dtype=int)
        self.libraryIntensities = np.array([],dtype=float)
        self.queryIntensities = np.array([],dtype=float)
        self.ppmMatches = np.array([],dtype=float)
        self.scores = np.array([],dtype=float)
        self.decoys = np.array([],dtype=float)
