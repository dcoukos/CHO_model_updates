#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 19 23:45:56 2018

@author: dimitricoukos
"""
import unittest
import json
import DataTreatment
from DataTreatment import openJson, write

class SampleData(unittest.TestCase):
    initial_input = {
    "GLNLASEer": {
        "N-octanoyl-DL-homoserine lactone": [],
        "5-butyl-4-methyldihydro-2(3H)-furanone": [],
        "gamma-undecanolactone": [
            {
                "wild-type": False,
                "organism": "Sulfolobus solfataricus",
                "turnoverNumber": "3.92",
                "ecNumber": "3.1.1.25"
            },
            {
                "wild-type": False,
                "organism": "Sulfolobus solfataricus",
                "turnoverNumber": "4.25",
                "ecNumber": "3.1.1.25"
            },
            {
                "wild-type": False,
                "organism": "Sulfolobus solfataricus",
                "turnoverNumber": "4.55",
                "ecNumber": "3.1.1.25"
            },
            {
                "wild-type": False,
                "organism": "Sulfolobus solfataricus",
                "turnoverNumber": "4.63",
                "ecNumber": "3.1.1.25"
            },
            {
                "wild-type": True,
                "organism": "Sulfolobus solfataricus",
                "turnoverNumber": "4.95",
                "ecNumber": "3.1.1.25"
            },
            {
                "wild-type": False,
                "organism": "Sulfolobus solfataricus",
                "turnoverNumber": "5.64",
                "ecNumber": "3.1.1.25"
            }
        ],
        "gamma-dodecanolactone": [],
        "N-(3-oxododecanoyl)-L-homoserine lactone": [
            {
                "wild-type": True,
                "organism": "Sulfolobus solfataricus",
                "turnoverNumber": "1.01",
                "ecNumber": "3.1.1.25"
            },
            {
                "wild-type": False,
                "organism": "Sulfolobus solfataricus",
                "turnoverNumber": "1.8",
                "ecNumber": "3.1.1.25"
            },
            {
                "wild-type": False,
                "organism": "Sulfolobus solfataricus",
                "turnoverNumber": "3",
                "ecNumber": "3.1.1.25"
            },
            {
                "wild-type": False,
                "organism": "Sulfolobus solfataricus",
                "turnoverNumber": "6.44",
                "ecNumber": "3.1.1.25"
            }
        ],
        "nonanoic-1,5-lactone": [],
        "gamma-dodecalactone": [],
        "N-(3-oxodecanoyl)-L-homoserine lactone": [
            {
                "wild-type": False,
                "organism": "Sulfolobus solfataricus",
                "turnoverNumber": "0.19",
                "ecNumber": "3.1.1.25"
            },
            {
                "wild-type": False,
                "organism": "Sulfolobus solfataricus",
                "turnoverNumber": "0.6",
                "ecNumber": "3.1.1.25"
            },
            {
                "wild-type": False,
                "organism": "Sulfolobus solfataricus",
                "turnoverNumber": "3.96",
                "ecNumber": "3.1.1.25"
            },
            {
                "wild-type": True,
                "organism": "Sulfolobus solfataricus",
                "turnoverNumber": "4.52",
                "ecNumber": "3.1.1.25"
            }
        ],
        "gamma-dodecanoic lactone": [
            {
                "organism": "Homo sapiens",
                "turnoverNumber": "101",
                "ecNumber": "3.1.1.25"
            }
        ],
        "gamma-heptalactone": [],
        "undecanoic-gamma-lactone": [],
        "N-(2-oxotetrahydrofuran-3-yl)pentanamide": [],
        "N-octanoylhomoserine lactone": [],
        "nonanoic-gamma-lactone": [
            {
                "wild-type": False,
                "organism": "Sulfolobus islandicus",
                "turnoverNumber": "2",
                "ecNumber": "3.1.1.25"
            },
            {
                "wild-type": True,
                "organism": "Sulfolobus islandicus",
                "turnoverNumber": "3.1",
                "ecNumber": "3.1.1.25"
            }
        ],
        "5-(thiobutyl)butyrolactone": [
            {
                "wild-type": False,
                "organism": "Homo sapiens",
                "turnoverNumber": "7.5",
                "ecNumber": "3.1.1.25"
            },
            {
                "wild-type": False,
                "organism": "Homo sapiens",
                "turnoverNumber": "19.4",
                "ecNumber": "3.1.1.25"
            },
            {
                "wild-type": True,
                "organism": "Homo sapiens",
                "turnoverNumber": "116",
                "ecNumber": "3.1.1.25"
            }
        ],
        "N-hexanoylhomoserine lactone": [],
        "N-(3-oxodecanoyl)-DL-homoserine lactone": [],
        "delta-undecalactone": [],
        "delta-dodecalactone": [],
        "gamma-(S)-valerolactone": [],
        "gamma-undecalactone": [],
        "gamma-(R)-valerolactone": [],
        "octanoyl-L-homoserine lactone": [],
        "N-(3-oxododecanoyl)-DL-homoserine lactone": [],
        "gamma-(S)-caprolactone": [],
        "dodecanoic-1,5-lactone": [],
        "gamma-nonanoic acid lactone": [],
        "gamma-heptanolactone": [],
        "Paraoxon": [
            {
                "wild-type": False,
                "organism": "Sulfolobus solfataricus",
                "turnoverNumber": "8.47",
                "ecNumber": "3.1.1.25"
            },
            {
                "wild-type": True,
                "organism": "Sulfolobus solfataricus",
                "turnoverNumber": "12.6",
                "ecNumber": "3.1.1.25"
            }
        ],
        "dodecanoic-gamma-lactone": [],
        "undecanoic-1,5-lactone": [],
        "gamma-heptanolide": [
            {
                "organism": "Sulfolobus acidocaldarius",
                "turnoverNumber": "10.25",
                "ecNumber": "3.1.1.25"
            },
            {
                "organism": "Homo sapiens",
                "turnoverNumber": "34",
                "ecNumber": "3.1.1.25"
            }
        ],
        "delta-undecanolactone": [
            {
                "wild-type": True,
                "organism": "Sulfolobus solfataricus",
                "turnoverNumber": "12.65",
                "ecNumber": "3.1.1.25"
            },
            {
                "wild-type": False,
                "organism": "Sulfolobus solfataricus",
                "turnoverNumber": "44.8",
                "ecNumber": "3.1.1.25"
            },
            {
                "wild-type": False,
                "organism": "Sulfolobus solfataricus",
                "turnoverNumber": "56.8",
                "ecNumber": "3.1.1.25"
            },
            {
                "wild-type": False,
                "organism": "Sulfolobus solfataricus",
                "turnoverNumber": "58",
                "ecNumber": "3.1.1.25"
            },
            {
                "wild-type": False,
                "organism": "Sulfolobus solfataricus",
                "turnoverNumber": "66.5",
                "ecNumber": "3.1.1.25"
            },
            {
                "wild-type": False,
                "organism": "Sulfolobus solfataricus",
                "turnoverNumber": "71.2",
                "ecNumber": "3.1.1.25"
            },
            {
                "wild-type": False,
                "organism": "Sulfolobus solfataricus",
                "turnoverNumber": "93.3",
                "ecNumber": "3.1.1.25"
            }
        ],
        "gamma-nonalactone": [
            {
                "wild-type": True,
                "organism": "Sulfolobus solfataricus",
                "turnoverNumber": "5.54",
                "ecNumber": "3.1.1.25"
            },
            {
                "wild-type": False,
                "organism": "Homo sapiens",
                "turnoverNumber": "5.57",
                "ecNumber": "3.1.1.25"
            },
            {
                "wild-type": True,
                "organism": "Homo sapiens",
                "turnoverNumber": "31",
                "ecNumber": "3.1.1.25"
            },
            {
                "wild-type": True,
                "organism": "Vulcanisaeta moutnovskia",
                "turnoverNumber": "44.49",
                "ecNumber": "3.1.1.25"
            }
        ],
        "N-(3-oxohexanoyl)-L-homoserine lactone": [],
        "N-(3-oxooctanoyl)-L-homoserine lactone": [],
        "3-oxo-octanoyl-L-homoserine lactone": [],
        "gamma-dodecanoic acid lactone": [],
        "gamma-(R)-caprolactone": [],
        "4-methoxy phenyl acetate": [],
        "epsilon-caprolactone": [
            {
                "wild-type": True,
                "organism": "Sulfolobus islandicus",
                "turnoverNumber": "7.27",
                "ecNumber": "3.1.1.25"
            },
            {
                "wild-type": True,
                "organism": "Sulfolobus acidocaldarius",
                "turnoverNumber": "15.04",
                "ecNumber": "3.1.1.25"
            }
        ],
        "Gamma-caprolactone": [
            {
                "wild-type": False,
                "organism": "Homo sapiens",
                "turnoverNumber": "25",
                "ecNumber": "3.1.1.25"
            },
            {
                "wild-type": False,
                "organism": "Homo sapiens",
                "turnoverNumber": "44",
                "ecNumber": "3.1.1.25"
            },
            {
                "wild-type": True,
                "organism": "Homo sapiens",
                "turnoverNumber": "44",
                "ecNumber": "3.1.1.25"
            },
            {
                "wild-type": True,
                "organism": "Vulcanisaeta moutnovskia",
                "turnoverNumber": "112.3",
                "ecNumber": "3.1.1.25"
            }
        ],
        "gamma-butyrolactone": [
            {
                "wild-type": True,
                "organism": "Sulfolobus islandicus",
                "turnoverNumber": "5.75",
                "ecNumber": "3.1.1.25"
            },
            {
                "wild-type": True,
                "organism": "Homo sapiens",
                "turnoverNumber": "111",
                "ecNumber": "3.1.1.25"
            },
            {
                "wild-type": True,
                "organism": "Homo sapiens",
                "turnoverNumber": "111",
                "ecNumber": "3.1.1.25"
            }
        ],
        "delta-valerolactone": [
            {
                "wild-type": False,
                "organism": "Homo sapiens",
                "turnoverNumber": "0.5",
                "ecNumber": "3.1.1.25"
            },
            {
                "wild-type": False,
                "organism": "Homo sapiens",
                "turnoverNumber": "0.9",
                "ecNumber": "3.1.1.25"
            },
            {
                "wild-type": True,
                "organism": "Homo sapiens",
                "turnoverNumber": "29.8",
                "ecNumber": "3.1.1.25"
            },
            {
                "wild-type": False,
                "organism": "Homo sapiens",
                "turnoverNumber": "40",
                "ecNumber": "3.1.1.25"
            },
            {
                "wild-type": False,
                "organism": "Homo sapiens",
                "turnoverNumber": "69.4",
                "ecNumber": "3.1.1.25"
            },
            {
                "wild-type": False,
                "organism": "Homo sapiens",
                "turnoverNumber": "94",
                "ecNumber": "3.1.1.25"
            },
            {
                "wild-type": False,
                "organism": "Homo sapiens",
                "turnoverNumber": "156",
                "ecNumber": "3.1.1.25"
            },
            {
                "wild-type": False,
                "organism": "Homo sapiens",
                "turnoverNumber": "210",
                "ecNumber": "3.1.1.25"
            },
            {
                "wild-type": True,
                "organism": "Homo sapiens",
                "turnoverNumber": "210",
                "ecNumber": "3.1.1.25"
            },
            {
                "wild-type": True,
                "organism": "Homo sapiens",
                "turnoverNumber": "210",
                "ecNumber": "3.1.1.25"
            },
            {
                "wild-type": True,
                "organism": "Homo sapiens",
                "turnoverNumber": "632",
                "ecNumber": "3.1.1.25"
            }
        ],
        "gamma-undecanoiclactone": [],
        "9-oxo-N-(2-oxotetrahydrofuran-3-yl)undecanamide": [],
        "N-(3-oxooctanoyl)-DL-homoserine lactone": [
            {
                "wild-type": False,
                "organism": "Sulfolobus islandicus",
                "turnoverNumber": "0.92",
                "ecNumber": "3.1.1.25"
            },
            {
                "wild-type": False,
                "organism": "Sulfolobus islandicus",
                "turnoverNumber": "0.97",
                "ecNumber": "3.1.1.25"
            },
            {
                "wild-type": True,
                "organism": "Sulfolobus islandicus",
                "turnoverNumber": "4.1",
                "ecNumber": "3.1.1.25"
            }
        ],
        "N-dodecanoylhomoserine lactone": [],
        "nonanoic-delta-lactone": [],
        "7-oxo-N-(2-oxotetrahydrofuran-3-yl)nonanamide": [],
        "dodecanoic-delta-lactone": [],
        "dihydrocoumarin": [
            {
                "organism": "Homo sapiens",
                "turnoverNumber": "152",
                "ecNumber": "3.1.1.25"
            }
        ],
        "N-dodecanoyl-DL-homoserine lactone": [],
        "dodecanoic-1,4-lactone": [],
        "gamma-undecanoic acid lactone": [],
        "delta-nonalactone": [
            {
                "organism": "Homo sapiens",
                "turnoverNumber": "48",
                "ecNumber": "3.1.1.25"
            },
            {
                "organism": "Vulcanisaeta moutnovskia",
                "turnoverNumber": "88.91",
                "ecNumber": "3.1.1.25"
            }
        ],
        "undecanoic-1,4-lactone": [],
        "pantoyl lactone": [],
        "nonanoic-1,4-lactone": [],
        "N-(3-oxohexanoyl)homoserine lactone": [],
        "undecanoic-delta-lactone": [
            {
                "wild-type": False,
                "organism": "Sulfolobus islandicus",
                "turnoverNumber": "12.9",
                "ecNumber": "3.1.1.25"
            },
            {
                "wild-type": False,
                "organism": "Sulfolobus islandicus",
                "turnoverNumber": "14.1",
                "ecNumber": "3.1.1.25"
            },
            {
                "wild-type": True,
                "organism": "Sulfolobus islandicus",
                "turnoverNumber": "17.65",
                "ecNumber": "3.1.1.25"
            }
        ],
        "3-oxo-decanoyl-L-homoserine lactone": [],
        "N-(3-oxooctanoyl)homoserine lactone": []
    },
    "CYSTS": {
        "L-Ser": [],
        "homocysteine": [
            {
                "wild-type": True,
                "organism": "Homo sapiens",
                "turnoverNumber": "6.2",
                "ecNumber": "4.2.1.22"
            },
            {
                "wild-type": True,
                "organism": "Saccharomyces cerevisiae",
                "turnoverNumber": "7.38",
                "ecNumber": "4.2.1.22"
            },
            {
                "wild-type": False,
                "organism": "Homo sapiens",
                "turnoverNumber": "15.5",
                "ecNumber": "4.2.1.22"
            },
            {
                "wild-type": False,
                "organism": "Homo sapiens",
                "turnoverNumber": "32.1",
                "ecNumber": "4.2.1.22"
            },
            {
                "wild-type": True,
                "organism": "Homo sapiens",
                "turnoverNumber": "34",
                "ecNumber": "4.2.1.22"
            }
        ],
        "L-homocysteine": [
            {
                "wild-type": True,
                "organism": "Homo sapiens",
                "turnoverNumber": "0.031",
                "ecNumber": "4.2.1.22"
            },
            {
                "wild-type": False,
                "organism": "Homo sapiens",
                "turnoverNumber": "0.04",
                "ecNumber": "4.2.1.22"
            },
            {
                "wild-type": False,
                "organism": "Homo sapiens",
                "turnoverNumber": "0.09",
                "ecNumber": "4.2.1.22"
            },
            {
                "wild-type": False,
                "organism": "Saccharomyces cerevisiae",
                "turnoverNumber": "0.85",
                "ecNumber": "4.2.1.22"
            },
            {
                "wild-type": True,
                "organism": "Homo sapiens",
                "turnoverNumber": "3.3",
                "ecNumber": "4.2.1.22"
            },
            {
                "wild-type": True,
                "organism": "Homo sapiens",
                "turnoverNumber": "4.66",
                "ecNumber": "4.2.1.22"
            },
            {
                "wild-type": True,
                "organism": "Homo sapiens",
                "turnoverNumber": "7.93",
                "ecNumber": "4.2.1.22"
            },
            {
                "wild-type": False,
                "organism": "Homo sapiens",
                "turnoverNumber": "9.06",
                "ecNumber": "4.2.1.22"
            },
            {
                "wild-type": True,
                "organism": "Homo sapiens",
                "turnoverNumber": "12.7",
                "ecNumber": "4.2.1.22"
            },
            {
                "wild-type": True,
                "organism": "Saccharomyces cerevisiae",
                "turnoverNumber": "17",
                "ecNumber": "4.2.1.22"
            },
            {
                "wild-type": True,
                "organism": "Saccharomyces cerevisiae",
                "turnoverNumber": "21.5",
                "ecNumber": "4.2.1.22"
            }
        ],
        "L-cystathionine": [
            {
                "wild-type": False,
                "organism": "Saccharomyces cerevisiae",
                "turnoverNumber": "0.083",
                "ecNumber": "4.2.1.22"
            },
            {
                "wild-type": False,
                "organism": "Saccharomyces cerevisiae",
                "turnoverNumber": "0.133",
                "ecNumber": "4.2.1.22"
            },
            {
                "wild-type": False,
                "organism": "Saccharomyces cerevisiae",
                "turnoverNumber": "0.418",
                "ecNumber": "4.2.1.22"
            },
            {
                "wild-type": False,
                "organism": "Saccharomyces cerevisiae",
                "turnoverNumber": "0.56",
                "ecNumber": "4.2.1.22"
            },
            {
                "wild-type": True,
                "organism": "Saccharomyces cerevisiae",
                "turnoverNumber": "0.56",
                "ecNumber": "4.2.1.22"
            },
            {
                "wild-type": True,
                "organism": "Saccharomyces cerevisiae",
                "turnoverNumber": "1.03",
                "ecNumber": "4.2.1.22"
            },
            {
                "wild-type": True,
                "organism": "Saccharomyces cerevisiae",
                "turnoverNumber": "6.08",
                "ecNumber": "4.2.1.22"
            },
            {
                "wild-type": True,
                "organism": "Saccharomyces cerevisiae",
                "turnoverNumber": "6.08",
                "ecNumber": "4.2.1.22"
            }
        ],
        "L-cysteine": [
            {
                "wild-type": False,
                "organism": "Homo sapiens",
                "turnoverNumber": "1.95",
                "ecNumber": "4.2.1.22"
            },
            {
                "wild-type": False,
                "organism": "Homo sapiens",
                "turnoverNumber": "3.13",
                "ecNumber": "4.2.1.22"
            },
            {
                "wild-type": True,
                "organism": "Homo sapiens",
                "turnoverNumber": "3.13",
                "ecNumber": "4.2.1.22"
            },
            {
                "wild-type": True,
                "organism": "Homo sapiens",
                "turnoverNumber": "4.39",
                "ecNumber": "4.2.1.22"
            },
            {
                "wild-type": True,
                "organism": "Homo sapiens",
                "turnoverNumber": "4.39",
                "ecNumber": "4.2.1.22"
            }
        ],
        "L-serine": [
            {
                "wild-type": False,
                "organism": "Saccharomyces cerevisiae",
                "turnoverNumber": "0.082",
                "ecNumber": "4.2.1.22"
            },
            {
                "wild-type": False,
                "organism": "Homo sapiens",
                "turnoverNumber": "0.15",
                "ecNumber": "4.2.1.22"
            },
            {
                "wild-type": False,
                "organism": "Saccharomyces cerevisiae",
                "turnoverNumber": "0.45",
                "ecNumber": "4.2.1.22"
            },
            {
                "wild-type": True,
                "organism": "Homo sapiens",
                "turnoverNumber": "0.52",
                "ecNumber": "4.2.1.22"
            },
            {
                "wild-type": False,
                "organism": "Saccharomyces cerevisiae",
                "turnoverNumber": "0.85",
                "ecNumber": "4.2.1.22"
            },
            {
                "wild-type": False,
                "organism": "Homo sapiens",
                "turnoverNumber": "1.3",
                "ecNumber": "4.2.1.22"
            },
            {
                "wild-type": False,
                "organism": "Homo sapiens",
                "turnoverNumber": "1.67",
                "ecNumber": "4.2.1.22"
            },
            {
                "wild-type": False,
                "organism": "Homo sapiens",
                "turnoverNumber": "2.5",
                "ecNumber": "4.2.1.22"
            },
            {
                "wild-type": False,
                "organism": "Homo sapiens",
                "turnoverNumber": "2.9",
                "ecNumber": "4.2.1.22"
            },
            {
                "wild-type": True,
                "organism": "Homo sapiens",
                "turnoverNumber": "3.67",
                "ecNumber": "4.2.1.22"
            },
            {
                "wild-type": False,
                "organism": "Saccharomyces cerevisiae",
                "turnoverNumber": "5.3",
                "ecNumber": "4.2.1.22"
            },
            {
                "wild-type": True,
                "organism": "Homo sapiens",
                "turnoverNumber": "5.4",
                "ecNumber": "4.2.1.22"
            },
            {
                "wild-type": False,
                "organism": "Homo sapiens",
                "turnoverNumber": "5.9",
                "ecNumber": "4.2.1.22"
            },
            {
                "wild-type": False,
                "organism": "Homo sapiens",
                "turnoverNumber": "7.5",
                "ecNumber": "4.2.1.22"
            },
            {
                "wild-type": False,
                "organism": "Saccharomyces cerevisiae",
                "turnoverNumber": "7.6",
                "ecNumber": "4.2.1.22"
            },
            {
                "wild-type": False,
                "organism": "Saccharomyces cerevisiae",
                "turnoverNumber": "8.2",
                "ecNumber": "4.2.1.22"
            },
            {
                "wild-type": True,
                "organism": "Homo sapiens",
                "turnoverNumber": "10.19",
                "ecNumber": "4.2.1.22"
            },
            {
                "wild-type": True,
                "organism": "Homo sapiens",
                "turnoverNumber": "10.2",
                "ecNumber": "4.2.1.22"
            },
            {
                "wild-type": False,
                "organism": "Saccharomyces cerevisiae",
                "turnoverNumber": "13.2",
                "ecNumber": "4.2.1.22"
            },
            {
                "wild-type": True,
                "organism": "Homo sapiens",
                "turnoverNumber": "14.01",
                "ecNumber": "4.2.1.22"
            },
            {
                "wild-type": False,
                "organism": "Homo sapiens",
                "turnoverNumber": "14.6",
                "ecNumber": "4.2.1.22"
            },
            {
                "wild-type": False,
                "organism": "Saccharomyces cerevisiae",
                "turnoverNumber": "14.7",
                "ecNumber": "4.2.1.22"
            },
            {
                "wild-type": False,
                "organism": "Homo sapiens",
                "turnoverNumber": "15.8",
                "ecNumber": "4.2.1.22"
            },
            {
                "wild-type": False,
                "organism": "Saccharomyces cerevisiae",
                "turnoverNumber": "16.8",
                "ecNumber": "4.2.1.22"
            },
            {
                "wild-type": True,
                "organism": "Saccharomyces cerevisiae",
                "turnoverNumber": "17",
                "ecNumber": "4.2.1.22"
            },
            {
                "wild-type": False,
                "organism": "Homo sapiens",
                "turnoverNumber": "19",
                "ecNumber": "4.2.1.22"
            },
            {
                "wild-type": True,
                "organism": "Homo sapiens",
                "turnoverNumber": "19.7",
                "ecNumber": "4.2.1.22"
            },
            {
                "wild-type": True,
                "organism": "Homo sapiens",
                "turnoverNumber": "21",
                "ecNumber": "4.2.1.22"
            },
            {
                "wild-type": True,
                "organism": "Saccharomyces cerevisiae",
                "turnoverNumber": "21.5",
                "ecNumber": "4.2.1.22"
            },
            {
                "wild-type": True,
                "organism": "Homo sapiens",
                "turnoverNumber": "39",
                "ecNumber": "4.2.1.22"
            },
            {
                "wild-type": False,
                "organism": "Homo sapiens",
                "turnoverNumber": "45",
                "ecNumber": "4.2.1.22"
            }
        ],
        "more": []
    },
    "BTNDe": {},
    "BTNDm": {},
    "GTHPe": {
        "cumene peroxide": [
            {
                "organism": "Lucilia cuprina",
                "turnoverNumber": "35.78",
                "ecNumber": "1.11.1.9"
            }
        ],
        "GSH": [
            {
                "organism": "Homo sapiens",
                "turnoverNumber": "19.5",
                "ecNumber": "1.11.1.9"
            },
            {
                "organism": "Homo sapiens",
                "turnoverNumber": "24.5",
                "ecNumber": "1.11.1.9"
            },
            {
                "wild-type": True,
                "organism": "Homo sapiens",
                "turnoverNumber": "221.7",
                "ecNumber": "1.11.1.9"
            },
            {
                "wild-type": True,
                "organism": "Homo sapiens",
                "turnoverNumber": "293.3",
                "ecNumber": "1.11.1.9"
            },
            {
                "wild-type": True,
                "organism": "Homo sapiens",
                "turnoverNumber": "361.7",
                "ecNumber": "1.11.1.9"
            },
            {
                "wild-type": True,
                "organism": "Homo sapiens",
                "turnoverNumber": "408.3",
                "ecNumber": "1.11.1.9"
            }
        ],
        "tert-butyl hydroperoxide": [],
        "H2O2": [
            {
                "organism": "Homo sapiens",
                "turnoverNumber": "10.55",
                "ecNumber": "1.11.1.9"
            },
            {
                "organism": "Homo sapiens",
                "turnoverNumber": "16.02",
                "ecNumber": "1.11.1.9"
            },
            {
                "organism": "Homo sapiens",
                "turnoverNumber": "20.83",
                "ecNumber": "1.11.1.9"
            },
            {
                "organism": "Lucilia cuprina",
                "turnoverNumber": "44.03",
                "ecNumber": "1.11.1.9"
            },
            {
                "wild-type": True,
                "organism": "Homo sapiens",
                "turnoverNumber": "316.7",
                "ecNumber": "1.11.1.9"
            },
            {
                "wild-type": True,
                "organism": "Homo sapiens",
                "turnoverNumber": "445",
                "ecNumber": "1.11.1.9"
            },
            {
                "wild-type": True,
                "organism": "Homo sapiens",
                "turnoverNumber": "560",
                "ecNumber": "1.11.1.9"
            },
            {
                "wild-type": True,
                "organism": "Homo sapiens",
                "turnoverNumber": "670",
                "ecNumber": "1.11.1.9"
            }
        ]
    },
    "GTHPm": {
        "cumene peroxide": [
            {
                "organism": "Lucilia cuprina",
                "turnoverNumber": "35.78",
                "ecNumber": "1.11.1.9"
            }
        ],
        "GSH": [
            {
                "organism": "Homo sapiens",
                "turnoverNumber": "19.5",
                "ecNumber": "1.11.1.9"
            },
            {
                "organism": "Homo sapiens",
                "turnoverNumber": "24.5",
                "ecNumber": "1.11.1.9"
            },
            {
                "wild-type": True,
                "organism": "Homo sapiens",
                "turnoverNumber": "221.7",
                "ecNumber": "1.11.1.9"
            },
            {
                "wild-type": True,
                "organism": "Homo sapiens",
                "turnoverNumber": "293.3",
                "ecNumber": "1.11.1.9"
            },
            {
                "wild-type": True,
                "organism": "Homo sapiens",
                "turnoverNumber": "361.7",
                "ecNumber": "1.11.1.9"
            },
            {
                "wild-type": True,
                "organism": "Homo sapiens",
                "turnoverNumber": "408.3",
                "ecNumber": "1.11.1.9"
            }
        ],
        "tert-butyl hydroperoxide": [],
        "H2O2": [
            {
                "organism": "Homo sapiens",
                "turnoverNumber": "10.55",
                "ecNumber": "1.11.1.9"
            },
            {
                "organism": "Homo sapiens",
                "turnoverNumber": "16.02",
                "ecNumber": "1.11.1.9"
            },
            {
                "organism": "Homo sapiens",
                "turnoverNumber": "20.83",
                "ecNumber": "1.11.1.9"
            },
            {
                "organism": "Lucilia cuprina",
                "turnoverNumber": "44.03",
                "ecNumber": "1.11.1.9"
            },
            {
                "wild-type": True,
                "organism": "Homo sapiens",
                "turnoverNumber": "316.7",
                "ecNumber": "1.11.1.9"
            },
            {
                "wild-type": True,
                "organism": "Homo sapiens",
                "turnoverNumber": "445",
                "ecNumber": "1.11.1.9"
            },
            {
                "wild-type": True,
                "organism": "Homo sapiens",
                "turnoverNumber": "560",
                "ecNumber": "1.11.1.9"
            },
            {
                "wild-type": True,
                "organism": "Homo sapiens",
                "turnoverNumber": "670",
                "ecNumber": "1.11.1.9"
            }
        ]
    },
    "FA120ACPHi": {},
    "RE1845C": {
        "more": [
            {
                "organism": "Bos taurus",
                "turnoverNumber": "-999",
                "ecNumber": "2.3.1.65"
            }
        ]
    },
    "TRDRm": {
        "GSSG": [],
        "methaneseleninic acid": [],
        "NADH": [
            {
                "organism": "Entamoeba histolytica",
                "turnoverNumber": "0.2",
                "ecNumber": "1.8.1.9"
            },
            {
                "organism": "Methanosarcina acetivorans",
                "turnoverNumber": "0.817",
                "ecNumber": "1.8.1.9"
            }
        ],
        "alloxan": [],
        "Hordeum vulgare thioredoxin disulfide h2": [
            {
                "organism": "Hordeum vulgare",
                "turnoverNumber": "0.8",
                "ecNumber": "1.8.1.9"
            },
            {
                "organism": "Hordeum vulgare",
                "turnoverNumber": "1.31",
                "ecNumber": "1.8.1.9"
            },
            {
                "organism": "Hordeum vulgare",
                "turnoverNumber": "2.98",
                "ecNumber": "1.8.1.9"
            }
        ],
        "protein disulfide-isomerase": [],
        "Hordeum vulgare thioredoxin disulfide h1": [
            {
                "organism": "Hordeum vulgare",
                "turnoverNumber": "3.26",
                "ecNumber": "1.8.1.9"
            }
        ],
        "hydrogen peroxide": [
            {
                "organism": "Mus musculus",
                "turnoverNumber": "6.08",
                "ecNumber": "1.8.1.9"
            }
        ],
        "thioredoxin-CAC": [],
        "thioredoxin P34S": [],
        "thioredoxin disulfide 41": [],
        "thioredoxin": [
            {
                "wild-type": False,
                "organism": "Mus musculus",
                "turnoverNumber": "0.02",
                "ecNumber": "1.8.1.9"
            },
            {
                "wild-type": True,
                "organism": "Caenorhabditis elegans",
                "turnoverNumber": "0.052",
                "ecNumber": "1.8.1.9"
            },
            {
                "wild-type": False,
                "organism": "Rattus norvegicus",
                "turnoverNumber": "0.243",
                "ecNumber": "1.8.1.9"
            },
            {
                "wild-type": False,
                "organism": "Solanum lycopersicum",
                "turnoverNumber": "0.38",
                "ecNumber": "1.8.1.9"
            },
            {
                "wild-type": False,
                "organism": "Entamoeba histolytica",
                "turnoverNumber": "0.5",
                "ecNumber": "1.8.1.9"
            },
            {
                "wild-type": False,
                "organism": "Entamoeba histolytica",
                "turnoverNumber": "1.25",
                "ecNumber": "1.8.1.9"
            },
            {
                "wild-type": False,
                "organism": "Sulfolobus solfataricus",
                "turnoverNumber": "1.3",
                "ecNumber": "1.8.1.9"
            },
            {
                "wild-type": False,
                "organism": "Escherichia coli",
                "turnoverNumber": "2.38",
                "ecNumber": "1.8.1.9"
            },
            {
                "wild-type": False,
                "organism": "Mus musculus",
                "turnoverNumber": "3.5",
                "ecNumber": "1.8.1.9"
            },
            {
                "wild-type": False,
                "organism": "Caenorhabditis elegans",
                "turnoverNumber": "4.03",
                "ecNumber": "1.8.1.9"
            },
            {
                "wild-type": False,
                "organism": "Caenorhabditis elegans",
                "turnoverNumber": "5.3",
                "ecNumber": "1.8.1.9"
            },
            {
                "wild-type": False,
                "organism": "Homo sapiens",
                "turnoverNumber": "5.58",
                "ecNumber": "1.8.1.9"
            },
            {
                "wild-type": False,
                "organism": "Homo sapiens",
                "turnoverNumber": "8.1",
                "ecNumber": "1.8.1.9"
            },
            {
                "wild-type": True,
                "organism": "Caenorhabditis elegans",
                "turnoverNumber": "10.17",
                "ecNumber": "1.8.1.9"
            },
            {
                "wild-type": False,
                "organism": "Escherichia coli",
                "turnoverNumber": "13.2",
                "ecNumber": "1.8.1.9"
            },
            {
                "wild-type": False,
                "organism": "Anopheles gambiae",
                "turnoverNumber": "14.3",
                "ecNumber": "1.8.1.9"
            },
            {
                "wild-type": False,
                "organism": "Anopheles gambiae",
                "turnoverNumber": "15.4",
                "ecNumber": "1.8.1.9"
            },
            {
                "wild-type": False,
                "organism": "Anopheles gambiae",
                "turnoverNumber": "15.7",
                "ecNumber": "1.8.1.9"
            },
            {
                "wild-type": False,
                "organism": "Mus musculus",
                "turnoverNumber": "19.97",
                "ecNumber": "1.8.1.9"
            },
            {
                "wild-type": True,
                "organism": "Escherichia coli",
                "turnoverNumber": "22",
                "ecNumber": "1.8.1.9"
            },
            {
                "wild-type": True,
                "organism": "Escherichia coli",
                "turnoverNumber": "22.8",
                "ecNumber": "1.8.1.9"
            },
            {
                "wild-type": True,
                "organism": "Mus musculus",
                "turnoverNumber": "25",
                "ecNumber": "1.8.1.9"
            },
            {
                "wild-type": True,
                "organism": "Homo sapiens",
                "turnoverNumber": "25.78",
                "ecNumber": "1.8.1.9"
            },
            {
                "wild-type": True,
                "organism": "Homo sapiens",
                "turnoverNumber": "27.4",
                "ecNumber": "1.8.1.9"
            },
            {
                "wild-type": True,
                "organism": "Mus musculus",
                "turnoverNumber": "29.5",
                "ecNumber": "1.8.1.9"
            },
            {
                "wild-type": True,
                "organism": "Mus musculus",
                "turnoverNumber": "37",
                "ecNumber": "1.8.1.9"
            },
            {
                "wild-type": True,
                "organism": "Mus musculus",
                "turnoverNumber": "37.88",
                "ecNumber": "1.8.1.9"
            },
            {
                "wild-type": True,
                "organism": "Rattus norvegicus",
                "turnoverNumber": "41.7",
                "ecNumber": "1.8.1.9"
            },
            {
                "wild-type": True,
                "organism": "Homo sapiens",
                "turnoverNumber": "46.57",
                "ecNumber": "1.8.1.9"
            },
            {
                "wild-type": True,
                "organism": "Rattus norvegicus",
                "turnoverNumber": "50",
                "ecNumber": "1.8.1.9"
            },
            {
                "wild-type": True,
                "organism": "Aeropyrum pernix",
                "turnoverNumber": "63.2",
                "ecNumber": "1.8.1.9"
            },
            {
                "wild-type": True,
                "organism": "Bos taurus",
                "turnoverNumber": "1030",
                "ecNumber": "1.8.1.9"
            },
            {
                "wild-type": True,
                "organism": "Bos taurus",
                "turnoverNumber": "1200",
                "ecNumber": "1.8.1.9"
            },
            {
                "wild-type": True,
                "organism": "Bos taurus",
                "turnoverNumber": "1300",
                "ecNumber": "1.8.1.9"
            }
        ],
        "thioredoxin disulfide 8": [],
        "thioredoxin-R": [],
        "methylseleninate": [
            {
                "organism": "Mus musculus",
                "turnoverNumber": "14",
                "ecNumber": "1.8.1.9"
            },
            {
                "organism": "Homo sapiens",
                "turnoverNumber": "23",
                "ecNumber": "1.8.1.9"
            },
            {
                "organism": "Plasmodium falciparum",
                "turnoverNumber": "31",
                "ecNumber": "1.8.1.9"
            }
        ],
        "thioredoxin disulfide": [
            {
                "organism": "Mus musculus",
                "turnoverNumber": "0.13",
                "ecNumber": "1.8.1.9"
            },
            {
                "organism": "Methanosarcina acetivorans",
                "turnoverNumber": "1.175",
                "ecNumber": "1.8.1.9"
            },
            {
                "organism": "Drosophila melanogaster",
                "turnoverNumber": "4.99",
                "ecNumber": "1.8.1.9"
            },
            {
                "organism": "Drosophila melanogaster",
                "turnoverNumber": "5.8",
                "ecNumber": "1.8.1.9"
            },
            {
                "organism": "Taenia crassiceps",
                "turnoverNumber": "19.2",
                "ecNumber": "1.8.1.9"
            },
            {
                "organism": "Schistosoma mansoni",
                "turnoverNumber": "30",
                "ecNumber": "1.8.1.9"
            },
            {
                "organism": "Mus musculus",
                "turnoverNumber": "37",
                "ecNumber": "1.8.1.9"
            },
            {
                "wild-type": False,
                "organism": "Rattus norvegicus",
                "turnoverNumber": "40.98",
                "ecNumber": "1.8.1.9"
            },
            {
                "wild-type": False,
                "organism": "Rattus norvegicus",
                "turnoverNumber": "44.85",
                "ecNumber": "1.8.1.9"
            },
            {
                "wild-type": False,
                "organism": "Rattus norvegicus",
                "turnoverNumber": "47",
                "ecNumber": "1.8.1.9"
            },
            {
                "wild-type": True,
                "organism": "Rattus norvegicus",
                "turnoverNumber": "94.17",
                "ecNumber": "1.8.1.9"
            },
            {
                "wild-type": True,
                "organism": "Rattus norvegicus",
                "turnoverNumber": "114.1",
                "ecNumber": "1.8.1.9"
            }
        ],
        "FAD": [],
        "more": [
            {
                "organism": "Escherichia coli",
                "turnoverNumber": "-999",
                "ecNumber": "1.8.1.9"
            },
            {
                "wild-type": True,
                "organism": "Escherichia coli",
                "turnoverNumber": "-999",
                "ecNumber": "1.8.1.9"
            },
            {
                "wild-type": True,
                "organism": "Escherichia coli",
                "turnoverNumber": "-999",
                "ecNumber": "1.8.1.9"
            },
            {
                "wild-type": True,
                "organism": "Drosophila melanogaster",
                "turnoverNumber": "-999",
                "ecNumber": "1.8.1.9"
            },
            {
                "wild-type": True,
                "organism": "Mus musculus",
                "turnoverNumber": "-999",
                "ecNumber": "1.8.1.9"
            }
        ],
        "Lipoamide": [
            {
                "wild-type": True,
                "organism": "Rattus norvegicus",
                "turnoverNumber": "2",
                "ecNumber": "1.8.1.9"
            },
            {
                "wild-type": True,
                "organism": "Homo sapiens",
                "turnoverNumber": "3.3",
                "ecNumber": "1.8.1.9"
            },
            {
                "wild-type": False,
                "organism": "Rattus norvegicus",
                "turnoverNumber": "27.6",
                "ecNumber": "1.8.1.9"
            },
            {
                "wild-type": False,
                "organism": "Rattus norvegicus",
                "turnoverNumber": "31.2",
                "ecNumber": "1.8.1.9"
            }
        ],
        "thioredoxin 41": [
            {
                "organism": "Entamoeba histolytica",
                "turnoverNumber": "2.2",
                "ecNumber": "1.8.1.9"
            }
        ],
        "selenocysteine": [],
        "NADPH": [
            {
                "organism": "Solanum lycopersicum",
                "turnoverNumber": "0.35",
                "ecNumber": "1.8.1.9"
            },
            {
                "organism": "Sulfolobus solfataricus",
                "turnoverNumber": "0.61",
                "ecNumber": "1.8.1.9"
            },
            {
                "organism": "Methanosarcina acetivorans",
                "turnoverNumber": "0.65",
                "ecNumber": "1.8.1.9"
            },
            {
                "organism": "Saccharomyces cerevisiae",
                "turnoverNumber": "33.3",
                "ecNumber": "1.8.1.9"
            }
        ],
        "DTNB": [
            {
                "wild-type": False,
                "organism": "Plasmodium falciparum",
                "turnoverNumber": "0.233",
                "ecNumber": "1.8.1.9"
            },
            {
                "wild-type": True,
                "organism": "Plasmodium falciparum",
                "turnoverNumber": "4.58",
                "ecNumber": "1.8.1.9"
            },
            {
                "wild-type": True,
                "organism": "Homo sapiens",
                "turnoverNumber": "29.5",
                "ecNumber": "1.8.1.9"
            },
            {
                "wild-type": False,
                "organism": "Escherichia coli",
                "turnoverNumber": "50.3",
                "ecNumber": "1.8.1.9"
            },
            {
                "wild-type": False,
                "organism": "Rattus norvegicus",
                "turnoverNumber": "66.7",
                "ecNumber": "1.8.1.9"
            }
        ],
        "lipoic acid": [],
        "thioredoxin 1": [],
        "thioredoxin 2": [
            {
                "wild-type": False,
                "organism": "Saccharomyces cerevisiae",
                "turnoverNumber": "47.1",
                "ecNumber": "1.8.1.9"
            }
        ],
        "thioredoxin 3": [],
        "glutaredoxin 4": [],
        "thioredoxin 8": [
            {
                "organism": "Entamoeba histolytica",
                "turnoverNumber": "2.7",
                "ecNumber": "1.8.1.9"
            }
        ],
        "rat thioredoxin": [],
        "5,5'-dithiobis(2-nitrobenzoic acid)": [
            {
                "wild-type": False,
                "organism": "Homo sapiens",
                "turnoverNumber": "0.018",
                "ecNumber": "1.8.1.9"
            },
            {
                "wild-type": False,
                "organism": "Homo sapiens",
                "turnoverNumber": "0.075",
                "ecNumber": "1.8.1.9"
            },
            {
                "wild-type": False,
                "organism": "Entamoeba histolytica",
                "turnoverNumber": "0.23",
                "ecNumber": "1.8.1.9"
            },
            {
                "wild-type": False,
                "organism": "Entamoeba histolytica",
                "turnoverNumber": "0.25",
                "ecNumber": "1.8.1.9"
            },
            {
                "wild-type": False,
                "organism": "Homo sapiens",
                "turnoverNumber": "0.52",
                "ecNumber": "1.8.1.9"
            },
            {
                "wild-type": False,
                "organism": "Homo sapiens",
                "turnoverNumber": "0.55",
                "ecNumber": "1.8.1.9"
            },
            {
                "wild-type": False,
                "organism": "Medicago truncatula",
                "turnoverNumber": "0.62",
                "ecNumber": "1.8.1.9"
            },
            {
                "wild-type": False,
                "organism": "Schistosoma mansoni",
                "turnoverNumber": "1.2",
                "ecNumber": "1.8.1.9"
            },
            {
                "wild-type": False,
                "organism": "Drosophila melanogaster",
                "turnoverNumber": "1.6",
                "ecNumber": "1.8.1.9"
            },
            {
                "wild-type": False,
                "organism": "Solanum lycopersicum",
                "turnoverNumber": "1.77",
                "ecNumber": "1.8.1.9"
            },
            {
                "wild-type": True,
                "organism": "Caenorhabditis elegans",
                "turnoverNumber": "2.23",
                "ecNumber": "1.8.1.9"
            },
            {
                "wild-type": True,
                "organism": "Caenorhabditis elegans",
                "turnoverNumber": "2.23",
                "ecNumber": "1.8.1.9"
            },
            {
                "wild-type": True,
                "organism": "Drosophila melanogaster",
                "turnoverNumber": "2.4",
                "ecNumber": "1.8.1.9"
            },
            {
                "wild-type": True,
                "organism": "Caenorhabditis elegans",
                "turnoverNumber": "2.53",
                "ecNumber": "1.8.1.9"
            },
            {
                "wild-type": True,
                "organism": "Drosophila melanogaster",
                "turnoverNumber": "2.62",
                "ecNumber": "1.8.1.9"
            },
            {
                "wild-type": True,
                "organism": "Drosophila melanogaster",
                "turnoverNumber": "2.62",
                "ecNumber": "1.8.1.9"
            },
            {
                "wild-type": True,
                "organism": "Anopheles gambiae",
                "turnoverNumber": "5.5",
                "ecNumber": "1.8.1.9"
            },
            {
                "wild-type": True,
                "organism": "Homo sapiens",
                "turnoverNumber": "8.28",
                "ecNumber": "1.8.1.9"
            },
            {
                "wild-type": True,
                "organism": "Aeropyrum pernix",
                "turnoverNumber": "9",
                "ecNumber": "1.8.1.9"
            },
            {
                "wild-type": False,
                "organism": "Mus musculus",
                "turnoverNumber": "15.6",
                "ecNumber": "1.8.1.9"
            },
            {
                "wild-type": False,
                "organism": "Schistosoma mansoni",
                "turnoverNumber": "16",
                "ecNumber": "1.8.1.9"
            },
            {
                "wild-type": False,
                "organism": "Homo sapiens",
                "turnoverNumber": "18.73",
                "ecNumber": "1.8.1.9"
            },
            {
                "wild-type": False,
                "organism": "Mus musculus",
                "turnoverNumber": "20.83",
                "ecNumber": "1.8.1.9"
            },
            {
                "wild-type": True,
                "organism": "Mus musculus",
                "turnoverNumber": "20.85",
                "ecNumber": "1.8.1.9"
            },
            {
                "wild-type": True,
                "organism": "Drosophila melanogaster",
                "turnoverNumber": "21.6",
                "ecNumber": "1.8.1.9"
            },
            {
                "wild-type": True,
                "organism": "Homo sapiens",
                "turnoverNumber": "30.02",
                "ecNumber": "1.8.1.9"
            },
            {
                "wild-type": True,
                "organism": "Rattus norvegicus",
                "turnoverNumber": "33.08",
                "ecNumber": "1.8.1.9"
            },
            {
                "wild-type": True,
                "organism": "Homo sapiens",
                "turnoverNumber": "33.33",
                "ecNumber": "1.8.1.9"
            },
            {
                "wild-type": False,
                "organism": "Rattus norvegicus",
                "turnoverNumber": "47.72",
                "ecNumber": "1.8.1.9"
            },
            {
                "wild-type": False,
                "organism": "Mus musculus",
                "turnoverNumber": "48.42",
                "ecNumber": "1.8.1.9"
            },
            {
                "wild-type": False,
                "organism": "Rattus norvegicus",
                "turnoverNumber": "49.87",
                "ecNumber": "1.8.1.9"
            },
            {
                "wild-type": True,
                "organism": "Rattus norvegicus",
                "turnoverNumber": "70.3",
                "ecNumber": "1.8.1.9"
            },
            {
                "wild-type": True,
                "organism": "Rattus norvegicus",
                "turnoverNumber": "106.3",
                "ecNumber": "1.8.1.9"
            }
        ],
        "thioredoxin K36E": [],
        "5-hydroxy-1,4-naphthoquinone": [
            {
                "wild-type": False,
                "organism": "Rattus norvegicus",
                "turnoverNumber": "52.75",
                "ecNumber": "1.8.1.9"
            },
            {
                "wild-type": False,
                "organism": "Rattus norvegicus",
                "turnoverNumber": "174.3",
                "ecNumber": "1.8.1.9"
            }
        ]
    },
    "MCD": {
        "malonyl-CoA": [
            {
                "wild-type": False,
                "organism": "Homo sapiens",
                "turnoverNumber": "13.3",
                "ecNumber": "4.1.1.9"
            },
            {
                "wild-type": False,
                "organism": "Homo sapiens",
                "turnoverNumber": "47.3",
                "ecNumber": "4.1.1.9"
            },
            {
                "wild-type": False,
                "organism": "Homo sapiens",
                "turnoverNumber": "94.6",
                "ecNumber": "4.1.1.9"
            },
            {
                "wild-type": True,
                "organism": "Homo sapiens",
                "turnoverNumber": "109.2",
                "ecNumber": "4.1.1.9"
            },
            {
                "wild-type": False,
                "organism": "Homo sapiens",
                "turnoverNumber": "114.2",
                "ecNumber": "4.1.1.9"
            },
            {
                "wild-type": False,
                "organism": "Homo sapiens",
                "turnoverNumber": "117.1",
                "ecNumber": "4.1.1.9"
            },
            {
                "wild-type": True,
                "organism": "Homo sapiens",
                "turnoverNumber": "128.3",
                "ecNumber": "4.1.1.9"
            },
            {
                "wild-type": True,
                "organism": "Homo sapiens",
                "turnoverNumber": "135",
                "ecNumber": "4.1.1.9"
            },
            {
                "wild-type": False,
                "organism": "Homo sapiens",
                "turnoverNumber": "137.5",
                "ecNumber": "4.1.1.9"
            },
            {
                "wild-type": True,
                "organism": "Homo sapiens",
                "turnoverNumber": "141.2",
                "ecNumber": "4.1.1.9"
            },
            {
                "wild-type": False,
                "organism": "Homo sapiens",
                "turnoverNumber": "162.5",
                "ecNumber": "4.1.1.9"
            },
            {
                "wild-type": False,
                "organism": "Homo sapiens",
                "turnoverNumber": "167.1",
                "ecNumber": "4.1.1.9"
            },
            {
                "wild-type": False,
                "organism": "Homo sapiens",
                "turnoverNumber": "175.4",
                "ecNumber": "4.1.1.9"
            },
            {
                "wild-type": False,
                "organism": "Homo sapiens",
                "turnoverNumber": "208.3",
                "ecNumber": "4.1.1.9"
            }
        ],
        "N-hydroxy-L-ornithine": []
    }}


    initial_brenda_reaction = {
    "GLNLASEer": {
        "N-octanoyl-DL-homoserine lactone": [],
        "5-butyl-4-methyldihydro-2(3H)-furanone": [],
        "gamma-undecanolactone": [
            {
                "wild-type": False,
                "organism": "Sulfolobus solfataricus",
                "turnoverNumber": "3.92",
                "ecNumber": "3.1.1.25"
            },
            {
                "wild-type": False,
                "organism": "Sulfolobus solfataricus",
                "turnoverNumber": "4.25",
                "ecNumber": "3.1.1.25"
            },
            {
                "wild-type": False,
                "organism": "Sulfolobus solfataricus",
                "turnoverNumber": "4.55",
                "ecNumber": "3.1.1.25"
            },
            {
                "wild-type": False,
                "organism": "Sulfolobus solfataricus",
                "turnoverNumber": "4.63",
                "ecNumber": "3.1.1.25"
            },
            {
                "wild-type": True,
                "organism": "Sulfolobus solfataricus",
                "turnoverNumber": "4.95",
                "ecNumber": "3.1.1.25"
            },
            {
                "wild-type": False,
                "organism": "Sulfolobus solfataricus",
                "turnoverNumber": "5.64",
                "ecNumber": "3.1.1.25"
            }
        ],
        "gamma-dodecanolactone": [],
        "N-(3-oxododecanoyl)-L-homoserine lactone": [
            {
                "wild-type": True,
                "organism": "Sulfolobus solfataricus",
                "turnoverNumber": "1.01",
                "ecNumber": "3.1.1.25"
            },
            {
                "wild-type": False,
                "organism": "Sulfolobus solfataricus",
                "turnoverNumber": "1.8",
                "ecNumber": "3.1.1.25"
            },
            {
                "wild-type": False,
                "organism": "Sulfolobus solfataricus",
                "turnoverNumber": "3",
                "ecNumber": "3.1.1.25"
            },
            {
                "wild-type": False,
                "organism": "Sulfolobus solfataricus",
                "turnoverNumber": "6.44",
                "ecNumber": "3.1.1.25"
            }
        ],
        "nonanoic-1,5-lactone": [],
        "gamma-dodecalactone": [],
        "N-(3-oxodecanoyl)-L-homoserine lactone": [
            {
                "wild-type": False,
                "organism": "Sulfolobus solfataricus",
                "turnoverNumber": "0.19",
                "ecNumber": "3.1.1.25"
            },
            {
                "wild-type": False,
                "organism": "Sulfolobus solfataricus",
                "turnoverNumber": "0.6",
                "ecNumber": "3.1.1.25"
            },
            {
                "wild-type": False,
                "organism": "Sulfolobus solfataricus",
                "turnoverNumber": "3.96",
                "ecNumber": "3.1.1.25"
            },
            {
                "wild-type": True,
                "organism": "Sulfolobus solfataricus",
                "turnoverNumber": "4.52",
                "ecNumber": "3.1.1.25"
            }
        ],
        "gamma-dodecanoic lactone": [
            {
                "organism": "Homo sapiens",
                "turnoverNumber": "101",
                "ecNumber": "3.1.1.25"
            }
        ],
        "gamma-heptalactone": [],
        "undecanoic-gamma-lactone": [],
        "N-(2-oxotetrahydrofuran-3-yl)pentanamide": [],
        "N-octanoylhomoserine lactone": [],
        "nonanoic-gamma-lactone": [
            {
                "wild-type": False,
                "organism": "Sulfolobus islandicus",
                "turnoverNumber": "2",
                "ecNumber": "3.1.1.25"
            },
            {
                "wild-type": True,
                "organism": "Sulfolobus islandicus",
                "turnoverNumber": "3.1",
                "ecNumber": "3.1.1.25"
            }
        ],
        "5-(thiobutyl)butyrolactone": [
            {
                "wild-type": False,
                "organism": "Homo sapiens",
                "turnoverNumber": "7.5",
                "ecNumber": "3.1.1.25"
            },
            {
                "wild-type": False,
                "organism": "Homo sapiens",
                "turnoverNumber": "19.4",
                "ecNumber": "3.1.1.25"
            },
            {
                "wild-type": True,
                "organism": "Homo sapiens",
                "turnoverNumber": "116",
                "ecNumber": "3.1.1.25"
            }
        ],
        "N-hexanoylhomoserine lactone": [],
        "N-(3-oxodecanoyl)-DL-homoserine lactone": [],
        "delta-undecalactone": [],
        "delta-dodecalactone": [],
        "gamma-(S)-valerolactone": [],
        "gamma-undecalactone": [],
        "gamma-(R)-valerolactone": [],
        "octanoyl-L-homoserine lactone": [],
        "N-(3-oxododecanoyl)-DL-homoserine lactone": [],
        "gamma-(S)-caprolactone": [],
        "dodecanoic-1,5-lactone": [],
        "gamma-nonanoic acid lactone": [],
        "gamma-heptanolactone": [],
        "Paraoxon": [
            {
                "wild-type": False,
                "organism": "Sulfolobus solfataricus",
                "turnoverNumber": "8.47",
                "ecNumber": "3.1.1.25"
            },
            {
                "wild-type": True,
                "organism": "Sulfolobus solfataricus",
                "turnoverNumber": "12.6",
                "ecNumber": "3.1.1.25"
            }
        ],
        "dodecanoic-gamma-lactone": [],
        "undecanoic-1,5-lactone": [],
        "gamma-heptanolide": [
            {
                "organism": "Sulfolobus acidocaldarius",
                "turnoverNumber": "10.25",
                "ecNumber": "3.1.1.25"
            },
            {
                "organism": "Homo sapiens",
                "turnoverNumber": "34",
                "ecNumber": "3.1.1.25"
            }
        ],
        "delta-undecanolactone": [
            {
                "wild-type": True,
                "organism": "Sulfolobus solfataricus",
                "turnoverNumber": "12.65",
                "ecNumber": "3.1.1.25"
            },
            {
                "wild-type": False,
                "organism": "Sulfolobus solfataricus",
                "turnoverNumber": "44.8",
                "ecNumber": "3.1.1.25"
            },
            {
                "wild-type": False,
                "organism": "Sulfolobus solfataricus",
                "turnoverNumber": "56.8",
                "ecNumber": "3.1.1.25"
            },
            {
                "wild-type": False,
                "organism": "Sulfolobus solfataricus",
                "turnoverNumber": "58",
                "ecNumber": "3.1.1.25"
            },
            {
                "wild-type": False,
                "organism": "Sulfolobus solfataricus",
                "turnoverNumber": "66.5",
                "ecNumber": "3.1.1.25"
            },
            {
                "wild-type": False,
                "organism": "Sulfolobus solfataricus",
                "turnoverNumber": "71.2",
                "ecNumber": "3.1.1.25"
            },
            {
                "wild-type": False,
                "organism": "Sulfolobus solfataricus",
                "turnoverNumber": "93.3",
                "ecNumber": "3.1.1.25"
            }
        ],
        "gamma-nonalactone": [
            {
                "wild-type": True,
                "organism": "Sulfolobus solfataricus",
                "turnoverNumber": "5.54",
                "ecNumber": "3.1.1.25"
            },
            {
                "wild-type": False,
                "organism": "Homo sapiens",
                "turnoverNumber": "5.57",
                "ecNumber": "3.1.1.25"
            },
            {
                "wild-type": True,
                "organism": "Homo sapiens",
                "turnoverNumber": "31",
                "ecNumber": "3.1.1.25"
            },
            {
                "wild-type": True,
                "organism": "Vulcanisaeta moutnovskia",
                "turnoverNumber": "44.49",
                "ecNumber": "3.1.1.25"
            }
        ],
        "N-(3-oxohexanoyl)-L-homoserine lactone": [],
        "N-(3-oxooctanoyl)-L-homoserine lactone": [],
        "3-oxo-octanoyl-L-homoserine lactone": [],
        "gamma-dodecanoic acid lactone": [],
        "gamma-(R)-caprolactone": [],
        "4-methoxy phenyl acetate": [],
        "epsilon-caprolactone": [
            {
                "wild-type": True,
                "organism": "Sulfolobus islandicus",
                "turnoverNumber": "7.27",
                "ecNumber": "3.1.1.25"
            },
            {
                "wild-type": True,
                "organism": "Sulfolobus acidocaldarius",
                "turnoverNumber": "15.04",
                "ecNumber": "3.1.1.25"
            }
        ],
        "Gamma-caprolactone": [
            {
                "wild-type": False,
                "organism": "Homo sapiens",
                "turnoverNumber": "25",
                "ecNumber": "3.1.1.25"
            },
            {
                "wild-type": False,
                "organism": "Homo sapiens",
                "turnoverNumber": "44",
                "ecNumber": "3.1.1.25"
            },
            {
                "wild-type": True,
                "organism": "Homo sapiens",
                "turnoverNumber": "44",
                "ecNumber": "3.1.1.25"
            },
            {
                "wild-type": True,
                "organism": "Vulcanisaeta moutnovskia",
                "turnoverNumber": "112.3",
                "ecNumber": "3.1.1.25"
            }
        ],
        "gamma-butyrolactone": [
            {
                "wild-type": True,
                "organism": "Sulfolobus islandicus",
                "turnoverNumber": "5.75",
                "ecNumber": "3.1.1.25"
            },
            {
                "wild-type": True,
                "organism": "Homo sapiens",
                "turnoverNumber": "111",
                "ecNumber": "3.1.1.25"
            },
            {
                "wild-type": True,
                "organism": "Homo sapiens",
                "turnoverNumber": "111",
                "ecNumber": "3.1.1.25"
            }
        ],
        "delta-valerolactone": [
            {
                "wild-type": False,
                "organism": "Homo sapiens",
                "turnoverNumber": "0.5",
                "ecNumber": "3.1.1.25"
            },
            {
                "wild-type": False,
                "organism": "Homo sapiens",
                "turnoverNumber": "0.9",
                "ecNumber": "3.1.1.25"
            },
            {
                "wild-type": True,
                "organism": "Homo sapiens",
                "turnoverNumber": "29.8",
                "ecNumber": "3.1.1.25"
            },
            {
                "wild-type": False,
                "organism": "Homo sapiens",
                "turnoverNumber": "40",
                "ecNumber": "3.1.1.25"
            },
            {
                "wild-type": False,
                "organism": "Homo sapiens",
                "turnoverNumber": "69.4",
                "ecNumber": "3.1.1.25"
            },
            {
                "wild-type": False,
                "organism": "Homo sapiens",
                "turnoverNumber": "94",
                "ecNumber": "3.1.1.25"
            },
            {
                "wild-type": False,
                "organism": "Homo sapiens",
                "turnoverNumber": "156",
                "ecNumber": "3.1.1.25"
            },
            {
                "wild-type": False,
                "organism": "Homo sapiens",
                "turnoverNumber": "210",
                "ecNumber": "3.1.1.25"
            },
            {
                "wild-type": True,
                "organism": "Homo sapiens",
                "turnoverNumber": "210",
                "ecNumber": "3.1.1.25"
            },
            {
                "wild-type": True,
                "organism": "Homo sapiens",
                "turnoverNumber": "210",
                "ecNumber": "3.1.1.25"
            },
            {
                "wild-type": True,
                "organism": "Homo sapiens",
                "turnoverNumber": "632",
                "ecNumber": "3.1.1.25"
            }
        ],
        "gamma-undecanoiclactone": [],
        "9-oxo-N-(2-oxotetrahydrofuran-3-yl)undecanamide": [],
        "N-(3-oxooctanoyl)-DL-homoserine lactone": [
            {
                "wild-type": False,
                "organism": "Sulfolobus islandicus",
                "turnoverNumber": "0.92",
                "ecNumber": "3.1.1.25"
            },
            {
                "wild-type": False,
                "organism": "Sulfolobus islandicus",
                "turnoverNumber": "0.97",
                "ecNumber": "3.1.1.25"
            },
            {
                "wild-type": True,
                "organism": "Sulfolobus islandicus",
                "turnoverNumber": "4.1",
                "ecNumber": "3.1.1.25"
            }
        ],
        "N-dodecanoylhomoserine lactone": [],
        "nonanoic-delta-lactone": [],
        "7-oxo-N-(2-oxotetrahydrofuran-3-yl)nonanamide": [],
        "dodecanoic-delta-lactone": [],
        "dihydrocoumarin": [
            {
                "organism": "Homo sapiens",
                "turnoverNumber": "152",
                "ecNumber": "3.1.1.25"
            }
        ],
        "N-dodecanoyl-DL-homoserine lactone": [],
        "dodecanoic-1,4-lactone": [],
        "gamma-undecanoic acid lactone": [],
        "delta-nonalactone": [
            {
                "organism": "Homo sapiens",
                "turnoverNumber": "48",
                "ecNumber": "3.1.1.25"
            },
            {
                "organism": "Vulcanisaeta moutnovskia",
                "turnoverNumber": "88.91",
                "ecNumber": "3.1.1.25"
            }
        ],
        "undecanoic-1,4-lactone": [],
        "pantoyl lactone": [],
        "nonanoic-1,4-lactone": [],
        "N-(3-oxohexanoyl)homoserine lactone": [],
        "undecanoic-delta-lactone": [
            {
                "wild-type": False,
                "organism": "Sulfolobus islandicus",
                "turnoverNumber": "12.9",
                "ecNumber": "3.1.1.25"
            },
            {
                "wild-type": False,
                "organism": "Sulfolobus islandicus",
                "turnoverNumber": "14.1",
                "ecNumber": "3.1.1.25"
            },
            {
                "wild-type": True,
                "organism": "Sulfolobus islandicus",
                "turnoverNumber": "17.65",
                "ecNumber": "3.1.1.25"
            }
        ],
        "3-oxo-decanoyl-L-homoserine lactone": [],
        "N-(3-oxooctanoyl)homoserine lactone": []
    },
    "CYSTS": {
        "L-Ser": [],
        "homocysteine": [
            {
                "wild-type": True,
                "organism": "Homo sapiens",
                "turnoverNumber": "6.2",
                "ecNumber": "4.2.1.22"
            },
            {
                "wild-type": True,
                "organism": "Saccharomyces cerevisiae",
                "turnoverNumber": "7.38",
                "ecNumber": "4.2.1.22"
            },
            {
                "wild-type": False,
                "organism": "Homo sapiens",
                "turnoverNumber": "15.5",
                "ecNumber": "4.2.1.22"
            },
            {
                "wild-type": False,
                "organism": "Homo sapiens",
                "turnoverNumber": "32.1",
                "ecNumber": "4.2.1.22"
            },
            {
                "wild-type": True,
                "organism": "Homo sapiens",
                "turnoverNumber": "34",
                "ecNumber": "4.2.1.22"
            }
        ],
        "L-homocysteine": [
            {
                "wild-type": True,
                "organism": "Homo sapiens",
                "turnoverNumber": "0.031",
                "ecNumber": "4.2.1.22"
            },
            {
                "wild-type": False,
                "organism": "Homo sapiens",
                "turnoverNumber": "0.04",
                "ecNumber": "4.2.1.22"
            },
            {
                "wild-type": False,
                "organism": "Homo sapiens",
                "turnoverNumber": "0.09",
                "ecNumber": "4.2.1.22"
            },
            {
                "wild-type": False,
                "organism": "Saccharomyces cerevisiae",
                "turnoverNumber": "0.85",
                "ecNumber": "4.2.1.22"
            },
            {
                "wild-type": True,
                "organism": "Homo sapiens",
                "turnoverNumber": "3.3",
                "ecNumber": "4.2.1.22"
            },
            {
                "wild-type": True,
                "organism": "Homo sapiens",
                "turnoverNumber": "4.66",
                "ecNumber": "4.2.1.22"
            },
            {
                "wild-type": True,
                "organism": "Homo sapiens",
                "turnoverNumber": "7.93",
                "ecNumber": "4.2.1.22"
            },
            {
                "wild-type": False,
                "organism": "Homo sapiens",
                "turnoverNumber": "9.06",
                "ecNumber": "4.2.1.22"
            },
            {
                "wild-type": True,
                "organism": "Homo sapiens",
                "turnoverNumber": "12.7",
                "ecNumber": "4.2.1.22"
            },
            {
                "wild-type": True,
                "organism": "Saccharomyces cerevisiae",
                "turnoverNumber": "17",
                "ecNumber": "4.2.1.22"
            },
            {
                "wild-type": True,
                "organism": "Saccharomyces cerevisiae",
                "turnoverNumber": "21.5",
                "ecNumber": "4.2.1.22"
            }
        ],
        "L-cystathionine": [
            {
                "wild-type": False,
                "organism": "Saccharomyces cerevisiae",
                "turnoverNumber": "0.083",
                "ecNumber": "4.2.1.22"
            },
            {
                "wild-type": False,
                "organism": "Saccharomyces cerevisiae",
                "turnoverNumber": "0.133",
                "ecNumber": "4.2.1.22"
            },
            {
                "wild-type": False,
                "organism": "Saccharomyces cerevisiae",
                "turnoverNumber": "0.418",
                "ecNumber": "4.2.1.22"
            },
            {
                "wild-type": False,
                "organism": "Saccharomyces cerevisiae",
                "turnoverNumber": "0.56",
                "ecNumber": "4.2.1.22"
            },
            {
                "wild-type": True,
                "organism": "Saccharomyces cerevisiae",
                "turnoverNumber": "0.56",
                "ecNumber": "4.2.1.22"
            },
            {
                "wild-type": True,
                "organism": "Saccharomyces cerevisiae",
                "turnoverNumber": "1.03",
                "ecNumber": "4.2.1.22"
            },
            {
                "wild-type": True,
                "organism": "Saccharomyces cerevisiae",
                "turnoverNumber": "6.08",
                "ecNumber": "4.2.1.22"
            },
            {
                "wild-type": True,
                "organism": "Saccharomyces cerevisiae",
                "turnoverNumber": "6.08",
                "ecNumber": "4.2.1.22"
            }
        ],
        "L-cysteine": [
            {
                "wild-type": False,
                "organism": "Homo sapiens",
                "turnoverNumber": "1.95",
                "ecNumber": "4.2.1.22"
            },
            {
                "wild-type": False,
                "organism": "Homo sapiens",
                "turnoverNumber": "3.13",
                "ecNumber": "4.2.1.22"
            },
            {
                "wild-type": True,
                "organism": "Homo sapiens",
                "turnoverNumber": "3.13",
                "ecNumber": "4.2.1.22"
            },
            {
                "wild-type": True,
                "organism": "Homo sapiens",
                "turnoverNumber": "4.39",
                "ecNumber": "4.2.1.22"
            },
            {
                "wild-type": True,
                "organism": "Homo sapiens",
                "turnoverNumber": "4.39",
                "ecNumber": "4.2.1.22"
            }
        ],
        "L-serine": [
            {
                "wild-type": False,
                "organism": "Saccharomyces cerevisiae",
                "turnoverNumber": "0.082",
                "ecNumber": "4.2.1.22"
            },
            {
                "wild-type": False,
                "organism": "Homo sapiens",
                "turnoverNumber": "0.15",
                "ecNumber": "4.2.1.22"
            },
            {
                "wild-type": False,
                "organism": "Saccharomyces cerevisiae",
                "turnoverNumber": "0.45",
                "ecNumber": "4.2.1.22"
            },
            {
                "wild-type": True,
                "organism": "Homo sapiens",
                "turnoverNumber": "0.52",
                "ecNumber": "4.2.1.22"
            },
            {
                "wild-type": False,
                "organism": "Saccharomyces cerevisiae",
                "turnoverNumber": "0.85",
                "ecNumber": "4.2.1.22"
            },
            {
                "wild-type": False,
                "organism": "Homo sapiens",
                "turnoverNumber": "1.3",
                "ecNumber": "4.2.1.22"
            },
            {
                "wild-type": False,
                "organism": "Homo sapiens",
                "turnoverNumber": "1.67",
                "ecNumber": "4.2.1.22"
            },
            {
                "wild-type": False,
                "organism": "Homo sapiens",
                "turnoverNumber": "2.5",
                "ecNumber": "4.2.1.22"
            },
            {
                "wild-type": False,
                "organism": "Homo sapiens",
                "turnoverNumber": "2.9",
                "ecNumber": "4.2.1.22"
            },
            {
                "wild-type": True,
                "organism": "Homo sapiens",
                "turnoverNumber": "3.67",
                "ecNumber": "4.2.1.22"
            },
            {
                "wild-type": False,
                "organism": "Saccharomyces cerevisiae",
                "turnoverNumber": "5.3",
                "ecNumber": "4.2.1.22"
            },
            {
                "wild-type": True,
                "organism": "Homo sapiens",
                "turnoverNumber": "5.4",
                "ecNumber": "4.2.1.22"
            },
            {
                "wild-type": False,
                "organism": "Homo sapiens",
                "turnoverNumber": "5.9",
                "ecNumber": "4.2.1.22"
            },
            {
                "wild-type": False,
                "organism": "Homo sapiens",
                "turnoverNumber": "7.5",
                "ecNumber": "4.2.1.22"
            },
            {
                "wild-type": False,
                "organism": "Saccharomyces cerevisiae",
                "turnoverNumber": "7.6",
                "ecNumber": "4.2.1.22"
            },
            {
                "wild-type": False,
                "organism": "Saccharomyces cerevisiae",
                "turnoverNumber": "8.2",
                "ecNumber": "4.2.1.22"
            },
            {
                "wild-type": True,
                "organism": "Homo sapiens",
                "turnoverNumber": "10.19",
                "ecNumber": "4.2.1.22"
            },
            {
                "wild-type": True,
                "organism": "Homo sapiens",
                "turnoverNumber": "10.2",
                "ecNumber": "4.2.1.22"
            },
            {
                "wild-type": False,
                "organism": "Saccharomyces cerevisiae",
                "turnoverNumber": "13.2",
                "ecNumber": "4.2.1.22"
            },
            {
                "wild-type": True,
                "organism": "Homo sapiens",
                "turnoverNumber": "14.01",
                "ecNumber": "4.2.1.22"
            },
            {
                "wild-type": False,
                "organism": "Homo sapiens",
                "turnoverNumber": "14.6",
                "ecNumber": "4.2.1.22"
            },
            {
                "wild-type": False,
                "organism": "Saccharomyces cerevisiae",
                "turnoverNumber": "14.7",
                "ecNumber": "4.2.1.22"
            },
            {
                "wild-type": False,
                "organism": "Homo sapiens",
                "turnoverNumber": "15.8",
                "ecNumber": "4.2.1.22"
            },
            {
                "wild-type": False,
                "organism": "Saccharomyces cerevisiae",
                "turnoverNumber": "16.8",
                "ecNumber": "4.2.1.22"
            },
            {
                "wild-type": True,
                "organism": "Saccharomyces cerevisiae",
                "turnoverNumber": "17",
                "ecNumber": "4.2.1.22"
            },
            {
                "wild-type": False,
                "organism": "Homo sapiens",
                "turnoverNumber": "19",
                "ecNumber": "4.2.1.22"
            },
            {
                "wild-type": True,
                "organism": "Homo sapiens",
                "turnoverNumber": "19.7",
                "ecNumber": "4.2.1.22"
            },
            {
                "wild-type": True,
                "organism": "Homo sapiens",
                "turnoverNumber": "21",
                "ecNumber": "4.2.1.22"
            },
            {
                "wild-type": True,
                "organism": "Saccharomyces cerevisiae",
                "turnoverNumber": "21.5",
                "ecNumber": "4.2.1.22"
            },
            {
                "wild-type": True,
                "organism": "Homo sapiens",
                "turnoverNumber": "39",
                "ecNumber": "4.2.1.22"
            },
            {
                "wild-type": False,
                "organism": "Homo sapiens",
                "turnoverNumber": "45",
                "ecNumber": "4.2.1.22"
            }
        ],
        "more": []
    },
    "BTNDe": {},
    "BTNDm": {},
    "GTHPe": {
        "cumene peroxide": [
            {
                "organism": "Lucilia cuprina",
                "turnoverNumber": "35.78",
                "ecNumber": "1.11.1.9"
            }
        ],
        "GSH": [
            {
                "organism": "Homo sapiens",
                "turnoverNumber": "19.5",
                "ecNumber": "1.11.1.9"
            },
            {
                "organism": "Homo sapiens",
                "turnoverNumber": "24.5",
                "ecNumber": "1.11.1.9"
            },
            {
                "wild-type": True,
                "organism": "Homo sapiens",
                "turnoverNumber": "221.7",
                "ecNumber": "1.11.1.9"
            },
            {
                "wild-type": True,
                "organism": "Homo sapiens",
                "turnoverNumber": "293.3",
                "ecNumber": "1.11.1.9"
            },
            {
                "wild-type": True,
                "organism": "Homo sapiens",
                "turnoverNumber": "361.7",
                "ecNumber": "1.11.1.9"
            },
            {
                "wild-type": True,
                "organism": "Homo sapiens",
                "turnoverNumber": "408.3",
                "ecNumber": "1.11.1.9"
            }
        ],
        "tert-butyl hydroperoxide": [],
        "H2O2": [
            {
                "organism": "Homo sapiens",
                "turnoverNumber": "10.55",
                "ecNumber": "1.11.1.9"
            },
            {
                "organism": "Homo sapiens",
                "turnoverNumber": "16.02",
                "ecNumber": "1.11.1.9"
            },
            {
                "organism": "Homo sapiens",
                "turnoverNumber": "20.83",
                "ecNumber": "1.11.1.9"
            },
            {
                "organism": "Lucilia cuprina",
                "turnoverNumber": "44.03",
                "ecNumber": "1.11.1.9"
            },
            {
                "wild-type": True,
                "organism": "Homo sapiens",
                "turnoverNumber": "316.7",
                "ecNumber": "1.11.1.9"
            },
            {
                "wild-type": True,
                "organism": "Homo sapiens",
                "turnoverNumber": "445",
                "ecNumber": "1.11.1.9"
            },
            {
                "wild-type": True,
                "organism": "Homo sapiens",
                "turnoverNumber": "560",
                "ecNumber": "1.11.1.9"
            },
            {
                "wild-type": True,
                "organism": "Homo sapiens",
                "turnoverNumber": "670",
                "ecNumber": "1.11.1.9"
            }
        ]
    },
    "GTHPm": {
        "cumene peroxide": [
            {
                "organism": "Lucilia cuprina",
                "turnoverNumber": "35.78",
                "ecNumber": "1.11.1.9"
            }
        ],
        "GSH": [
            {
                "organism": "Homo sapiens",
                "turnoverNumber": "19.5",
                "ecNumber": "1.11.1.9"
            },
            {
                "organism": "Homo sapiens",
                "turnoverNumber": "24.5",
                "ecNumber": "1.11.1.9"
            },
            {
                "wild-type": True,
                "organism": "Homo sapiens",
                "turnoverNumber": "221.7",
                "ecNumber": "1.11.1.9"
            },
            {
                "wild-type": True,
                "organism": "Homo sapiens",
                "turnoverNumber": "293.3",
                "ecNumber": "1.11.1.9"
            },
            {
                "wild-type": True,
                "organism": "Homo sapiens",
                "turnoverNumber": "361.7",
                "ecNumber": "1.11.1.9"
            },
            {
                "wild-type": True,
                "organism": "Homo sapiens",
                "turnoverNumber": "408.3",
                "ecNumber": "1.11.1.9"
            }
        ],
        "tert-butyl hydroperoxide": [],
        "H2O2": [
            {
                "organism": "Homo sapiens",
                "turnoverNumber": "10.55",
                "ecNumber": "1.11.1.9"
            },
            {
                "organism": "Homo sapiens",
                "turnoverNumber": "16.02",
                "ecNumber": "1.11.1.9"
            },
            {
                "organism": "Homo sapiens",
                "turnoverNumber": "20.83",
                "ecNumber": "1.11.1.9"
            },
            {
                "organism": "Lucilia cuprina",
                "turnoverNumber": "44.03",
                "ecNumber": "1.11.1.9"
            },
            {
                "wild-type": True,
                "organism": "Homo sapiens",
                "turnoverNumber": "316.7",
                "ecNumber": "1.11.1.9"
            },
            {
                "wild-type": True,
                "organism": "Homo sapiens",
                "turnoverNumber": "445",
                "ecNumber": "1.11.1.9"
            },
            {
                "wild-type": True,
                "organism": "Homo sapiens",
                "turnoverNumber": "560",
                "ecNumber": "1.11.1.9"
            },
            {
                "wild-type": True,
                "organism": "Homo sapiens",
                "turnoverNumber": "670",
                "ecNumber": "1.11.1.9"
            }
        ]
    },
    "FA120ACPHi": {},
    "RE1845C": {
        "more": [
            {
                "organism": "Bos taurus",
                "turnoverNumber": "-999",
                "ecNumber": "2.3.1.65"
            }
        ]
    },
    "TRDRm": {
        "GSSG": [],
        "methaneseleninic acid": [],
        "NADH": [
            {
                "organism": "Entamoeba histolytica",
                "turnoverNumber": "0.2",
                "ecNumber": "1.8.1.9"
            },
            {
                "organism": "Methanosarcina acetivorans",
                "turnoverNumber": "0.817",
                "ecNumber": "1.8.1.9"
            }
        ],
        "alloxan": [],
        "Hordeum vulgare thioredoxin disulfide h2": [
            {
                "organism": "Hordeum vulgare",
                "turnoverNumber": "0.8",
                "ecNumber": "1.8.1.9"
            },
            {
                "organism": "Hordeum vulgare",
                "turnoverNumber": "1.31",
                "ecNumber": "1.8.1.9"
            },
            {
                "organism": "Hordeum vulgare",
                "turnoverNumber": "2.98",
                "ecNumber": "1.8.1.9"
            }
        ],
        "protein disulfide-isomerase": [],
        "Hordeum vulgare thioredoxin disulfide h1": [
            {
                "organism": "Hordeum vulgare",
                "turnoverNumber": "3.26",
                "ecNumber": "1.8.1.9"
            }
        ],
        "hydrogen peroxide": [
            {
                "organism": "Mus musculus",
                "turnoverNumber": "6.08",
                "ecNumber": "1.8.1.9"
            }
        ],
        "thioredoxin-CAC": [],
        "thioredoxin P34S": [],
        "thioredoxin disulfide 41": [],
        "thioredoxin": [
            {
                "wild-type": False,
                "organism": "Mus musculus",
                "turnoverNumber": "0.02",
                "ecNumber": "1.8.1.9"
            },
            {
                "wild-type": True,
                "organism": "Caenorhabditis elegans",
                "turnoverNumber": "0.052",
                "ecNumber": "1.8.1.9"
            },
            {
                "wild-type": False,
                "organism": "Rattus norvegicus",
                "turnoverNumber": "0.243",
                "ecNumber": "1.8.1.9"
            },
            {
                "wild-type": False,
                "organism": "Solanum lycopersicum",
                "turnoverNumber": "0.38",
                "ecNumber": "1.8.1.9"
            },
            {
                "wild-type": False,
                "organism": "Entamoeba histolytica",
                "turnoverNumber": "0.5",
                "ecNumber": "1.8.1.9"
            },
            {
                "wild-type": False,
                "organism": "Entamoeba histolytica",
                "turnoverNumber": "1.25",
                "ecNumber": "1.8.1.9"
            },
            {
                "wild-type": False,
                "organism": "Sulfolobus solfataricus",
                "turnoverNumber": "1.3",
                "ecNumber": "1.8.1.9"
            },
            {
                "wild-type": False,
                "organism": "Escherichia coli",
                "turnoverNumber": "2.38",
                "ecNumber": "1.8.1.9"
            },
            {
                "wild-type": False,
                "organism": "Mus musculus",
                "turnoverNumber": "3.5",
                "ecNumber": "1.8.1.9"
            },
            {
                "wild-type": False,
                "organism": "Caenorhabditis elegans",
                "turnoverNumber": "4.03",
                "ecNumber": "1.8.1.9"
            },
            {
                "wild-type": False,
                "organism": "Caenorhabditis elegans",
                "turnoverNumber": "5.3",
                "ecNumber": "1.8.1.9"
            },
            {
                "wild-type": False,
                "organism": "Homo sapiens",
                "turnoverNumber": "5.58",
                "ecNumber": "1.8.1.9"
            },
            {
                "wild-type": False,
                "organism": "Homo sapiens",
                "turnoverNumber": "8.1",
                "ecNumber": "1.8.1.9"
            },
            {
                "wild-type": True,
                "organism": "Caenorhabditis elegans",
                "turnoverNumber": "10.17",
                "ecNumber": "1.8.1.9"
            },
            {
                "wild-type": False,
                "organism": "Escherichia coli",
                "turnoverNumber": "13.2",
                "ecNumber": "1.8.1.9"
            },
            {
                "wild-type": False,
                "organism": "Anopheles gambiae",
                "turnoverNumber": "14.3",
                "ecNumber": "1.8.1.9"
            },
            {
                "wild-type": False,
                "organism": "Anopheles gambiae",
                "turnoverNumber": "15.4",
                "ecNumber": "1.8.1.9"
            },
            {
                "wild-type": False,
                "organism": "Anopheles gambiae",
                "turnoverNumber": "15.7",
                "ecNumber": "1.8.1.9"
            },
            {
                "wild-type": False,
                "organism": "Mus musculus",
                "turnoverNumber": "19.97",
                "ecNumber": "1.8.1.9"
            },
            {
                "wild-type": True,
                "organism": "Escherichia coli",
                "turnoverNumber": "22",
                "ecNumber": "1.8.1.9"
            },
            {
                "wild-type": True,
                "organism": "Escherichia coli",
                "turnoverNumber": "22.8",
                "ecNumber": "1.8.1.9"
            },
            {
                "wild-type": True,
                "organism": "Mus musculus",
                "turnoverNumber": "25",
                "ecNumber": "1.8.1.9"
            },
            {
                "wild-type": True,
                "organism": "Homo sapiens",
                "turnoverNumber": "25.78",
                "ecNumber": "1.8.1.9"
            },
            {
                "wild-type": True,
                "organism": "Homo sapiens",
                "turnoverNumber": "27.4",
                "ecNumber": "1.8.1.9"
            },
            {
                "wild-type": True,
                "organism": "Mus musculus",
                "turnoverNumber": "29.5",
                "ecNumber": "1.8.1.9"
            },
            {
                "wild-type": True,
                "organism": "Mus musculus",
                "turnoverNumber": "37",
                "ecNumber": "1.8.1.9"
            },
            {
                "wild-type": True,
                "organism": "Mus musculus",
                "turnoverNumber": "37.88",
                "ecNumber": "1.8.1.9"
            },
            {
                "wild-type": True,
                "organism": "Rattus norvegicus",
                "turnoverNumber": "41.7",
                "ecNumber": "1.8.1.9"
            },
            {
                "wild-type": True,
                "organism": "Homo sapiens",
                "turnoverNumber": "46.57",
                "ecNumber": "1.8.1.9"
            },
            {
                "wild-type": True,
                "organism": "Rattus norvegicus",
                "turnoverNumber": "50",
                "ecNumber": "1.8.1.9"
            },
            {
                "wild-type": True,
                "organism": "Aeropyrum pernix",
                "turnoverNumber": "63.2",
                "ecNumber": "1.8.1.9"
            },
            {
                "wild-type": True,
                "organism": "Bos taurus",
                "turnoverNumber": "1030",
                "ecNumber": "1.8.1.9"
            },
            {
                "wild-type": True,
                "organism": "Bos taurus",
                "turnoverNumber": "1200",
                "ecNumber": "1.8.1.9"
            },
            {
                "wild-type": True,
                "organism": "Bos taurus",
                "turnoverNumber": "1300",
                "ecNumber": "1.8.1.9"
            }
        ],
        "thioredoxin disulfide 8": [],
        "thioredoxin-R": [],
        "methylseleninate": [
            {
                "organism": "Mus musculus",
                "turnoverNumber": "14",
                "ecNumber": "1.8.1.9"
            },
            {
                "organism": "Homo sapiens",
                "turnoverNumber": "23",
                "ecNumber": "1.8.1.9"
            },
            {
                "organism": "Plasmodium falciparum",
                "turnoverNumber": "31",
                "ecNumber": "1.8.1.9"
            }
        ],
        "thioredoxin disulfide": [
            {
                "organism": "Mus musculus",
                "turnoverNumber": "0.13",
                "ecNumber": "1.8.1.9"
            },
            {
                "organism": "Methanosarcina acetivorans",
                "turnoverNumber": "1.175",
                "ecNumber": "1.8.1.9"
            },
            {
                "organism": "Drosophila melanogaster",
                "turnoverNumber": "4.99",
                "ecNumber": "1.8.1.9"
            },
            {
                "organism": "Drosophila melanogaster",
                "turnoverNumber": "5.8",
                "ecNumber": "1.8.1.9"
            },
            {
                "organism": "Taenia crassiceps",
                "turnoverNumber": "19.2",
                "ecNumber": "1.8.1.9"
            },
            {
                "organism": "Schistosoma mansoni",
                "turnoverNumber": "30",
                "ecNumber": "1.8.1.9"
            },
            {
                "organism": "Mus musculus",
                "turnoverNumber": "37",
                "ecNumber": "1.8.1.9"
            },
            {
                "wild-type": False,
                "organism": "Rattus norvegicus",
                "turnoverNumber": "40.98",
                "ecNumber": "1.8.1.9"
            },
            {
                "wild-type": False,
                "organism": "Rattus norvegicus",
                "turnoverNumber": "44.85",
                "ecNumber": "1.8.1.9"
            },
            {
                "wild-type": False,
                "organism": "Rattus norvegicus",
                "turnoverNumber": "47",
                "ecNumber": "1.8.1.9"
            },
            {
                "wild-type": True,
                "organism": "Rattus norvegicus",
                "turnoverNumber": "94.17",
                "ecNumber": "1.8.1.9"
            },
            {
                "wild-type": True,
                "organism": "Rattus norvegicus",
                "turnoverNumber": "114.1",
                "ecNumber": "1.8.1.9"
            }
        ],
        "FAD": [],
        "more": [
            {
                "organism": "Escherichia coli",
                "turnoverNumber": "-999",
                "ecNumber": "1.8.1.9"
            },
            {
                "wild-type": True,
                "organism": "Escherichia coli",
                "turnoverNumber": "-999",
                "ecNumber": "1.8.1.9"
            },
            {
                "wild-type": True,
                "organism": "Escherichia coli",
                "turnoverNumber": "-999",
                "ecNumber": "1.8.1.9"
            },
            {
                "wild-type": True,
                "organism": "Drosophila melanogaster",
                "turnoverNumber": "-999",
                "ecNumber": "1.8.1.9"
            },
            {
                "wild-type": True,
                "organism": "Mus musculus",
                "turnoverNumber": "-999",
                "ecNumber": "1.8.1.9"
            }
        ],
        "Lipoamide": [
            {
                "wild-type": True,
                "organism": "Rattus norvegicus",
                "turnoverNumber": "2",
                "ecNumber": "1.8.1.9"
            },
            {
                "wild-type": True,
                "organism": "Homo sapiens",
                "turnoverNumber": "3.3",
                "ecNumber": "1.8.1.9"
            },
            {
                "wild-type": False,
                "organism": "Rattus norvegicus",
                "turnoverNumber": "27.6",
                "ecNumber": "1.8.1.9"
            },
            {
                "wild-type": False,
                "organism": "Rattus norvegicus",
                "turnoverNumber": "31.2",
                "ecNumber": "1.8.1.9"
            }
        ],
        "thioredoxin 41": [
            {
                "organism": "Entamoeba histolytica",
                "turnoverNumber": "2.2",
                "ecNumber": "1.8.1.9"
            }
        ],
        "selenocysteine": [],
        "NADPH": [
            {
                "organism": "Solanum lycopersicum",
                "turnoverNumber": "0.35",
                "ecNumber": "1.8.1.9"
            },
            {
                "organism": "Sulfolobus solfataricus",
                "turnoverNumber": "0.61",
                "ecNumber": "1.8.1.9"
            },
            {
                "organism": "Methanosarcina acetivorans",
                "turnoverNumber": "0.65",
                "ecNumber": "1.8.1.9"
            },
            {
                "organism": "Saccharomyces cerevisiae",
                "turnoverNumber": "33.3",
                "ecNumber": "1.8.1.9"
            }
        ],
        "DTNB": [
            {
                "wild-type": False,
                "organism": "Plasmodium falciparum",
                "turnoverNumber": "0.233",
                "ecNumber": "1.8.1.9"
            },
            {
                "wild-type": True,
                "organism": "Plasmodium falciparum",
                "turnoverNumber": "4.58",
                "ecNumber": "1.8.1.9"
            },
            {
                "wild-type": True,
                "organism": "Homo sapiens",
                "turnoverNumber": "29.5",
                "ecNumber": "1.8.1.9"
            },
            {
                "wild-type": False,
                "organism": "Escherichia coli",
                "turnoverNumber": "50.3",
                "ecNumber": "1.8.1.9"
            },
            {
                "wild-type": False,
                "organism": "Rattus norvegicus",
                "turnoverNumber": "66.7",
                "ecNumber": "1.8.1.9"
            }
        ],
        "lipoic acid": [],
        "thioredoxin 1": [],
        "thioredoxin 2": [
            {
                "wild-type": False,
                "organism": "Saccharomyces cerevisiae",
                "turnoverNumber": "47.1",
                "ecNumber": "1.8.1.9"
            }
        ],
        "thioredoxin 3": [],
        "glutaredoxin 4": [],
        "thioredoxin 8": [
            {
                "organism": "Entamoeba histolytica",
                "turnoverNumber": "2.7",
                "ecNumber": "1.8.1.9"
            }
        ],
        "rat thioredoxin": [],
        "5,5'-dithiobis(2-nitrobenzoic acid)": [
            {
                "wild-type": False,
                "organism": "Homo sapiens",
                "turnoverNumber": "0.018",
                "ecNumber": "1.8.1.9"
            },
            {
                "wild-type": False,
                "organism": "Homo sapiens",
                "turnoverNumber": "0.075",
                "ecNumber": "1.8.1.9"
            },
            {
                "wild-type": False,
                "organism": "Entamoeba histolytica",
                "turnoverNumber": "0.23",
                "ecNumber": "1.8.1.9"
            },
            {
                "wild-type": False,
                "organism": "Entamoeba histolytica",
                "turnoverNumber": "0.25",
                "ecNumber": "1.8.1.9"
            },
            {
                "wild-type": False,
                "organism": "Homo sapiens",
                "turnoverNumber": "0.52",
                "ecNumber": "1.8.1.9"
            },
            {
                "wild-type": False,
                "organism": "Homo sapiens",
                "turnoverNumber": "0.55",
                "ecNumber": "1.8.1.9"
            },
            {
                "wild-type": False,
                "organism": "Medicago truncatula",
                "turnoverNumber": "0.62",
                "ecNumber": "1.8.1.9"
            },
            {
                "wild-type": False,
                "organism": "Schistosoma mansoni",
                "turnoverNumber": "1.2",
                "ecNumber": "1.8.1.9"
            },
            {
                "wild-type": False,
                "organism": "Drosophila melanogaster",
                "turnoverNumber": "1.6",
                "ecNumber": "1.8.1.9"
            },
            {
                "wild-type": False,
                "organism": "Solanum lycopersicum",
                "turnoverNumber": "1.77",
                "ecNumber": "1.8.1.9"
            },
            {
                "wild-type": True,
                "organism": "Caenorhabditis elegans",
                "turnoverNumber": "2.23",
                "ecNumber": "1.8.1.9"
            },
            {
                "wild-type": True,
                "organism": "Caenorhabditis elegans",
                "turnoverNumber": "2.23",
                "ecNumber": "1.8.1.9"
            },
            {
                "wild-type": True,
                "organism": "Drosophila melanogaster",
                "turnoverNumber": "2.4",
                "ecNumber": "1.8.1.9"
            },
            {
                "wild-type": True,
                "organism": "Caenorhabditis elegans",
                "turnoverNumber": "2.53",
                "ecNumber": "1.8.1.9"
            },
            {
                "wild-type": True,
                "organism": "Drosophila melanogaster",
                "turnoverNumber": "2.62",
                "ecNumber": "1.8.1.9"
            },
            {
                "wild-type": True,
                "organism": "Drosophila melanogaster",
                "turnoverNumber": "2.62",
                "ecNumber": "1.8.1.9"
            },
            {
                "wild-type": True,
                "organism": "Anopheles gambiae",
                "turnoverNumber": "5.5",
                "ecNumber": "1.8.1.9"
            },
            {
                "wild-type": True,
                "organism": "Homo sapiens",
                "turnoverNumber": "8.28",
                "ecNumber": "1.8.1.9"
            },
            {
                "wild-type": True,
                "organism": "Aeropyrum pernix",
                "turnoverNumber": "9",
                "ecNumber": "1.8.1.9"
            },
            {
                "wild-type": False,
                "organism": "Mus musculus",
                "turnoverNumber": "15.6",
                "ecNumber": "1.8.1.9"
            },
            {
                "wild-type": False,
                "organism": "Schistosoma mansoni",
                "turnoverNumber": "16",
                "ecNumber": "1.8.1.9"
            },
            {
                "wild-type": False,
                "organism": "Homo sapiens",
                "turnoverNumber": "18.73",
                "ecNumber": "1.8.1.9"
            },
            {
                "wild-type": False,
                "organism": "Mus musculus",
                "turnoverNumber": "20.83",
                "ecNumber": "1.8.1.9"
            },
            {
                "wild-type": True,
                "organism": "Mus musculus",
                "turnoverNumber": "20.85",
                "ecNumber": "1.8.1.9"
            },
            {
                "wild-type": True,
                "organism": "Drosophila melanogaster",
                "turnoverNumber": "21.6",
                "ecNumber": "1.8.1.9"
            },
            {
                "wild-type": True,
                "organism": "Homo sapiens",
                "turnoverNumber": "30.02",
                "ecNumber": "1.8.1.9"
            },
            {
                "wild-type": True,
                "organism": "Rattus norvegicus",
                "turnoverNumber": "33.08",
                "ecNumber": "1.8.1.9"
            },
            {
                "wild-type": True,
                "organism": "Homo sapiens",
                "turnoverNumber": "33.33",
                "ecNumber": "1.8.1.9"
            },
            {
                "wild-type": False,
                "organism": "Rattus norvegicus",
                "turnoverNumber": "47.72",
                "ecNumber": "1.8.1.9"
            },
            {
                "wild-type": False,
                "organism": "Mus musculus",
                "turnoverNumber": "48.42",
                "ecNumber": "1.8.1.9"
            },
            {
                "wild-type": False,
                "organism": "Rattus norvegicus",
                "turnoverNumber": "49.87",
                "ecNumber": "1.8.1.9"
            },
            {
                "wild-type": True,
                "organism": "Rattus norvegicus",
                "turnoverNumber": "70.3",
                "ecNumber": "1.8.1.9"
            },
            {
                "wild-type": True,
                "organism": "Rattus norvegicus",
                "turnoverNumber": "106.3",
                "ecNumber": "1.8.1.9"
            }
        ],
        "thioredoxin K36E": [],
        "5-hydroxy-1,4-naphthoquinone": [
            {
                "wild-type": False,
                "organism": "Rattus norvegicus",
                "turnoverNumber": "52.75",
                "ecNumber": "1.8.1.9"
            },
            {
                "wild-type": False,
                "organism": "Rattus norvegicus",
                "turnoverNumber": "174.3",
                "ecNumber": "1.8.1.9"
            }
        ]
    },
    "MCD": {
        "malonyl-CoA": [
            {
                "wild-type": False,
                "organism": "Homo sapiens",
                "turnoverNumber": "13.3",
                "ecNumber": "4.1.1.9"
            },
            {
                "wild-type": False,
                "organism": "Homo sapiens",
                "turnoverNumber": "47.3",
                "ecNumber": "4.1.1.9"
            },
            {
                "wild-type": False,
                "organism": "Homo sapiens",
                "turnoverNumber": "94.6",
                "ecNumber": "4.1.1.9"
            },
            {
                "wild-type": True,
                "organism": "Homo sapiens",
                "turnoverNumber": "109.2",
                "ecNumber": "4.1.1.9"
            },
            {
                "wild-type": False,
                "organism": "Homo sapiens",
                "turnoverNumber": "114.2",
                "ecNumber": "4.1.1.9"
            },
            {
                "wild-type": False,
                "organism": "Homo sapiens",
                "turnoverNumber": "117.1",
                "ecNumber": "4.1.1.9"
            },
            {
                "wild-type": True,
                "organism": "Homo sapiens",
                "turnoverNumber": "128.3",
                "ecNumber": "4.1.1.9"
            },
            {
                "wild-type": True,
                "organism": "Homo sapiens",
                "turnoverNumber": "135",
                "ecNumber": "4.1.1.9"
            },
            {
                "wild-type": False,
                "organism": "Homo sapiens",
                "turnoverNumber": "137.5",
                "ecNumber": "4.1.1.9"
            },
            {
                "wild-type": True,
                "organism": "Homo sapiens",
                "turnoverNumber": "141.2",
                "ecNumber": "4.1.1.9"
            },
            {
                "wild-type": False,
                "organism": "Homo sapiens",
                "turnoverNumber": "162.5",
                "ecNumber": "4.1.1.9"
            },
            {
                "wild-type": False,
                "organism": "Homo sapiens",
                "turnoverNumber": "167.1",
                "ecNumber": "4.1.1.9"
            },
            {
                "wild-type": False,
                "organism": "Homo sapiens",
                "turnoverNumber": "175.4",
                "ecNumber": "4.1.1.9"
            },
            {
                "wild-type": False,
                "organism": "Homo sapiens",
                "turnoverNumber": "208.3",
                "ecNumber": "4.1.1.9"
            }
        ],
        "N-hydroxy-L-ornithine": []
        }}

    def test_file_correction(self):
        '''correctJson() should be able to reverse KEGG codes and BRENDA names
            that have been reversed.
        '''
        brenda_keggs = DataTreatment.correctJson('Unit Tests/incorrect_json.json')
        with open('Unit Tests/correct_json.json') as infile:
            correct_file = json.load(infile)
        self.assertEqual(brenda_keggs, correct_file)

    def test_load_brenda(self):
        '''BRENDA parameters must be loaded correctly for program to work.
        '''
        treated_brenda_output = DataTreatment.openJson('Unit Tests/sample_brenda_output.json')
        self.assertEqual(treated_brenda_output, SampleData.initial_input)
        
    def test_convert_brenda_to_data_structure(self):
        '''test that brenda is converted to an Enzyme and MetaboliteCandidate - 
            based structure'''
        #Where does this happen in DataTreatment? 

class FileCorrectionBadInput(unittest.TestCase):
    def test_no_code(self):
        '''correctJson() must have a kegg code in a key:value pair'''
        self.assertRaises(DataTreatment.BadDataError,
                          DataTreatment.correctJson, 'Unit Tests/no_code.json')

    def test_no_file(self):
        '''Throw FileNotFoundError if no file'''
        self.assertRaises(FileNotFoundError, DataTreatment.correctJson,
                          'Unit Tests/no_file_here.json')

    def test_incomplete(self):
        '''File must be populated and be proper JSON.'''
        self.assertRaises(json.decoder.JSONDecodeError, DataTreatment.correctJson,
                          'Unit Tests/incomplete.json')

    def test_empty(self):
        '''File must be populated and be proper JSON.'''
        self.assertRaises(json.decoder.JSONDecodeError, DataTreatment.correctJson,
                          'Unit Tests/empty.json')

class IO(unittest.TestCase):
    '''Meant to test write(). openJson is already tested in TestDataPassing.py
    '''
    def test_write(self):
        file_readout = openJson('Unit Tests/sample_brenda_output.json')
        write('Unit Tests/sample_write_output.json', file_readout)
        self.assertEqual(file_readout, openJson('Unit Tests/sample_write_output.json'))
    
    
if __name__ == '__main__':
    unittest.main()