#!/usr/bin/env python3
# vim: tabstop=8 expandtab shiftwidth=4 softtabstop=4
    
"""Class to store / uptake predicted uptake

Assumed uptake from input csv file
"""

__version__ = '2.0'
__author__ = 'Adam Brentnall'

import numpy as np
import pandas as pd

class InviteUptake:
    
    """Invitation projected uptake"""

    def __init__(self, nGP=4, nAGESEX=12, upt=1.0):

        uptake_file = pd.read_csv("input/uptake.csv")
        self.uptake = upt * uptake_file.iloc[:, 1:(nAGESEX+1) ].to_numpy()



