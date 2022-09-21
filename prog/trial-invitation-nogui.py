#!/usr/bin/env python3
# vim: tabstop=8 expandtab shiftwidth=4 softtabstop=4
    
"""Script to determine number of invitations to send in a trial

Uses a linear program to determine number to send.
"""

__version__ = '2.0'
__author__ = 'Adam Brentnall'

from invitedata import InviteData
from inviteuptake import InviteUptake
from invitesolve import InviteSolve

import numpy as np
import pandas as pd

def runmodel(): # Run the model:
	
    # Load the data for this run:
    myinput = InviteData()   

    # Update the uptake model:
    myuptake = InviteUptake(nGP=myinput.nGP, nAGESEX=myinput.nAGESEX, upt=myinput.upt)

    # Call the solver:
    mysolver = InviteSolve(myinput, myuptake)

    # writeoutput
    if(mysolver.solved):
        mysolver.write_csv(output_file='output/output-proto.csv')
    else:
        print("!!! Not feasible !!! -> check constraints ")

    mysolver.nhsd_output3()
    mysolver.summary_uptake()



if __name__ == "__main__":
    runmodel()
