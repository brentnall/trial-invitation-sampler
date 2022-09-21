#!/usr/bin/env python3
# vim: tabstop=8 expandtab shiftwidth=4 softtabstop=4
    
"""Class to keep data on invitation problem

Reads and stores input / output data
"""

__version__ = '2.0'
__author__ = 'Adam Brentnall'

import numpy as np
import configparser

class InviteData:
    
    """Invitation data"""

    def __init__(self, input_path='input/input.csv', nAGESEX = 12, config_path='input/runpars.ini', utility_path='input/utility.csv'):
        self.unit = 1
        self.input_file_path = input_path
        self.input_config_path = config_path
        self.input_utility_path = utility_path
        self.nAGESEX = nAGESEX
        self.load_data()
        self.load_utility()

    def load_data(self): #Load the data for this run

        inputdatagp = np.genfromtxt(
           self.input_file_path, delimiter=',', usecols=(0), skip_header=1,
           dtype='str')

        self.GPID = inputdatagp

        inputdata = np.genfromtxt(
           self.input_file_path, delimiter=',', skip_header=1,
           dtype='f8')

        self.nGP = inputdata.shape[0]
        self.priority = inputdata[:, 1].astype(float)
        self.deprivation = inputdata[:, 2]
        self.bowel = inputdata[:,3]
        self.nage = inputdata[:, range(4, 4 + self.nAGESEX) ].astype(float)
        self.aage = inputdata[:, range(4 + self.nAGESEX, 4 + self.nAGESEX*2)].astype(float)
        self.tage = inputdata[:, range(4 + self.nAGESEX*2, 4 + self.nAGESEX*3)].astype(float)

        self.config = configparser.ConfigParser()
        self.config.read(self.input_config_path)
        self.nMAXBOOK = float(self.config['number_slots']['n'])
        self.inviteround = float(self.config['invitation_round']['r'])
        self.upt = float(self.config['uptake']['upt'])
        self.invitetarget = float(self.config['invitation_target']['btarget'])

    def load_utility(self): #Load utility data corresponding to invitation problem

        utilities = np.genfromtxt(
            self.input_utility_path, delimiter=',', skip_header=1, usecols=(3,4), dtype='f8')
        self.utility = utilities[:, 0].astype(float)
        self.rate =  utilities[:, 1].astype(float)
        self.utilityarr = np.empty([self.nGP, self.nAGESEX ])
        counter = 0
        for idx in range(0, self.nGP):
            for idy in range(0, self.nAGESEX):
                self.utilityarr[idx, idy] = self.utility[counter]
                counter = counter + 1
        

