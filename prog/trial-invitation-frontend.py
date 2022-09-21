#!/usr/bin/env python3
    	
"""
Module to run invitation algorithm via form 
"""
    
# vim: tabstop=8 expandtab shiftwidth=4 softtabstop=4

__version__ = '2.0'
__author__ = 'Adam Brentnall'


from tkinter import filedialog
from tkinter import *
import configparser
from invitedata import InviteData
from inviteuptake import InviteUptake
from invitesolve import InviteSolve


PROGNAME="Trial invitation tool, version 2.0"


class Application(Frame):

    def createWidgets(self):
	
        # N bookings 
        Label(self, text="Number slots \n (total - other gps booked)").grid(row=2, column=1)
        self.innbookings = Entry(self, textvariable=self.nbooking)
        self.innbookings.grid(row=2, column=2)
        
        # Invite level (v1.4)
        Label(self, text="Invitation target \n proportion of remaining capacity \n eg. 0.5 = half remaining slots, 1.1 = 110%").grid(row=3, column=1)
        self.ininvlevel = Entry(self, textvariable=self.invlevel)
        self.ininvlevel.grid(row=3, column=2)

        # Uptake
        Label(self, text="Uptake relative to expected \n (1.0=as expected, 1.2=20% higher, etc)").grid(row=4, column=1)
        self.inupt = Entry(self, textvariable=self.exupt)
        self.inupt.grid(row=4, column=2)
 
        self.RUN = Button(self)
        self.RUN["text"] = "RUN"
        self.RUN["command"] =  self.runmodel
        self.RUN.grid(row=7, column=4)

        self.QUIT = Button(self)
        self.QUIT["text"] = "QUIT"
        self.QUIT["command"] =  self.quit
        self.QUIT.grid(row=8, column=4)

    

    def __load_inpar(self):    
        self.config = configparser.ConfigParser()
        self.config.read(self.configpath)
        self.nMAXBOOK = int(self.config['number_slots']['n'])
        self.iround = int(self.config['invitation_round']['r'])
        self.ininvlevel = float(self.config['invitation_target']['btarget'])
        self.upt = float(self.config['uptake']['upt'])


    def __init__(self, master=None):
        self.configpath = 'input/runpars.ini'
        self.__load_inpar()
        self.nbooking=IntVar()
        self.nbooking.set(self.nMAXBOOK)
        self.invlevel = DoubleVar() 
        self.invlevel.set(self.ininvlevel) 
        self.exupt= DoubleVar()
        self.exupt.set(self.upt)

        Frame.__init__(self, master)
        self.grid()
        self.createWidgets()

    def __updatefromform(self):
        self.ininvlevel = self.invlevel.get() 
        self.nMAXBOOK = self.nbooking.get()
        self.upt = self.exupt.get()


    def __update_config(self):
        self.config['number_slots']['n'] = str(self.nMAXBOOK)
        self.config['invitation_target']['btarget'] = str(self.ininvlevel) 

        self.config['uptake']['upt'] = str(self.upt)

        with open(self.configpath, 'w') as configfile:
            self.config.write(configfile)

    def runmodel(self): # Run the model:

        #update pars
        self.__updatefromform()

        # update config file w/ number bookings etc
        self.__update_config()

        # Load the data for this run:
        self.myinput = InviteData(config_path=self.configpath)

        # Update the uptake model:
        self.myuptake = InviteUptake(nGP=self.myinput.nGP, nAGESEX=self.myinput.nAGESEX, upt=self.myinput.upt)

        # Call the solver:
        self.mysolver = InviteSolve(self.myinput, self.myuptake)

        #write output
        if(self.mysolver.solved):
            self.mysolver.nhsd_output3()
            self.mysolver.summary_uptake()

        else:
            print("!!! Not feasible !!! -> check constraints in config.ini")


       

root = Tk()
root.title(PROGNAME)
app = Application(master=root)
app.mainloop()
root.destroy()
