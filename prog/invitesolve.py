#!/usr/bin/env python3
# vim: tabstop=8 expandtab shiftwidth=4 softtabstop=4
    
"""Class to solve invitation problem

Reads and stores input / output data
"""

__version__ = '2.0'
__author__ = 'Adam Brentnall'

import numpy as np
import pandas as pd
import configparser
from cvxopt import matrix, spmatrix, spdiag
from cvxopt import solvers
from datetime import date


class InviteSolve:
    
    """Invitation problem to solve"""

    def __init__(self, input_invitedata, input_inviteuptake):
        # Invitation data 
        self.invitedata = input_invitedata
        # Projected uptake
        self.inviteuptake = input_inviteuptake
        # Invitation round
        self.inviteround = int(self.invitedata.inviteround)
        # Solver problem parameters
        self.config = configparser.ConfigParser()
        self.config.read('config.ini')
        # Setup and run
        self._setup()
        # Print out summary
        if(self.solved):
            self.summary() 

    def _setup(self):
        # constants
        self.nGP = self.invitedata.nGP
        self.nAGESEX = self.invitedata.nAGESEX
        self.nX = self.nGP * self.nAGESEX

        # input variables for the LP
        # total number already sent gp/age/sex categories for solver
        self.input_a = np.array(self.invitedata.aage).astype('float').reshape(1, self.nX)
        
        # total number already booked gp/age/sex categories for solver
        self.input_t = np.array(self.invitedata.tage).astype('float').reshape(1, self.nX)
        self.input_tp = np.sum(np.array(self.invitedata.tage).astype('float'), axis=0)
        self.input_tpp = np.sum(self.input_t)

        # total number gp/age/sex categories for solver
        self.input_n =  self.invitedata.nage
        self.input_n_vec = self.invitedata.nage.reshape(1, self.nX)
        #print(self.input_n)

        # uptake
        self.input_u = self.inviteuptake.uptake.reshape(1, self.nX)
        #print("update assump", self.input_u)
        
        ## bounds
        # total number uptake expected / targetted: based on number available, and max capacity
        if(self.inviteround>0): #if invite round = 0 then 
            self.bound_b = float(self.config['invbound']['b'].split(',')[self.inviteround-1]) * ( self.invitedata.nMAXBOOK - self.input_tpp ) ## if using presets
        else: # specify level of booking manually
            self.bound_b = float(self.invitedata.invitetarget) * ( self.invitedata.nMAXBOOK - self.input_tpp )

        # max proportion to invite each gp
        self.bound_g_pc = float(self.config['maxgp']['g'])
        # total number each GP
        self.bound_n_gp= np.sum(self.invitedata.nage, axis=1)
        # constraint in LP on max num to GP
        self.bound_g = (self.bound_g_pc * self.bound_n_gp)
        # min proportion expected to accept in age/sex groups
        self.bound_d_pc = np.array(self.config['uptakeagesex']['d'].split(',')).astype('float')
        ## proportion expected invited male
        self.bound_e = float(self.config['propmen']['e'])

        ## utility function weights
        self.utility = self.invitedata.utility
        self.rate = self.invitedata.rate

        # control on expected advanced cancer rate
        self.bound_cstar = float(self.config['advcancer']['cstar'])

        self.soln = self.solve(self.nGP, self.nAGESEX, self.nX, self.input_n, self.input_n_vec, self.input_a, self.input_u, 
                self.bound_b, self.bound_g, self.bound_d_pc, self.bound_e, self.utility, self.input_t, self.input_tp, self.input_tpp, self.rate, self.bound_cstar)


    def solve(self, nGP, nAGESEX, nX, myn, myn_vec, mya, myu, mybound_b, mybound_g, mybound_d_pc, mybound_e, myutility, myt, mytp, mytpp, myrate, mycstar):
        # Set up and solve the LP
        # n*u, calc'd based on u and n, and used below
        mynu = (myn_vec * myu).reshape(1, nX)
        # expected rate advanced cancer
        mycnu = (myrate * mynu).reshape(1, nX)
        # utility depends on n and u! (forgot until v 0.3)
        mycost = matrix( (myutility * mynu).reshape(nX,1))


        ## d depends on b and min proportion in each cat
        mybound_d = np.maximum((mybound_d_pc * (mytpp + mybound_b) - mytp), 0) 

        # primary endpoint rate control
        mybound_C = mycstar * (mytpp + mybound_b) - np.sum(myt * myrate)

        ## Constraints c
        ##1. Define 0<=x<=1
        # Xjk as real number less than 1
        G1a = spmatrix(1.0, range(nX), range(nX))
        h1a = matrix(np.repeat(1,nX)) 
        # Xik as real number more than 0
        G1b = spmatrix(-1.0, range(nX), range(nX))
        h1b = matrix(np.repeat(0,nX))

        #2. No more than 100% patients in group invited
        G2 = spmatrix(1.0, range(nX), range(nX))
        h2 = matrix( (np.repeat(1,nX) - mya).reshape(nX,1) )

        #3. total invitations yields number of bookings
        A3 = matrix(mynu)
        b3 = matrix(mybound_b)

        #4. Proportion to each GP less than constant
        # Gen A matrix with 1's each GP
        mygp = myn[0,:]
        mynot = np.repeat(0, nX)
        G4m= np.concatenate([mygp, mynot])
        for idx in range(0,(nGP-2)):
            mygp = myn[idx+1, :]
            G4m = np.concatenate([G4m, mygp, mynot])

        G4m = np.concatenate([G4m, myn[nGP-1, :]] )
        G4m = G4m.reshape(nGP, nX)
        G4 = matrix(G4m)
        h4 = matrix(mybound_g)

        #5. Expected uptake each xk greater than min bound
        def agesexvector(startidx, nX, nAGESEX):
            myagesex = np.repeat(0, nX)
            myagesex[np.arange(startidx,nX,nAGESEX)] = 1
            return(myagesex)

        G5m = agesexvector(0, nX, nAGESEX)
        for idx in range(1, nAGESEX):
            G5m = np.concatenate([G5m, agesexvector(idx, nX, nAGESEX)])

        G5m = G5m.reshape(nAGESEX,nX)
        mynuG5 = np.repeat(mynu, nAGESEX).reshape(nX, nAGESEX).transpose()
        G5 = matrix(-G5m * mynuG5)
        h5 = matrix(-mybound_d)
        # (Not used in LP, but for printing info - number invites sent):
        mynG5 = np.repeat(myn_vec, nAGESEX).reshape(nX, nAGESEX).transpose()
        G5_n = matrix(G5m * mynG5)

        #6. Proportion men invited is within a determined range
        mymen = np.concatenate([np.repeat(1, nAGESEX/2), np.repeat(0, nAGESEX/2)])
        mymen2 = np.tile(mymen, nGP).reshape(1,nX)
        mywomen = np.concatenate([np.repeat(0, nAGESEX/2), np.repeat(1, nAGESEX/2)])
        mywomen2 = np.tile(mywomen, nGP).reshape(1,nX)
        # Not used in LP, but for printing info - number invites sent
        A6_nmen = matrix(mymen2 * myn_vec)
        A6_nwomen = matrix(mywomen2 * myn_vec)

        #6b. alternative - require balance by age group
        #nAGE = int(nAGESEX/2)
        #A6alt1 = spmatrix(1.0, range(nAGE), range(nAGE))
        #A6alt2 = matrix([[A6alt1], [-A6alt1]])
        #A6alt3 = np.tile(A6alt2, nGP)
        #A6alt = matrix(A6alt3 * mynu)
        ##print("mynu", mynu)
        ##print("A6alt", A6alt.size)
        #b6alt = matrix( (np.repeat(mybound_e,nAGE) - 0.5).reshape(nAGE,1) )
        ##print("b6alt", b6alt)

        #6c. require balance by gp and age group
        mybound_e_male = (1 - mybound_e) / mybound_e
        nAGE = int(nAGESEX/2)
        A6alt1 = spmatrix(1.0, range(nAGE), range(nAGE))
        ## for gp to balance, this is matrix to use.
        A6alt2 = matrix([[A6alt1 * mybound_e_male], [-A6alt1]])
        ## need a block diagonal matrix with these elements (balance each age/gp group)
        A6alt3 = np.zeros((nGP*nAGE, nGP*nAGESEX))
        rr=nAGE
        cc=nAGESEX
        myblock = np.array(A6alt2)
        r, c = 0, 0
        for i in range(nGP):
            A6alt3[r:r + rr, c:c + cc] = myblock
            r += rr
            c += cc

        A6alt = matrix(A6alt3 * mynu)
        b6alt = matrix( (np.repeat(0, nGP*nAGE)).reshape(nGP*nAGE,1) )
    
        ##7 expected yield advanced cancer greater than a bound
        G7 = matrix(-mycnu)
        h7 = matrix(-mybound_C)
       
        ## set up complete LP
        A = matrix([A3, A6alt])
        b = matrix([b3, b6alt])

        G = matrix([G1a, G1b, G2, G4, G5, G7   ])
        h = matrix([h1a, h1b, h2, h4, h5, h7  ])
       
        c = matrix(mycost)
        ## solve
        ## if using gsl solvers
#        self.sol=solvers.lp(c,G,h,A,b, solver='glpk')
        solvers.options['show_progress'] = False
        self.sol=solvers.lp(c,G,h,A,b)
        
        self.solved = (self.sol['status']=="optimal")

        self.__internals = {
            "A": A,
            "b": b,
            "G": G,
            "h": h,
            "c": c,
            "G1a": G1a,
            "G1b": G1b,
            "G2": G2,
            "A3": A3,
            "G4": G4,
            "G5": G5,
            "G5_n": G5_n,
            "G7": G7,
            "A6": A6alt,
            "A6_nmen": A6_nmen,
            "A6_nwomen": A6_nwomen,
            "h1a": h1a,
            "h1b": h1b,
            "h2": h2,
            "b3": b3,
            "h4": h4,
            "h5": h5,
            "b6": b6alt,
            "b7": h7
            }

        if(self.solved):
            self.xfit = self.sol['x']
            self.xfit_matrix = np.array(self.xfit).reshape(nGP, nAGESEX)
            self.ninvite = (np.array(self.xfit).transpose() * myn_vec).reshape(nGP, nAGESEX)



    def summary(self):
        print("==Invitation algorithm summary stats==")
    
        ## number bookings
        nbook = float((self.__internals["A3"] * self.xfit)[0])
        ninvite = float((matrix(self.input_n_vec) * self.xfit)[0])

        ## number endpoint
        rateadv = float(sum(-(self.__internals["G7"] * self.xfit))) / nbook * 1000

        ##number expected men/women age
        temp = -self.__internals["G5"] * self.xfit
        temp_n = self.__internals["G5_n"] * self.xfit

        print("n expted men %.0f" % float(sum(temp[0:6])))
        print("n invited men %.0f" % sum(self.__internals["A6_nmen"] * self.xfit))
        print("n expted women %.0f" % sum(temp[6:12]))
        print("n invited women %.0f" % sum(self.__internals["A6_nwomen"] * self.xfit))
        print("n expted 50-59y %.0f" % float(temp[0]+temp[1]+temp[6]+temp[7]))
        print("n expted 60-69y %.0f" % float(temp[2]+temp[3]+temp[8]+temp[9]))
        print("n expted 70-77y %.0f" % float(temp[4]+temp[5]+temp[10]+temp[11]))
        print("n invite 50-59y %.0f" % float(temp_n[0]+temp_n[1]+temp_n[6]+temp_n[7]))
        print("n invite 60-69y %.0f" % float(temp_n[2]+temp_n[3]+temp_n[8]+temp_n[9]))
        print("n invite 70-77y %.0f" % float(temp_n[4]+temp_n[5]+temp_n[10]+temp_n[11]))
        
        print("Expected rate endpoint, this round: %.1f" % rateadv) 
        print("Expected number invitations: n= %.0f" % ninvite) 
        print("Expected number bookings: n= %.0f" % nbook) 
        expuptake = nbook/ninvite *100
        print("Expected uptake:  %.1f" % expuptake) 
        


    def write_csv(self, output_file='output/invite.csv', round_digits=0):
        df = pd.DataFrame(np.round(self.ninvite, round_digits).astype('int'), columns=['m_age_1','m_age_2','m_age_3','m_age_4','m_age_5', 'm_age_6', 'f_age_1','f_age_2','f_age_3','f_age_4','f_age_5', 'f_age_6'])
        df['gp'] = self.invitedata.GPID 
        df.to_csv(output_file, index=False)

    def nhsd_output(self, output_file='output/nhsd-invite.csv'):

        df = pd.DataFrame(self.ninvite, columns=['m_age_1','m_age_2','m_age_3','m_age_4','m_age_5', 'm_age_6', 'f_age_1','f_age_2','f_age_3','f_age_4','f_age_5', 'f_age_6'])
        
        df['gp'] = self.invitedata.GPID 

        df2 = pd.melt(df, id_vars='gp')

        df2['sex'] = df2.variable.str.slice(0,1)
        df2['age'] = df2.variable.str.slice(6,7)
        lookmeup={"1": 50,
                "2":55,
                "3":60,
                "4":65,
                "5":70,
                "6":75}
        lookmeup2={"1": 54,
                "2":59,
                "3":64,
                "4":69,
                "5":74,
                "6":77}

        thisage=np.repeat(-99, df2.shape[0])
        thisage2=np.repeat(-99, df2.shape[0])

        for idx in range(0, df2.shape[0]):
            thisage[idx] = lookmeup[df2.age[idx]]
            thisage2[idx] = lookmeup2[df2.age[idx]]

        df2['lower_age_limit'] = thisage
        df2['upper_age_limit'] = thisage2

        df3 = df2.groupby(['gp','age']).sum().loc[:,['value']]

        df4 = df2.loc[df2['sex']=='f'].loc[:,['gp', 'age', 'lower_age_limit', 'upper_age_limit', 'value']]
        df4['value_female'] = np.round(np.maximum(0,df4['value']),0).astype('int')
        
        df4b = df4.loc[:,['gp', 'age', 'lower_age_limit', 'upper_age_limit', 'value_female']]

        df5 = df4b.merge(df3, on=['gp', 'age'])
        df5['number_invitations_requested'] = np.round(np.maximum(0,df5['value']),0).astype('int')

        newf = 100 * df5['value_female'] / df5['number_invitations_requested']
        newf.loc[df5['number_invitations_requested']==0] = '' #np.nan
        df5['percent_female'] = newf
        df5['participant_postcode'] = ''

        df6=df5.loc[:,['gp', 'lower_age_limit', 'upper_age_limit', 'percent_female','participant_postcode', 'number_invitations_requested']].sort_values(by=['gp', 'lower_age_limit', 'upper_age_limit'])

        df6.to_csv(output_file, index=False)


    ## one row male, one row female, % female=0 or 100, NHSD format with today's date
    def nhsd_output2(self, output_file='output/nhsd-invite.csv'):

        df = pd.DataFrame(self.ninvite, columns=['m_age_1','m_age_2','m_age_3','m_age_4','m_age_5', 'm_age_6', 'f_age_1','f_age_2','f_age_3','f_age_4','f_age_5', 'f_age_6'])
        
        df['GP_CODE'] = self.invitedata.GPID 

        df2 = pd.melt(df, id_vars='GP_CODE')

        df2['sex'] = df2.variable.str.slice(0,1)
        df2['age'] = df2.variable.str.slice(6,7)
        lookmeup={"1": 50,
                "2":55,
                "3":60,
                "4":65,
                "5":70,
                "6":75}
        lookmeup2={"1": 54,
                "2":59,
                "3":64,
                "4":69,
                "5":74,
                "6":77}
        lookmeup3={"m": 0,
                "f":100}


        thisage=np.repeat(-99, df2.shape[0])
        thisage2=np.repeat(-99, df2.shape[0])
        thisfemale=np.repeat(-99, df2.shape[0])

        for idx in range(0, df2.shape[0]):
            thisage[idx] = lookmeup[df2.age[idx]]
            thisage2[idx] = lookmeup2[df2.age[idx]]
            thisfemale[idx] = lookmeup3[df2.sex[idx]]

        df2['LOWER_AGE_LIMIT'] = thisage
        df2['UPPER_AGE_LIMIT'] = thisage2
        df2['PERCENT_FEMALE'] = thisfemale
        df2['PARTICIPANT_POSTCODE'] = ''
        df2['NUMBER_OF_INVITATIONS_REQUESTED']=np.round(np.maximum(0, df2['value']),0).astype('int')

        df6=df2.sort_values(by=['GP_CODE', 'LOWER_AGE_LIMIT', 'UPPER_AGE_LIMIT'])

        df7 = df6.copy().loc[df6['NUMBER_OF_INVITATIONS_REQUESTED']>0]

        df7['REQUEST_ID'] = range(1, df7.shape[0]+1)

        df7['NUMBER_OF_INVITEES_SELECTED'] = ''
        df7['NUMBER_OF_INVITEES_REMAINING'] = ''

        today = date.today()

        today_format = today.strftime("%d/%m/%Y")

        df7['DATE_OF_REQUEST'] = today_format
 
        df8=df7.loc[:,['REQUEST_ID', 'DATE_OF_REQUEST', 'GP_CODE', 'LOWER_AGE_LIMIT', 'UPPER_AGE_LIMIT', 'PERCENT_FEMALE',
            'PARTICIPANT_POSTCODE', 'NUMBER_OF_INVITATIONS_REQUESTED', 'NUMBER_OF_INVITEES_SELECTED', 'NUMBER_OF_INVITEES_REMAINING']]
       
        df8.to_csv(output_file, index=False)


    ## adjust for multiplication factor due to opt out (need to request more than want - rather than number want due to process constraint)
    ## requires input/upfactor.csv to lookup upscale factors
    def nhsd_output3(self, output_file='output/nhsd-invite-inflated.csv'):
        
        #run and output unadjusted
        self.nhsd_output2()
       
        #now adjust and output
        scalefile = pd.read_csv("input/upfactor.csv")
        scale_long = pd.melt(scalefile, id_vars="ORG_CODE")
        scale_long['sex'] = scale_long['variable'].str.slice(0,1)
        scale_long['age'] = scale_long['variable'].str.slice(1,3)
        scale_long['id'] = scale_long['ORG_CODE'] + '-' + scale_long['sex'] + '-' + scale_long['age']
        lookupsex={0: "M", 100:"F"}
        
        outputfile = pd.read_csv("output/nhsd-invite.csv")
        outputfile['sex'] = ""
        for idx in range(0, outputfile.shape[0]):
            outputfile.loc[idx, 'sex'] = lookupsex[outputfile.PERCENT_FEMALE[idx]]

        outputfile['id'] = outputfile['GP_CODE'] + '-' + outputfile['sex'] + '-' + outputfile['LOWER_AGE_LIMIT'].astype(str)

        combofile = pd.merge(outputfile, scale_long, on='id', how='left')
        newfile = combofile[['REQUEST_ID', 'DATE_OF_REQUEST', 'GP_CODE', 'LOWER_AGE_LIMIT', 'UPPER_AGE_LIMIT', 'PERCENT_FEMALE']].copy()
        newfile['NUMBER_OF_INVITATIONS_REQUESTED'] = np.ceil(combofile['value'] * combofile['NUMBER_OF_INVITATIONS_REQUESTED']).astype('int') 
        newfile['PARTICIPANT_POSTCODE']=''
        newfile['NUMBER_OF_INVITEES_SELECTED'] = ''
        newfile['NUMBER_OF_INVITEES_REMAINING'] = ''

        newfile_out=newfile.loc[:,['REQUEST_ID', 'DATE_OF_REQUEST', 'GP_CODE', 'LOWER_AGE_LIMIT', 'UPPER_AGE_LIMIT', 'PERCENT_FEMALE',
            'PARTICIPANT_POSTCODE', 'NUMBER_OF_INVITATIONS_REQUESTED', 'NUMBER_OF_INVITEES_SELECTED', 'NUMBER_OF_INVITEES_REMAINING']]

        newfile_out.to_csv(output_file, index=False)

        print("Number of NHSD invitations requested n= %.0f" % newfile['NUMBER_OF_INVITATIONS_REQUESTED'].sum())

    # use predicted uptake 
    def summary_uptake(self):
         
        print("==Predicted trial uptake==")
        
        outputfile = pd.read_csv("output/nhsd-invite.csv")

        uptakefile = pd.read_csv("../inputgen/uptakefiles/uptake-gp-age-sex.csv")

        lookupsex={0: "M", 100:"F"}

        outputfile['sex'] = ""
        for idx in range(0, outputfile.shape[0]):
            outputfile.loc[idx, 'sex'] = lookupsex[outputfile.PERCENT_FEMALE[idx]]
            
        outputfile['id'] = outputfile['GP_CODE'] + '-' + outputfile['sex'] + '-' + outputfile['LOWER_AGE_LIMIT'].astype(str)

        uptake_long = uptakefile.iloc[:,[0,4,5,6,7,8,9,10,11,12,13,14,15]].melt(id_vars="GP_CODE")

        uptake_long["sex"] = uptake_long['variable'].str.slice(0,1)

        uptake_long["ageidx"] = uptake_long['variable'].str.slice(1,2)

        uptake_long["age"] = 45 + uptake_long.loc[:,"ageidx"].astype(int) *5

        uptake_long['id'] = uptake_long['GP_CODE'] + '-' + uptake_long['sex'] + '-' + uptake_long['age'].astype(str)

        combofile = pd.merge(outputfile, uptake_long, on='id', how='left')

        nestup = np.sum(combofile.loc[:,"NUMBER_OF_INVITATIONS_REQUESTED"] * combofile.loc[:,"value"]/100)

        nmissing = np.sum(combofile.loc[combofile["value"].isna(),"NUMBER_OF_INVITATIONS_REQUESTED"])

        ntotal = np.sum(combofile.loc[:,"NUMBER_OF_INVITATIONS_REQUESTED"])

        missfactor = ntotal/(ntotal-nmissing)

        print("Expected bookings n= %.0f" % np.round(nestup*missfactor,0))
        print(f"Expected uptake: {np.round( (nestup*missfactor) / ntotal*100,1)}%")

         
