#!/usr/bin/env python3

import pandas as pd
import numpy as np


class TargetList(object):
    def __init__(self,targets=None):

        """
        Initialize TargetList object.
        targets - A Pandas dataframe of targets with columns:

        Star: Name of star (ideally Simbad resovable)
        RA: RA of target
        Dec: Dec of target
        eqnx: J2000. or otherwise
        vmag: Host V mag
        exptime: expected exporsure time
        max_exptime: maximum exposure time 
        num_exp: number of exposures to take
        cts: exposure meter counts target
        comment: any appropriate comments
        """

        if targets is not None:
            if not isinstance(targets, pd.DataFrame):
                print("Targets must be a Pandas dataframe")
                exit()
            else:
                self.targets = targets
                #self.targets['comment'] = ['']*len(self.targets.index)


    def remove_target(self,star=None):

        """
        Function to remove target from starlist based on Star name.
        """

        if not star:
            print("Provide the name of target to remove.\n")
            exit()

        else:
            self.targets = self.targets[self.targets['Star'] != star].reset_index(drop=True)


    def make_targetlist(self,target_list_name='expres_target_list',save=False):

        """
        Function to create starlist in a Keck/Palomar readable form.
        """

        df = self.targets.copy(deep=True) #make copy of working DataFrame

        #format Star name,RA/Dec, vmag, exposure times
        df.Star = df.Star.astype(str) + "          " 
        df.RA = df.RA.str.replace(":"," ")
        df.Dec = df.Dec.str.replace(":"," ")
        df.vmag = 'vmag='+df.vmag.astype(str)
        df['exposure_column'] = df.exptime.astype(str)+"/"+df.max_exptime.astype(str)
        df['num_exp'] = df.num_exp.astype(str)+"x"
        
        #fill nan comments with empty string
        df.comment = df.comment.fillna('')

        #make new dataframe
        f = pd.DataFrame({'Target':df.Star,'RA':df.RA,'Dec':df.Dec,'Equinox':df.eqnx,'vmag':df.vmag,'exp':df.exposure_column,'num_exp':df.num_exp,'cts':df.cts,'comment':df.comment}) #create new dataframe

        if save == True:
            np.savetxt(f"{target_list_name}.tbl",f,delimiter="  ",fmt='%s',header='\n\n',comments='')
        return(f)
