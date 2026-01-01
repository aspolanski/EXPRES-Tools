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
        eprv: whether this is an EPRV observation (True/False). will add a block of calibration observations after target line.
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


    def make_targetlist(self,target_list_name='expres_target_list',cal_block="./eprv_cal_request.csv",target_list_header="./target_list_header.csv"):

        """
        Function to create targetlist in a Keck/Palomar readable form.
        """

        df = self.targets.copy(deep=True) #make copy of working DataFrame

        #format Star name,RA/Dec, vmag, exposure times 
        df.Star = df.Star.str.pad(14,side='right')
        df.RA = df.RA.str.replace(":"," ")
        df.Dec = df.Dec.str.replace(":"," ")
        df.vmag = 'vmag='+df.vmag.astype(str)
        df['exposure_column'] = df.exptime.astype(str)+"/"+df.max_exptime.astype(str)
        df['num_exp'] = df.num_exp.astype(str)+"x"
        #fill nan comments with empty string
        df.comment = df.comment.fillna('')

        #same for calibrations block and header
        cal_block = pd.read_csv(f"{cal_block}")
        cal_block.Star = cal_block.Star.str.pad(14,side='right')
        cal_block.RA = cal_block.RA.str.replace(":"," ")
        cal_block.Dec = cal_block.Dec.str.replace(":"," ")
        cal_block.vmag = 'vmag='+cal_block.vmag.astype(str)
        cal_block['exposure_column'] = cal_block.exptime.astype(str)+"/"+cal_block.max_exptime.astype(str)
        cal_block['num_exp'] = cal_block.num_exp.astype(str)+"x"
        
        target_list_header = pd.read_csv(f"{target_list_header}")
        target_list_header.Star = target_list_header.Star.str.pad(14,side='right')
        target_list_header.RA = target_list_header.RA.str.replace(":"," ")
        target_list_header.Dec = target_list_header.Dec.str.replace(":"," ")
        target_list_header.vmag = 'vmag='+target_list_header.vmag.astype(str)
        target_list_header['exposure_column'] = target_list_header.exptime.astype(str)+"/"+target_list_header.max_exptime.astype(str)
        target_list_header['num_exp'] = target_list_header.num_exp.astype(str)+"x"

        with open(f"{target_list_name}.tbl", "w") as f:
            f.write('\n\n\n')
            np.savetxt(f,target_list_header[['Star','RA','Dec','eqnx','vmag','exposure_column','num_exp','cts','comment']],delimiter="  ",fmt='%s')
            f.write('\n\n\n')
            for i in range(len(df.index)):
                np.savetxt(f, df.iloc[[i]][['Star','RA','Dec','eqnx','vmag','exposure_column','num_exp','cts','comment']],delimiter="  ",fmt='%s')

                if df.iloc[i].eprv:
                    np.savetxt(f, cal_block[['Star','RA','Dec','eqnx','vmag','exposure_column','num_exp','cts','comment']],delimiter="  ",fmt='%s')
                    
                f.write('\n')


