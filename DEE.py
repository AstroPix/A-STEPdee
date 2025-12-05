
import numpy as np
import h5py as h5
import argparse

from ASTEP_RevCal import ASTEP_RevCal
from ASTEP_Add_BG import ASTEP_Add_BG
from ASTEP_Effects import ASTEP_Effects

"""

This script executes an A-STEP detector effects engine.

The usage is 
>> DEE.py <.sim file name> <Calibration & Resolution file name> <optional -h5 flag> 
<optional --ASTEP_BG_filename flag followed by the path to a csv ASTEP background data>
<optional --seed flag followed by a seed number for smearing>

The .sim file name is a Cosima-generated simulation file, where the detectors lie in the xz plane.
(x,y,z)=(0,0,0) is defined as the center of the inner-most layer. Layer 2 is the outer-most layer.

The calibration & resolution file is an h5 file with the following group structure:

Layer0
|
|--> Chip0
|    |
|    |--> Calibration Array
|    |--> Calibration Column Definitions
|    |--> Resolution Array
|    ---> Resolution Column Definitions
|
|--> Chip1
|    |
|    |--> Calibration Array
|    |--> Calibration Column Definitions
|    |--> Resolution Array
|    ---> Resolution Column Definitions
|
|--> Chip2
|    .
|    .
|    .
|
Layer1
|
|--> Chip0
.
.
.

The calibration and resolution arrays' formats are both Row # | Column # | Slope | Intercept
    * These column names are given in the <> Column Definitions datasets
    * Calibration: Energy = slope x ToT_us + intercept
    * Resolution:  sigma(ToT_us) = slope x ToT_us + intercept
    * This file can also include a Threshold array with Row # | Column # | ToT threshold (default 0 us)
    * This file can also include an Offset array with Row # | Column # | ToT row-column (default 0 us. Defined Column - Row)
    
Including the -h5 flag forces the output of the RevCal step to be an h5 file. Otherwise, it writes to a csv.

Including the --ASTEP_BG_filename flag needs an accompanying path to a .csv file with ASTEP data. 
If given, this script combines the post-reverse calibration simulated data with this dataset, truncating
at the earlier of the ends of each dataset.

Including the --seed flag provides a seed for the random number generator in ASTEP_RevCal that is used for smearing

"""

def parseargs():

    """
    This function handles getting arguments when this function is called

    input 1: path to the TKR sim file
    input 2: path to the TKR calibration file

    """

    parser = argparse.ArgumentParser()
    parser.add_argument("sim_filename", help = "Path to *.sim file from Cosima")
    parser.add_argument("TKR_calib_filename", help = "Path to tracker calibration & resolution file")
    parser.add_argument("--ASTEP_BG_filename", default = '', help = "Empirical A-STEP Background File")
    parser.add_argument("-h5", action = 'store_true', help = "Write outputs as h5?")
    parser.add_argument("--seed", default = -1, help = "Seed for random number generation during smearing")
    args = parser.parse_args()
    return args


def cli():

    args = parseargs()
    sim_filename = args.sim_filename
    TKR_calib_filename = args.TKR_calib_filename
    ASTEP_BG_filename = args.ASTEP_BG_filename
    is_h5 = args.h5
    seed = int(args.seed)
    
    ARC = ASTEP_RevCal(sim_filename,TKR_calib_filename,is_h5,seed)
    ARC.process()
    print('A-STEP .sim File Processed')
    
    if ASTEP_BG_filename != '':
        ABG = ASTEP_Add_BG(ARC,ASTEP_BG_filename,is_h5)
        ABG.process()
        print('A-STEP Background Added')
        with_BG = True
    else:
        print('No A-STEP Background Given')
        with_BG = False
    
    if with_BG:
        AE = ASTEP_Effects(sim_filename,ABG.combined_array_sorted,is_h5,with_BG,ARC)
        AE.process()
    else:
        AE = ASTEP_Effects(sim_filename,ARC.out_array,is_h5,with_BG,ARC)
        AE.process()
    print('A-STEP Instrument Effects Added')
    
    

if  __name__ == '__main__': cli()