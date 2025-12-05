# A-STEPdee

## Purpose

The A-STEP DEE transforms simulated interactions into a form resembling the output of the
detectors. These simulations should be made with MEGAlib Cosima using the geometry in 
AstroPix/A-STEPsimulation, specifically a geometry in which the rocket is aligned along
+z, the detector stack is separated along +y, and the origin is the center of the inner-
most detector.

## Usage

The individual pieces of the DEE is called from the DEE.py script as 

	python DEE.py <.sim filename> <.h5 calibration filename> (-h5) (--ASTEP_BG_filename 
	<.csv background filename>) (--seed <seed number>)

### Required Arguments

<.sim filename>:            The path to the Cosima-simulated .sim file

<.h5 calibration filename>: The path to a file containing calibration and resolution 
information, and optionally, threshold and row-column offset information. This file should
have the following group structure:

Layer0  
|  
|--> Chip0
|    |
|    |--> Calibration Array
|    |--> Resolution Array
|    |--> (Threshold Array)
|    |--> (Offset Array)
|
|--> Chip1
|    |
|    |--> Calibration Array
|    |--> Resolution Array
|    |--> (Threshold Array)
|    |--> (Offset Array)
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

The calibration arrays should be formatted as Row # | Column # | Slope | Intercept, where
(ToT [us])  = Slope x (Energy [keV]) + Intercept.

The resolution arrays should be formatted as Row # | Column # | Slope | Intercept, where
(sigma_ToT [us]) = Slope x (ToT [us]) + Intercept

The optional threshold arrays should be formatted as Row # | Column # | ToT Threshold [us]

The optional offset arrays should be formatted as as Row # | Column # | Offset [us], where
the offset is defined as (column ToT [us]) - (row ToT [us]).

### Optional Arguments

-h5:                                            toggles whether the output will be saved 
as a .csv file or a .h5 file.

--ASTEP_BG_filename <.csv background filename>: the path to empirical data from a no-source
run.

--seed <seed number>:                           The seed for the random number generators

## Process

DEE.py calls three more scripts - ASTEP_RevCal, ASTEP_Add_BG, and ASTEP_Effects. ASTEP_RevCal
takes the .sim file, and converts from (x,y,z,E) to (layer, chip, pixel/row/column, ToT). 
This script then populates a (<# of events> x 12) array in the same format as the quad chip
decoder. ASTEP_Add_BG is called only if a background .csv file is provided with the --ASTEP_BG_filename
keyword. This script performs some cleaning of the FPGA timestamps in both the background 
and simulated arrays. Then it combines the two arrays, sorting them by FPGA timestamps 
(accounting for rollover). Finally, ASTEP_Effects includes any other instrumental effects.
In its current form, this is only handling coincident entries in the output array. When 
entries are in the array within the time it takes the FPGA to read-out, then the entries are 
sorted by what would have been the arrival time at the FPGA and the FPGA timestamps are altered
so that they are at least separated by this read-out time.

## Outputs

There will be up to 3 output files in the same directory as the .sim file. The first is 
written immediately after ASTEP_RevCal and is named *.sim.ASTEP<.h5/.csv>. This file 
is source-only without the extra effects in ASTEP_Effects. The second output file is written
after ASTEP_Add_BG and is named *.sim.ASTEP_wBG_<.h5/.csv>. This file is source + background
but without the extra effects in ASTEP_Effects. The final output file is written after 
ASTEP_Effects and is named either *.sim.ASTEP_wEff<.h5/.csv> or *.sim.ASTEP_wBG_wEff<.h5/.csv>,
depending on if the --ASTEP_BG_filename keyword was used. 




