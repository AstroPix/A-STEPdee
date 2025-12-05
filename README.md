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

The arguments are as follows:

<.sim filename>:            The path to the Cosima-simulated .sim file

<.h5 calibration filename>: The path to a file containing calibration and resolution 
information