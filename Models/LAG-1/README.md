--------------------------------------------------------------------------------
# LAG-1 - Learning, Attention, and Gaze 1.0
A free and open source MATLAB implementation of the integration
of learning, attention and gaze during category learning, as outlined 
in the theory called LAG-1 (Barnes, Blair, Walshe & Tupper, submitted).
--------------------------------------------------------------------------------

Version 1.0

Written by Jordan Barnes, Mark Blair, Calen Walshe and Paul Tupper at
Simon Fraser University and University of Texas, CPS.

GPL license.

## OVERVIEW

LAG-1 is a simple package that can be run from Matlab or compiled and run as a
standalone application. It is so named to highlight the important connections between learning and looking in visual category learning experiments. It is known to exhibit stimulus responsive attention for 1 feature contingency within a category but could likely do more. The code is quite "experimental" so don't be surprised if you need to ask for help with running it out of the box.

You can download the latest version of LAG-1 from the git repository at:

https://github.com/SFU-Cognitive-Science-Lab/LAG-1


## INSTALLATION

Clone or decompress the archive in a folder of your choice and add this folder to to your MATLAB path. 

LAG-1 requires a Matlab installation, available for Windows, Mac OSX and Linux, or the installation of Matlab runtime libraries which are available here:

http://www.mathworks.com/products/compiler/mcr/

LAG-1, as first committed, has been tested on versions of Matlab dating as far back as 2007.


## QUICK START

* Begin by defining a category structure in experimentConfig.txt (and in several random spots as experiment specific conditions until cleaned up). 

* Set the desired values of parameter_source.txt

* Call the model using several different options:


	### Example usage

	**Default options**

	[accuracy, category, fixOrder, subjectNumber, Duration, recordSaccStart, recordSaccEnd, endFit, wtHist, accuracyLevels, Responses, ResponseDist] = TwoDSimulator(varargin)

	**Specific category structure** 
	* 1 = ET3 (code name for Meier2013)  
	* 2 = SSHRC_IF (code name for McColeman2011)  
	* 3 = 5/4 (code name for Rehder2005)  
	* 6 = Kruschke2005 blocking structure  
	* 7 = Simple blocking structure

	[accuracy, category, fixOrder, subjectNumber, Duration, recordSaccStart, recordSaccEnd, endFit, wtHist, accuracyLevels, Responses, ResponseDist] = TwoDSimulator(1)

 	**Specific category structure and learning type** 
 	
	* 1 = corrective Hebbian association
	* 2 = only correct Hebbian association 
	* 3 = unsupervised decoupling of chosen answer with features 
	* 4 = unsupervised decoupling of chosen answer with features and coupling of alternative features with alternative categories.

	[accuracy, category, fixOrder, subjectNumber, Duration, recordSaccStart, recordSaccEnd, endFit, wtHist, accuracyLevels, Responses, ResponseDist] = TwoDSimulator(3,1)


## Agenda

Integration with Cosivina based pupil model with reinforcement learning. (Almost complete. Will take care of a number of obvious shortcomings with the organization of the main functions/scripts).

