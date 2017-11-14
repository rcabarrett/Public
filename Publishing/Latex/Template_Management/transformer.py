# Author: Jordan Barnes
# Date: Spring, 2017
#
# Title: transformer.py
# Description: feeds a directory of .tex files that need their template changed
# through the template transformer named contentTransform.py
#
# Example usage:
#   input (bash):> python transformer.py



import glob
import os

texFiles = glob.glob('*.tex')

for texF in texFiles:
    os.system("python contentTransform.py --citationKey " + texF.replace('.tex',''))
