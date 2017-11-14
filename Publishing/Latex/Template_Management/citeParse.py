# Author: Jordan Barnes
# Date: Fall, 2016
#
# Title: citeParse.py
# Description: applies a regular expression that extracts citation keys
# from a latex document and calls the function templateGen.py to build a
# seperate template for each new citation.
#
# Example usage:
#   input (bash):> printf '\Parencite{Barnes2014}\n\textcite{Barnes2016}' >> parseText.txt
#                > python citeParse.py
#
#   output keys.txt, containing just a list of citation keys.


text2parse = "parseText.txt"
outputKeyList = "keys.txt"

reKeyExpression = "(cite{)(.*)(})"

import re
import os
f = open(text2parse, "r")
fo = open(outputKeyList, "w+")

citeKeys = []
for line in f:
    m = re.search(reKeyExpression, line)
    #print line
    if m is not None:
        print(format(m.group(2)))
        os.system("python templateGen.py --citationKey " + m.group(2))
        citeKeys.extend(m.group(2))
        fo.write(m.group(2) + "\n")
f.close()
fo.close()
