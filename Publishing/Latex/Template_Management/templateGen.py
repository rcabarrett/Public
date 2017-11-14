# Author: Jordan Barnes
# Date: Fall, 2016
#
# Title: templateGen.py
# Description: Accepts a citation key as an input and inserts it into a latex
# template which you will to define below.
#
# Additional credit: http://stackoverflow.com/questions/8085520/generating-pdf-latex-with-python-script
#
# Example usage:
#   input (bash):> python templateGen.py --citationKey "Barnes2016"


import argparse
import os
import subprocess

# Begin latex template. YOU WILL NEED TO MODIFY THIS TO SUIT YOUR NEEDS. What
# you're looking at is just what the initial author was using for his template.

content = r'''\documentclass[jou,american]{apa6}
\usepackage{xstring}
%%  Review specific parameters
\newcommand{\Braincoin}{OFF} %%ON vs. OFF
\newcommand{\mainpaper}{%(citationKey)s}


\IfStrEq{\Braincoin}{ON}{
    \input{../header.tex}
	\newcommand{\pcoin}{%(ppcoin)s}
	\newcommand{\bcoin}{%(bitcoin)s}
}{
    \input{../header_include_style.tex}
}
\newcommand{\concept}{Required.}
\newcommand{\brief}{Required.}

\IfStrEq{\Braincoin}{ON}{
    \input{../init.tex}
}{
    \input{../init_plain.tex}
}

\subsection{Background}

Required.

\subsection{Method}

Required.

\subsection{Results}

Required.

\subsection{Ideas/Citation}

Required.

\IfStrEq{\Braincoin}{ON}{
    \input{../footer.tex}
}{
    \input{../footer_include_style.tex}
}

\end{document}
'''

parser = argparse.ArgumentParser()
parser.add_argument('-c', '--citationKey', default='noKey')
parser.add_argument('-p', '--ppcoin', default='noPeercoin')
parser.add_argument('-b', '--bitcoin', default='noBitcoin')

args = parser.parse_args()

docTitle = args.citationKey + ".tex"
docPdf = args.citationKey + ".pdf"


with open(docTitle,'w') as f:
    f.write(content%args.__dict__)


# If automatic pdf compilation is additionally desired, the commented out commands
# below will get you started.

#cmd = ['pdflatex', '-synctex', '1', '-interaction', 'nonstopmode', docTitle]
#proc = subprocess.Popen(cmd)
#proc.communicate()
#cmd1 = ['biber',args.citationKey]
#proc = subprocess.Popen(cmd1)
#proc.communicate()
#cmd2 = ['pdflatex', '-synctex', '1', '-interaction', 'nonstopmode', docTitle]
#proc = subprocess.Popen(cmd2)
#proc.communicate()
#cmd3 = ['pdflatex', '-synctex', '1', '-interaction', 'nonstopmode', docTitle]
#proc = subprocess.Popen(cmd3)
#proc.communicate()


# retcode = proc.returncode
# if not retcode == 0:
#     os.unlink(docPdf)
#     raise ValueError('Error {} executing command: {}'.format(retcode, ' '.join(cmd)))
#
# #os.unlink(docTitle)
# os.unlink(args.citationKey + ".log")
# os.unlink(args.citationKey + ".aux")
# os.unlink(args.citationKey + ".bbl")
# os.unlink(args.citationKey + ".bcf")
# os.unlink(args.citationKey + ".blg")
# os.unlink(args.citationKey + ".out")
# os.unlink(args.citationKey + ".run.xml")
# os.unlink(args.citationKey + ".synctex.gz")
