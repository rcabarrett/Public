# Author: Jordan Barnes
# Date: Spring, 2017
#
# Title: contentTransform.py
# Description: Takes a citation key, e.g. Barnes2016, and converts an existing
# tex file of that name (i.e. Barnes2016.tex) to "Two column article,
# Source: http://www.howtotex.com/"
#
# Note:
#
# Example usage:
#   Precondition: ./Barnes2016.tex exists
#   input (bash):> mkdir Converted
#               :> python contentTransform.py --citationKey "Barnes2016"
#   output: ./Converted/Barnes2016.tex


import argparse
import os
import subprocess
import re

# Change these regular expressions according to what will capture your existing
# text. In this example body text had been written under the headings "Concept",
# "Brief", and "Content".

reConcept = "(ncept\}\{)(.*)([^a]\\\\n)"
reBrief = "(brief}{)(.*)(}\\\\IfStrEq)"
reContent = "(ground\})(.*)(\\\IfStrEq)"

parser = argparse.ArgumentParser()
parser.add_argument('-c', '--citationKey', default='noKey')
args = parser.parse_args()

fname = args.citationKey + ".tex"


with open(fname) as g:
    bodyText=''.join(line.rstrip() for line in g)

mConcept = re.search(reConcept, bodyText).group(2)
mBrief = re.search(reBrief, bodyText).group(2)
mContent = re.search(reContent, bodyText).group(2)



content = r'''
%% LaTeX Template: Two column article, Source: http://www.howtotex.com/
%% Feel free to distribute this template, but please keep to referal to http://www.howtotex.com/ here.
%% Date: February 2011
%% Modified: Jordan B., 2017. Removing comments for python generation.
\documentclass[	DIV=calc,paper=a4,fontsize=11pt,twocolumn]{scrartcl}

\usepackage[english]{babel}
\usepackage[protrusion=true,expansion=true]{microtype}
\usepackage{amsmath,amsfonts,amsthm}
\usepackage[pdftex]{graphicx}
\usepackage[svgnames]{xcolor}
\usepackage[hang, small,labelfont=bf,up,textfont=it,up]{caption}
\usepackage{epstopdf}
\usepackage{subfig}
\usepackage{booktabs}
\usepackage{fix-cm}
\usepackage{sectsty}
\allsectionsfont{\usefont{OT1}{phv}{b}{n}}
\sectionfont{\usefont{OT1}{phv}{b}{n}}
\usepackage{fancyhdr}
\pagestyle{fancy}
\usepackage{lastpage}

\lhead{}
\chead{}
\rhead{}

%%\lfoot{\footnotesize \texttt{HowToTeX.com} \textbullet ~Two column article template}
%%\cfoot{}
%%\rfoot{\footnotesize page \thepage\ of \pageref{LastPage}}

\lfoot{}
\cfoot{}
\rfoot{}
\renewcommand{\headrulewidth}{0.0pt}
\renewcommand{\footrulewidth}{0.4pt}


\usepackage{lettrine}
\newcommand{\initial}[1]{\lettrine[lines=3,lhang=0.3,nindent=0em]{\color{DarkGoldenrod}{\textsf{#1}}}{}}



\usepackage{titling}
\newcommand{\HorRule}{\color{DarkGoldenrod}\rule{\linewidth}{1pt}}

%% Custom addition - JB 27/02/2017
\usepackage[
    backend=biber,
    style=apa,
    sortcites=true,
    sorting=nyt,
    maxcitenames=99,
    mincitenames=4,
    uniquename=false,
    citestyle=authoryear
]{biblatex}
\usepackage{float}
\usepackage[autostyle]{csquotes}
\graphicspath{../images}

\DeclareLanguageMapping{american}{american-apa}
\AtEveryCitekey{\renewcommand*{\finalnamedelim}{\ifthenelse{\value{listcount}>\maxprtauth}
      {}
      {\ifthenelse{\value{liststop}>2}
         {\finalandcomma\addspace\bibstring{and}\space}
         {\addspace\bibstring{and}\space}}}}

\AtBeginBibliography{\renewcommand*{\finalnamedelim}{\ifthenelse{\value{listcount}>\maxprtauth}
      {}
      {\ifthenelse{\value{liststop}>2}
         {\finalandcomma\addspace\bibstring{and}\space}
         {\addspace\bibstring{and}\space}}}}

\addbibresource{../../library.bib}
\usepackage{filecontents}

\newcommand{\tempCiteKey}{%(key5)s}
\newcommand{\mainpaper}{\citetitle{\tempCiteKey}}
\newcommand{\auths}{\cite{\tempCiteKey}}

%%

\pretitle{\vspace{-30pt} \begin{flushleft} \HorRule
				\fontsize{15}{30} \usefont{OT1}{phv}{b}{n} \color{DarkRed} \selectfont
				}

\title{\mainpaper}
\posttitle{\par\end{flushleft}\vskip 0.5em}
\preauthor{\begin{flushleft}
					\large \lineskip 0.5em \usefont{OT1}{phv}{b}{sl} \color{DarkRed}}
\author{\auths}
\postauthor{\footnotesize \usefont{OT1}{phv}{m}{sl} \color{Black}\par\end{flushleft}\HorRule}
\date{}

\setcounter{secnumdepth}{-2}

%% Begin document
\begin{document}
\maketitle
\thispagestyle{fancy}




\initial{%(key1)s}\textbf{%(key2)s}


\section*{Brief}

%(key3)s

\subsection*{Background}

%(key4)s


\end{document}
''' % {'key1': mConcept[0], 'key2': mConcept[1:], 'key3': mBrief, 'key4': mContent, 'key5': args.citationKey}


docTitle = './Converted/' + fname
with open(docTitle,'w') as f:
    f.write(content)
