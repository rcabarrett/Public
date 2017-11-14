# Python scripts for working with/converting between, latex templates.

Example usage:

"I need to write a bunch of reviews for all the papers I cited in parseText.txt.
If only I could grab all the citations I used and create separate pdf reviews
for each of these papers without having to do a bunch of tedious copy and
pasting...oh citeParse.py can create a list of everything I cited? That gets me
started. The person that wrote that must be super cool."  

> python citeParse.py

"Cool, now I have a file called keys.txt that has all the citation keys I used
in my paper called parseText.txt. If only I had an automated way to insert those
keys into a latex literature review template I like. Wait, citeParse already did this
for me using templateGen.py??"  

> python templateGen.py --citationKey "Barnes2016"  

"Whoa. The guy who wrote this is amazing. But oh no, I don't like the template that
I defined in templateGen.py any more and I've already written a bunch of stuff
into the .tex files it created. Whatever shall I do? OMG transformer.py exists!
Just by putting my new latex template in contentTransform.py I can run:

> python transformer.py  

which will extract everything I've written so far in this directory and
repeatedly transform each file. It uses the filenames as citation keys:

> python contentTransform.py --citationKey "Barnes2016"

and puts the new files into a subdirectory called ./Converted. This is so crazy
awesome that I wont even complain about all the regular expressions and hard
coded file name parameters that I had to deal with throughout this process.
What an awesome person it must have been to have worked through this tedious
problem initially!"
