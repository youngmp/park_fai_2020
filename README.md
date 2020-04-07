# park_fai_2020
Code for Park and Fai, paper to be submitted 2020

We require the following programs/libraries:
Python3 + numpy, scipy, matplotlib.

If one wants to reproduce the XPPAUT .dat files, you will need 
XPPAUT version 8 (http://www.math.pitt.edu/~bard/xpp/xpp.html) and follow the mini tutorial here https://docs.google.com/document/d/1sHhuRor1b937iw7GiIRfmJ13KYhhsLeQu84tbvYQrqw/edit?usp=sharing

generate_figures.py creates all figures using scripts and data in this github directory. Run the script directly using Python3 with the dependencies listed above.

My code tends to be verbose. You will see all kinds of outputs and error messages. If the figures generate, there is no problem. Please contact the first author directly if you run into issues. The script generates largely uncolored images. Shaded regions were filled in by hand using inkscape. Those files are included as .pdf files in this repository. If you wish to edit them in inkscape, open them using inkscape and select Poppler/Cairo import or else the fonts will not look right. The script will also create a trash file called junk.png which is needed to extract the scientific notation for axis labels. it may be deleted at any time.
