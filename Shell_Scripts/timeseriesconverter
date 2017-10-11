import os
from decimal import *
from shutil import copy2

names = os.listdir("./knotplots");
cwd = os.getcwd()+"/knotplots";

Period = Decimal('11.2');

for name in names[:]:
    bits = name.split('_')
    bobs = bits[1].split('.vtk')
    number = Decimal(bobs[0]);

    divisor = number/Period;
    if divisor%1==0:
        oldfilename = cwd + '/' + name
        hackedformatnumber = "%.0f" % divisor
        newfilename = cwd+ '/velocityknotplot' + str(hackedformatnumber) +'.vtk'
        copy2(oldfilename,newfilename);





