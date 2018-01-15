import subprocess
import glob 
import os 
for crossingnumber in ["three","four","five","six","seven","eight"]:
	for tablenumber in range(1,30):
		filename = crossingnumber+str(tablenumber)
		stlfilepath ="/home/maths/matgqz/FNInitialisationFiles/stl/"+filename+".stl"

		if os.path.isfile(stlfilepath):
			# make an appropriate parameters file from the template
			jobstartparameters = open("jobstartparameters").read()
			jobstartparameters = jobstartparameters.replace("INSERT_SURFACE_FILENAME","INSERT_SURFACE_FILENAME="+filename)
			f = open("parameters","w")
			f.write(jobstartparameters)
			f.close()
			# make a directory
			directoryname=filename
			subprocess.call(["mkdir",directoryname])
			# copy everything in
			filenames = glob.glob("*")
			subprocess.call(["cp"] + filenames +[directoryname])
			# grab the right stl file
			subprocess.call(["cp"] + [stlfilepath] +[directoryname])
			# okay run it
			childdirectoryname = "/home/maths/matgqz/FNSimulation/Simulation/"+filename+"/"
			subprocess.call([childdirectoryname+"CompilationScript"],cwd=childdirectoryname)
			subprocess.call(["msub", "/home/maths/matgqz/FNSimulation/Simulation/"+filename+"/myscript.pbs"],cwd=childdirectoryname)
		



