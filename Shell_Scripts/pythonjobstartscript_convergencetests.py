import subprocess
import glob 
import math 
for h in [ 0.4]:
	for dt in [0.05]:

		NX = int(90/h) + int(90/h)%2
		NY=NX
		# how long should we run for? lets run the quickest for an hour and scale from there
		numhours = int(math.ceil((0.6/h)*(0.1/dt)))
		scriptparameters = open("jobstartmyscript.pbs").read()
		scriptparameters = scriptparameters.replace("NUMHOURS",str(numhours))
		f = open("myscript.pbs","w")
		f.write(scriptparameters)
		f.close()
		# make an appropriate parameters file from the template
		jobstartparameters = open("jobstartparameters").read()
		jobstartparameters = jobstartparameters.replace("INSERT_TIMESTEP","INSERT_TIMESTEP="+str(dt))
		jobstartparameters = jobstartparameters.replace("INSERT_GRIDSPACING","INSERT_GRIDSPACING="+str(h))
		jobstartparameters = jobstartparameters.replace("INSERT_NX","INSERT_NX="+str(NX))
		jobstartparameters = jobstartparameters.replace("INSERT_NY","INSERT_NY="+str(NY))
		f = open("parameters","w")
		f.write(jobstartparameters)
		f.close()
		# make a directory
		directoryname="dt"+str(dt)+"h"+str(h)
		subprocess.call(["mkdir",directoryname])
		# copy everything in
		filenames = glob.glob("*")
		subprocess.call(["cp"] + filenames +[directoryname])
		# okay run it
		subprocess.call([directoryname+"/CompilationScript"] )
	#	subprocess.call(["msub"] + [directoryname+"/"+"myscript.pbs"])



