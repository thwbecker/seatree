import os, subprocess

def uname():
	return os.popen("uname -m").readline()[:-1]

def sys_var(name):
	return os.popen("echo $"+name).readline()[:-1]

def testCommand(command, path=None):
	if (path):
		command = path + os.sep + command
	proc = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	
	output = proc.communicate()
	ret = proc.returncode
	if (ret > 0):
		return False
	else:
		return True

def isResponseYes(response):
	return response.lower().startswith('y')

def isResponseNo(response):
	return response.lower().startswith('n')

def svnCheckout(repository, path, svn="svn"):
	while not testCommand("svn help"):
		if svn == "svn":
			print("SVN could not be found in your system path.")
		svn = input("Path to SVN: (ex: /usr/bin/svn)? ")
		try:
			if not os.path.exists(svn):
				print("'" + svn + "' doesn't exist!")
				svn = "svn"
				continue
		except:
			continue
	
	retval = os.system(svn + " checkout " + repository + " " + path)
	
	return retval == 0

if __name__ == "__main__":
	conmanRepo = "http://geodynamics.org/svn/cig/mc/2D/ConMan/trunk"
	conmanPath = "/tmp/ConMan"
	
	svnCheckout(conmanRepo, conmanPath)
def buildWithMake(pythonRoot, target, extraArgs=""):
	"""
	Build one module through the top-level SEATREE Makefile, i.e.
	run 'make <target>' in the directory above pythonRoot. Returns
	True on success. On failure the make output is written to
	<seatree root>/<target>_build_error.log.
	"""
	seatreeRoot = os.path.abspath(os.path.join(pythonRoot, ".."))
	command = "make -C " + seatreeRoot + " " + target
	if extraArgs:
		command += " " + extraArgs
	print("COMMAND: " + command)
	proc = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	output = proc.communicate()
	if proc.returncode != 0:
		errorFile = os.path.join(seatreeRoot, target + "_build_error.log")
		print("Error building " + target + "...output written to " + errorFile)
		print("Possible causes:")
		print("  * make, gcc or gfortran are not installed")
		print("  * a git submodule is not checked out (git submodule update --init)")
		f = open(errorFile, 'w')
		f.write(output[0].decode("utf-8", "replace"))
		f.write(output[1].decode("utf-8", "replace"))
		f.close()
		return False
	print("Success!")
	return True

def relativeConfPath(pythonRoot, path):
	"""
	Return path relative to the python3 directory when it lies inside
	the SEATREE tree (so that the conf files stay valid when the tree
	is moved), the absolute path otherwise.
	"""
	path = os.path.abspath(path)
	rel = os.path.relpath(path, os.path.abspath(pythonRoot))
	if rel.startswith(".." + os.sep + ".."):
		return path
	return rel
