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