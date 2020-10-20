"""Infer read orientation from sample data.""" def infer(iputDir):
    """Main function coordinating the execution of all other functions.
    Should be imported/called from main app and return results to it.
    """
    # implement me
    
    #create .tmp
    mkTmpDir(inputDir)
    
    #Do something then change directory back to previous directory to allow .tmp to 
    #be deleted
    
    os.chdir(inputDir)
    print("Changed working directory to %s " % tmpPath)
    
    # delete .tmp
    rmTmpDir()
    
    
def mkTmpDir(inputDir):
  import os
  # define the name of the directory to be created
  tmpPath = inputDir + "/.tmp"
  
  try:
      os.mkdir(tmpPath)
  except OSError:
      print ("Creation of the directory %s failed" % tmpPath)
  else:
      print ("Successfully created the directory %s " % tmpPath)
  
  os.chdir()
  print("Changed working directory to %s " % tmpPath)
      
def rmTmpDir(tmpPath):
  import os
  
  #remove temporary directory
  try:
    os.rmdir(tmpPath) 
except OSError:
    print ("Deletion of the directory %s failed" % tmpPath)
  else:
    print ("Successfully deleted the directory %s" % tmpPath)
