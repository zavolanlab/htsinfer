"""Infer read orientation from sample data."""


def infer(inputDir):
    """Main function coordinating the execution of all other functions.
    Should be imported/called from main app and return results to it.
    """
    
    #create .tmp directory in "inputDir" and retuns .tmp path
    tmpDir = mkTmpDir(inputDir)
    
    #Call functions of issues #25, 26, 27 and 34
    
    
    # delete .tmp
    rmTmpDir(tmpDir)
    
    
def mkTmpDir(inputDir):
    """Function to create a temporary directory for the storage of 
    temporary data."""
  
    import os

    # define the name of the directory to be created
    tmpPath = inputDir + "/.tmp"
  
    try:
        os.mkdir(tmpPath)
    except OSError:
      print ("Creation of the directory %s failed" % tmpPath)
    else:
      print ("Successfully created the directory %s " % tmpPath)
  
    return (tmpPath)
      
      
def rmTmpDir(tmpPath):
    """Function to remove the temporary directory."""
  
    import os
    #remove temporary directory
    try:
        os.rmdir(tmpPath)
    except OSError:
        print ("Deletion of the directory %s failed" % tmpPath)
    else:
        print ("Successfully deleted the directory %s" % tmpPath)
