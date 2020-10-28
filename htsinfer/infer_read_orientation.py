"""Infer read orientation from sample data."""

from typing import Union

def infer(
    file_1: str,  # pylint: disable=unused-argument
    file_2: str = None,  # pylint: disable=unused-argument
    organism: Union[int, str] = "hsapiens",  # pylint: disable=unused-argument
) -> str:
    """Infers read orientation for single- or paired-ended sequencing libraries
    in FASTQ format.

    Args:
        file_1: File path to read/first mate library.
        file_2: File path to second mate library.
        organism: Source organism of the sequencing library; either a short
            name (string, e.g., `hsapiens`) or a taxon identifier (integer,
            e.g., `9606`).

    Returns:
        LIBTYPE string according to Salmon documentation, cf.
        https://salmon.readthedocs.io/en/latest/library_type.html
    """
    #create .tmp directory in "inputDir" and retuns .tmp path
    tmpDir = mkTmpDir(inputDir)
    
    # call functions of issues #25, 26, 27 and 34
    
    # delete .tmp
    rmTmpDir(tmpDir)
        
    # implement logic
    return "U"
    

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
