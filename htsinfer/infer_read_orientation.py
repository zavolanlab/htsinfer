"""Infer read orientation from sample data."""

import logging
import os
import shutil
from tempfile import mkdtemp

LOGGER = logging.getLogger(__name__)


def infer():
    
    """Creates a temporary directory (tmp_dir) in the input 
    directory (input_dir). Returns the path of tmp_dir.
    tmp_dir is used for the storage of temporary data 
    created when performing the main function. 
    After completion of main function, tmp_dir is deleted.
    """
    
    #creates temporary directory in the working directory
    try:
        tmp_dir = mkdtemp(dir=os.getcwd())
    except OSError as exc:
        raise OSError("Creation of temporary directory failed") from exc
    LOGGER.info(f"Created temporary directory {tmp_dir}.")
    
    # Main function in between creation and deletion of tmp_dir
    
    
    # delete temporary directory
    try:
        shutil.rmtree(tmp_dir)
    except OSError as exc:
        raise OSError("Deletion of temporary directory failed") from exc
    LOGGER.info(f"Deleted temporary directory {tmp_dir}.")

    
    

