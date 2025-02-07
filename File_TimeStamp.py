'''
Allows you to extract the timestamp (in s) at which a file was created.
Return is a list.
N.B. only works for the original file
'''

import os, glob
import pathlib

def TimeStamp(path_to_file):
    file_path = glob.glob(path_to_file)
    
    time = []
    
    for file in file_path:
        file_name = pathlib.Path(file)
        timestamp = file_name.stat().st_mtime
        time.append(timestamp)
        
    return(time)
