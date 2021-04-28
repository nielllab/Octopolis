"""
paths.py
"""
import os, fnmatch

def find(pattern, path):
    result = [] # initialize the list as empty
    for root, dirs, files in os.walk(path): # walk though the path directory, and files
        for name in files:  # walk to the file in the directory
            if fnmatch.fnmatch(name,pattern):  # if the file matches the filetype append to list
                result.append(os.path.join(root,name))
    return result # return full list of file of a given type

def list_subdirs(root_dir):
    dirnames = []
    for root, dirs, files in os.walk(root_dir):
        for rec_dir in dirs:
            dirnames.append(rec_dir)
    return dirnames