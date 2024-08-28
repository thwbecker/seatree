#! /user/bin/env python3

import os, sys

def add_subdirectories_to_path(root_dir):
    for dirpath, dirnames, filenames in os.walk(root_dir):
        sys.path.append(dirpath)

# Replace '/path/to/python' with the path to your 'python' folder
root_directory = os.getcwd()
add_subdirectories_to_path(root_directory)

for path in sys.path:
    print(path)

os.system('python3 '+root_directory+'/python/seatree/gui/SEATREE.py')
