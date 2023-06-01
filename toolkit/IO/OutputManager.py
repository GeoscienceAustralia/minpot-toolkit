"""
Output module
-----------------
This module is used to write results from the feature engineering and 
statitsics modules ot text files.


Contacts
--------
Lachlan Adams
    - Lachlan.Adams@ga.gov.au or Lachlan.Adams.1996@outlook.com
    
    
"""

import os, sys

def results_to_file(stat_obj, path):
    name = stat_obj.name
    write_path = os.path.join(path, name)
    if not os.path.exists(write_path):
        os.mkdir(write_path)
    #end if
    filename = os.path.join(write_path, 'output.txt')
    if os.path.exists(filename):
        os.remove(filename)
    #end if
    stdout = sys.stdout
    sys.stdout = open(filename, 'w')
    stat_obj.print_results()
    #sys.stdout.close() #Runs forever?
    sys.stdout = stdout
#end func
    
def print_to_file(feat_obj, path):
    name = feat_obj.name
    write_path = os.path.join(path, name)
    if not os.path.exists(write_path):
        os.mkdir(write_path)
    #end if
    filename = os.path.join(write_path, 'output.txt')
    if os.path.exists(filename):
        os.remove(filename)
    #end if
    stdout = sys.stdout
    sys.stdout = open(filename, 'w')
    feat_obj.print_description()
    #sys.stdout.close() #Runs forever?
    sys.stdout = stdout
#end func