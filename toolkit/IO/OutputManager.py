# Copyright (C) 2021-2023 Geoscience Australia
# 
# The minpot-toolkit is released under the Apache License, Version 2.0 
# (the "License");you may not use this software except in compliance with 
# the License. You may obtain a copy of the License at
# 
#     http://www.apache.org/licenses/LICENSE-2.0
# 
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
# 
# The project uses third party components which may have different licenses. 
# Please refer to individual components for more details.
"""
Output module
-----------------
This module is used to write results from the feature engineering and 
statitsics modules ot text files.
    
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
