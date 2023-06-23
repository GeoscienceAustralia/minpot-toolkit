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
Mineral Potential Toolkit
-------------------------
This repository is designed to facilitate statistical appraisal of the 
correlations between geological data / models / features, and classes of known 
mineral deposits. In the future the repository will be developed to enable the 
integration of multiple models together to generate mineral prospectivity 
assessments.

Usage
-----
python main.py --config_file .../config.xml

"""

import toolkit.functions as fn
import toolkit.plotting as pt
import numpy as np
import os, warnings, argparse
from matplotlib.image import imread
import xml.etree.ElementTree as ElementTree
from time import time

from toolkit.feat_eng import FeatureManager
from toolkit.IO.DataManager import DataManager, RandomPointsManager
from toolkit.IO.ModelsManager import ModelsManager
from toolkit.pros_mod import Prospectivity
from toolkit.IO.OutputManager import print_to_file, results_to_file
from toolkit.stat_corr import StatsManager

warnings.filterwarnings("ignore", category = RuntimeWarning)
warnings.filterwarnings("ignore", category = UserWarning)

class MineralPotential():
    def __init__(self, file):
        """
        A class to facilitate mineral prospectivity analysis. Workflow:
            
            1. Parse XML file
            
            2. Read geophysical models
            
            3. Read deposit dataset
            
            4. Generate "features" from geophysical models based on XML config
            
            5. Perform statistical analysis described in XML config
            
            6. Assess prospectivity (needs work)
            
            7. Generate outputs
            
        
        Parameters
        ----------
        file: string
            Name of XML configuration file to use.
        
        
        Calls
        -----
        toolkit.IO.ModelManager.ModelsManager
        
        toolkit.IO.DataManager.DataManager
        
        toolkit.IO.DataManager.RandomPointsManager
        
        
        """
        self.ModelsManager = ModelsManager()
        self.DataManager = DataManager()
        self.RPManager = RandomPointsManager()
        self.bgimage = ""
        self.generate_output = True
        self.make_plots = False
        self.log_plots = False
        self.model_from_savefile = False
        self.model_to_savefile = False
        self.data_from_savefile = False
        self.data_to_savefile = False
        self.stats_flag = False
        self.model_savefile="model_savefile.npy"
        self.domain_savefile="domain_savefile.npy"
        self.data_savefile="data_savefile.npy"
        self.input_path = '.'
        self.output_path = '.'
        
        self.ParseXMLNode(file)
    #end func
    
    def ParseXMLNode(self, file):
        """
        Initialise the problem workflow from xml specification. This function
        will read all commands from the config file and store them in a list to
        be used by all modules required by the workflow.
        
        
        Parameters
        ----------
        file: string
            Name of XML configuration file to use.
            
            
        Calls
        -----
        functions.changetype
        
        IO.ModelsManager.ModelsManager.ParseXMLNode
        
        IO.DataManager.DataManager.ParseXMLNode
        
        IO.DataManager.RandomPointsManager.ParseXMLNode
        
        
        """
        
        rootNode = ElementTree.parse(file).getroot()
        commandsNode = rootNode.find("commands")
        self.commandList = list()
        items = list(commandsNode.attrib.keys())
        for item in items:
            val = commandsNode.attrib[item]
            self.commandList.append([item, val])
        #end for
        for row in self.commandList:
            if hasattr(self, row[0]):
                tp = type(getattr(self, row[0]))
                val = fn.changetype(row[1], tp)
                setattr(self, row[0], val)
            #end if
        #end for
        
        # Input models
        modelsNode = rootNode.find("models")
        self.ModelsManager.ParseXMLNode(modelsNode, self.commandList)
        
        # Input data
        datasetNode = rootNode.find("dataset")
        self.DataManager.ParseXMLNode(datasetNode, self.commandList)
        
        # Input random data
        self.RPManager.ParseXMLNode(self.commandList)
        
        # Input commands for analysis
        featEngNode = rootNode.find("features")
        self.featureNodes = featEngNode.findall("feature")
        if( not (rootNode.find('statistics') is None) ):
            statCorrNode = rootNode.find("statistics")
            self.statnodes = statCorrNode.findall("statsnode")
            self.stats_flag = True
        #end if
    #end func
    
    def Prepare(self):
        """
        Prepares model and deposit data inputs. Based on XML config file, it 
        will either read models and deposit data from an existing numpy binary
        savefile, or will generate new ones.
        
        
        Calls
        -----
        IO.ModelsManager.ModelsManager.Read_Models
        
        IO.DataManager.DataManager.Read_Dataset
        
        IO.DataManager.RandomPointsManager.Prepare
        
        
        """
        print('Reading models')
        if self.model_from_savefile:
            self.ModelsManager.model_dict = \
                np.load(os.path.join(self.input_path, self.model_savefile), 
                        allow_pickle=True).item()
            self.ModelsManager.domain = \
                np.load(os.path.join(self.input_path, self.domain_savefile), 
                        allow_pickle=True).item()
        else:
            self.ModelsManager.Process()
        #end if
        if self.model_to_savefile:
            np.save(os.path.join(self.input_path, self.model_savefile), 
                    self.ModelsManager.model_dict, allow_pickle=True)
            np.save(os.path.join(self.input_path, self.domain_savefile), 
                    self.ModelsManager.domain, allow_pickle=True)
        #end if
        print('Completed reading models, reading deposit datasets')
        if self.data_from_savefile:
            self.DataManager.points_dict = \
                np.load(os.path.join(self.input_path, self.data_savefile), 
                        allow_pickle=True).item()
        else:
            self.DataManager.Read_Dataset(self.ModelsManager.domain)
        #end if
        if self.data_to_savefile:
            np.save(os.path.join(self.input_path, self.data_savefile), 
                    self.DataManager.points_dict, allow_pickle=True)
        #end if
        print('Completed reading deposit datasets, initialising random points')
        self.RPManager.Prepare(self.ModelsManager.domain,
                               self.DataManager.points_dict)
        print('Completed initialising random points')
    #end func
    
    def Run_Analysis(self):
        """
        Runs analysis by generating features and performing statistical 
        computations.
        
        Needs work: prospectivity module.
        
        
        Calls
        -----
        feat_eng.FeatureManager.Process
        
        stat_corr.StatsManager.Process
        
        pros_mod.Prospectivity.Process
        
        
        """
        print('Engineering features')
        self.FeatureManager = FeatureManager(self.commandList, 
                                             self.featureNodes)
        self.FeatureManager.Process(self.ModelsManager.model_dict)
        print('Completed engineering features,',)
        if self.stats_flag:
            print(' running analysis')
            self.StatsManager = StatsManager(self.commandList, self.statnodes)
            self.StatsManager.Process(self.FeatureManager.feature_dict, 
                                      self.DataManager.points_dict, self.RPManager)
            print('Completed running analysis, generating prospectivity maps')
            self.Prospectivity = Prospectivity(self.StatsManager.stat_obj_dict,
                                               self.FeatureManager.feature_dict)
            print('Completed generating prospectivity maps')
        #end if
        print('Finished')
    #end func
    
    def Plots(self):
        """
        Generates plots.
        
        
        Calls
        -----
        plotting.plot_coverage_region
        
        feat_eng.FeatureManager.feature_object.Plots
        
        stat_corr.StatsManager.(various).Plots
        
        pros_mod.Prospectivity.Plots
        
        
        """
        print('Generating plots')
        if len(self.bgimage):
            self.bgimage = imread(self.bgimage)
        #end if
        pt.plot_coverage_region(self.ModelsManager.domain, self.output_path,
                                bgimage=self.bgimage)
        if self.stats_flag:
            for key in self.StatsManager.stat_obj_dict.keys():
                test = self.StatsManager.stat_obj_dict[key]
                path = os.path.join(self.output_path, 'Stats')
                if not os.path.exists(path):
                    os.mkdir(path)
                #end if
                test.Plots(path)
            #end for
        #end if
        for key in self.FeatureManager.feature_dict.keys():
            feat = self.FeatureManager.feature_dict[key]
            path = os.path.join(self.output_path, 'Features')
            if not os.path.exists(path):
                os.mkdir(path)
            #end if
            feat.Plots(self.DataManager.points_dict, path, 
                       bgimage=self.bgimage)
        #end for
        self.Prospectivity.Plots(bgimage=self.bgimage, 
                                 wd_images=self.output_path)
        print('Completed generating plots')
    #end func
    
    def Output(self):
        """
        Function to generate the output of the feature engineering, statistics,
        and prospectivity modelling modules.
        
        
        Calls
        -----
        IO.OutputManager.print_to_file
        
        IO.OutputManager.results_to_file
        
        
        """
        if self.make_plots:
            self.Plots()
        #end if
        if self.stats_flag:
            keys = self.StatsManager.stat_obj_dict.keys()
            for key in keys:
                stat_obj = self.StatsManager.stat_obj_dict[key]
                path = os.path.join(self.output_path, 'Stats')
                if not os.path.exists(path):
                    os.mkdir(path)
                #end if
                results_to_file(stat_obj, path)
            #end for
        #end if
        keys = self.FeatureManager.feature_dict.keys()
        for key in keys:
            feat_obj = self.FeatureManager.feature_dict[key]
            path = os.path.join(self.output_path, 'Features')
            if not os.path.exists(path):
                os.mkdir(path)
            #end if
            print_to_file(feat_obj, path)
        #end for

        try:
            np.save(os.path.join(self.output_path, 'savefile.npy'), self, allow_pickle=True)
        except Exception as e:
            print('Warning: Failed to save model file with error ({})'.format(str(e)))
        # end try
    #end func
    
    def __call__(self):
        """
        Prepares data, runs analysis, and generates plots and output (if 
        flagged).
        
        
        """
        t1 = time()
        # Data and models are initialised during XML parsing
        self.Prepare()
        
        # Calculate statistical relations, and analyse prospectivity
        self.Run_Analysis()
        
        # Outputs and visualisation
        if self.generate_output:
            self.Output()
        #end if
        t2 = time()
        print('Completed in',t2-t1,'seconds.')
    #end func
#end class

if __name__=='__main__':
    parser = argparse.ArgumentParser(description='Mineral Potential Toolkit')
    parser.add_argument("-C", "--config_file", type=str, required=True,
                        default="../config/config.xml")
    args = parser.parse_args()
    file = args.config_file
    #file = os.path.join(r'..\config', 'config.xml')
    ProblemManager = MineralPotential(file)
    ProblemManager()
#end if
