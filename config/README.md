# Input data

All user inputs are facilitated using the config.xml file in this directory, which will overwrite any default values for any module in the toolkit.

A user must provide, along with any input commands:
- Input geophysical model.
- Input deposit dataset.
- Description of features to analyse.
- Description of statistical tests to perform.

## Input commands

Enter each of the below as a string. 
For boolean values, enter as "0" or "False" for false, or "1" or "True" for true.
For lists or arrays, enter as a string with a comma as separator.

- generate_output: boolean
	Command to tell the toolkit to write output to files, or simply store it in memory.
- input_path: string
	Directory for model and deposit dataset information to be saved in after initial processing if "model_to_savefile" and "data_to_savefile" are true, or directory to read model and deposit dataset information from if "model_from_savefile" and "data_from_savefile" are true.
- output_path: string
	Directory in which to save outputs.
- inside_polygon_file: string
	Name of shapefile to use as a polygon that deposits or random points must be within.
- outside_polygon_file: string
	Name of shapefile to use as a polygon that deposits or random points must be outside of.
- model_savefile: string
	File name for input geophysical models.
- domain_savefile: string
	Name of file to save the coverage region polygon to.
- data_savefile: string
	File name for input deposit datasets.
- bgimage: string
	File name for background image to use for spatial plots (if defined).
- make_plots: boolean
	Command to generate output plots.
- model_from_savefile: boolean
	Command to read geophysical models from savefile if they have been generated in a previous run of the algorithm.
- data_from_savefile: boolean
	Command to read deposit datasets from a savefile if they have been processed in a previous run of the algorithm.
- model_to_savefile: boolean
	Command to write models produced in the current run to a savefile.
- data_to_savefile: boolean
	Command to write deposit dataset processed in the current run to a savefile.
- read_columns: list
	Command to tell the algorithm which columns from deposit datasets are required.
- point_properties: list
	Command to tell the algorithm which columns from the deposit dataset are to be stored/saved in a savefile.
- bounds: array
	Command which gives a general area for the algorithm to expect to process models and deposit datasets for. Enter as minimum longitude, maximum longitude, minimum latitude, maximum latitude.
- cellsize: float
	Size of a grid square for the model interpolate to grid function.


## Input deposit datasets

The input/output module can currently handle deposit datasets from .csv, .xls, .xlsx, or .txt files (as long as .txt file is in the format output by numpy savetxt).
Each different input file is to be described by a different xml node. Separate sheets in an excel document also need separate nodes.

If the dataset is to be read from a .txt or .csv file, the following attributes are needed for the "data" node:
- name
- columns: list of the names of columns which are present in the file, in order. Must include "longitude", and "latitude".
- file_name
- path

If the dataset is to be read from an excel document or a shapefile, the following attributes are needed for the "data" node:
- name
- sheet name
- lon_col_name: the column name in which the longitude coordinate of each deposit is stored.
- lat_col_name: the column name in which the latitude coordinate of each deposit is stored.
- file_name
- path

The dataset will be filtered by various properties input by the user. These are each defined in a separate node.
Each "property" node must have the following attributes:
- name: a unique identifier for the truth array generated by the function which checks each property.
- description: list with three elements.
	- 1st element: name of column, or name of property.
	- 2nd element: condition relating the 1st and 2nd elements.
	- 3rd element: name of column, name of property, or value.
	
The algorithm will then filter the dataset based on another node. The "filter" node must have the following attributes:
- filters: list of properties to be used.
- condition: either "union" or "intersection".

Examples of conditions which may be applied are as follows:
- description = "Cu, >, 0": column "Cu" must have a value larger than 0.
- description = "Deposit group, ==, Iron oxide copper-gold": the entry in column "Deposit group" must match the string "Iron oxide copper-gold".
- description = "Cu, >, Au": column "Cu" must have a larger value than column "Au".
- description = "property1, and, property2": two previously defined properties "property1" and "property2" must both be true.


## Input models

The input/output module can currently handle geophysical models which are in various formats, explained below.
Each model is to be described by a different xml node. 
Every model node must have the following attributes:

- name: string
	Keyword identifying the name of the model.
- type: string
	Keyword identifying the type of geophysical model
- model_format: string
	Keyword identifying the format of the input model. Currently supported are 'ModEM', 'MagGrav', 'netcdf', 'geotiff', 'csv'
- model_file_name: string
	Name of file containing the model.
- path: string
	Directory the model is stored in.

Every model node, except the ones for geotiff models, must have one of the following:
- epsg: integer 
	Value representing the EPSG code for the model projection.
- proj_str: string 
	String defining the model projection, in the format required by pyproj.
	
If the model is saved in the ModEM format, the model node must also have the following attributes:
- centre: array
	(x, y) coordinates on which the model is centred.
- station_buffer: float 
	Value describing a threshold distance in degrees allowable from a recording station. Points outside of the threshold will be assigned NaN values.
- One of the following:
	- stationxy_fn: string
		Name of file containing list of station coordinates in the projection specified by epsg or proj_str.
	- stationll_fn: string
		Name of file containing list of station coordinates as latitude and longitude.

If the model is saved in the comma separated values format, the model node must also have the following attributes:
- dims: array 
	Array of dimensions of the model. (m, n, o) is interpreted as 'm' x samples, 'n' y samples, and 'o' z samples.
- cellsize: array 
	Array representing the cell spacing in metres between grid points. Format is (dx, dy).
- bottomleft: array 
	Coordinates representing the bottom left corner of the model. Format is (x, y, z), all in metres. Positive z is a depth below the surface.
- depth_interval: float
	Spacing between depth slices in metres.
- nanval: float 
	Value which represents missing data or NaN values.
	
If the model is saved in the netcdf format, the model node must also have the following attributes:
- nanval: float 
	Value which represents missing data or NaN values.
- x_var_name: string
	The name of the variable which represents x values.
- y_var_name: string
	The name of the variable which represents y values.
- values_var_name: string
	The name of the variable which represents data values.
- Optional:
	- z_var_name: string
		If the model is 3D, this is the name of the variable which represents z values.
	- swap_xy: boolean
		If true, the data values array will be transposed in the (x, y) plane.
	- interval: integer
		If the model has a fine resolution, this variable tells the toolkit to take every 'i'th data point. For example, if i=2, the toolkit takes every second x, y, z, and data value, decreasing the memory usage by 8 times.

If the model is saved in the csv format specific to James Goodwin's Magnetic/Gravity inversions, the model node must also have the following attributes:
- dims: array 
	Array of dimensions of the model. (m, n) is interpreted as 'm' x samples, and 'n' y samples.
- cellsize: array 
	Array representing the cell spacing in metres between grid points. Format is (dx, dy).
- bottomleft: array 
	Coordinates representing the bottom left corner of the model. Format is (x, y, z), all in metres. Positive z is a depth below the surface.
- nzlst: list
	List of the number of depth slices at the corresponding spacing in dzlst.
- dzlst: list
	List of depth intervals (metres) for irregular depth spacing.
- nanval: float 
	Value which represents missing data or NaN values.


## Features

The user must define a series of "features" to be generated from the geophysical models, in order to perform statistical analysis on the models. A separate node must be created for each feature the user wishes to generate.

The "feature" node must have the following attributes:
- name: string
	Unique identifier for the feature.
- model_types
	List of model types/existing features required to generate the feature.
- Several "criterion" nodes.
- Optional:
	- use_existing_features: boolean
		Flag to tell toolkit that at least one of the criteria is to be fulfilled using an already existing model, for example, we may wish to find the area which is less than a threshold distance of an existing feature.
	- generate_plots: boolean
		Flag to tell the toolkit whether or not to generate an image of this feature.

Each "criterion" node must contain the following:
- name: string
	Unique identifier for the criterion.
- description: list
	List of 3 elements describing the criterion to be satisfied.
- One of the following:
	- model_type: string
		The type of model to use to fulfil the criterion.
	- feature_name: string
		The name of the already generated feature to be used to fulfil the criterion.
- Optional:
	- project_to_surface: list
		List in the form (min_z, max_z). If this attribute exists then the toolkit looks down the z axis and projects any true values between these depth thresholds to a surface, and creates a 2D feature.
		If there is at least one criterion which does not have the attribute project_to_surface, the toolkit will then project this surface back onto the original grid to generate a 3D feature.

Some example descriptions for the criteria are as follows.
- description = "Value, less than, 100"
- description = "Gradient, greater than, 10"
- description = "Distance, less than, 1000"


## Statistical calculations

A separate node must be created for each statistical calculation the user wishes to perform.
Currently supported statistical tests are:
- Binomial p test - calculates the probability that an observed distribution of deposits could occur if points generated at random are assumed to follow the binomial distribution.
- Poisson p test - as with the binomial p test, for the Poisson distribution.
- Kolmogorov-Smirnov p test - Calculates the probability that an observed difference between the probability density function of the deposit datasets and random locations with respect to some variable is statistically significant.
For the Binomial and Poisson p tests, the chi-squared goodness of fit test it also performed to determine if the randomly generated points do follow the assumed distribution.

Each "statistics" node must have the following attributes:
- name: string
	Unique identifier for the test being performed.
- test: string
	The name of the test to perform, from the list above.
- depth_type: string
	Either "values", representing that the "depths" attribute represents a list of values, or "range", indicating that the "depths" attribute represents the parameters to use for the numpy.arange function.
- depths: array
	Array of values. Either in the form of an array of depth values to perform the statistical test at (depth_type = values), or in the form (minimum depth, maximum depth, spacing) for the numpy arange function (depth_type = range).
- n_repeats: integer
	The number of iterations of the analysis to perform for a randomly generates set of locations, for comparison to real deposit locations.

If the test to be performed is the Kolmogorov-Smirnov p test, then the statistics node must also have the following attributes:
- binsize: float
	The width of a bin for computation of the histogram, to generate probability density functions and cumulative distribution functions.
- correlation: string
	Either "Positive", "Negative", or "Neutral", to flag whether we expect a positive or negative correlation, or neither.
- feature: string
	Unique identifier for the feature analysis is to be performed on.
- Optional:
	- pros_map: boolean
		If true, for each grid point in the geophysical model, the distance to the feature identified above will be computed.
		Then, a contour plot is generated, with values representing the proportion of deposits which fell at a distance less than or equal to that of the grid point from the feature.
		
If the test to be performed is the Binomial or Poisson p test, then each statistics node must also have the following:
- Several criterion nodes

Each criterion node must have the following attributes:
- name: string
	Unique identifier for the criterion.
- feature: string
	Name of the feature to be considered.
- description: list
	A list with three elements describing the criterion.
	
Example descriptions are:
description = "Distance, less than, 0": the toolkit will find the points which fall within the specified feature.
description = "Distance, less than, 10000": the toolkit will find the points which fall less than 10km from the specified feature.
description = "Value, greater than, 0.5": as the feature is represented as a list of boolean values, this will also find the points which fall within the specified feature.