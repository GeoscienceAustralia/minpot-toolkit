# minpot-toolkit
Tools to facilitate mineral potential analysis, from spatial associations to feature engineering and fully integrated mineral potential mapping.

## Setup

Currently, the environment list to run the minpot-toolkit is difficult to setup. To assist, the repository has been packaged with a preferred environment which can be activated using Conda:

`conda env create -f minpot_env.yml`
`conda activate minpot_env`

Ensure to add the minpot-toolkit directory into the PYTHONPATH user environment variable. For instance, from Windows call (updating "<path_to_minpot-toolkit>"):

`set PYTHONPATH=%PYTHONPATH%;C:\<path_to_minpot-toolkit>`

or for Unix call (updating "<path_to_minpot-toolkit>"):

`export PYTHONPATH=$PYTHONPATH:<path_to_minpot-toolkit>`

To use the environment within the Jupyter Notebook tutorials it will need to be installed as a kernel. With the environment activated, run the following command:

`ipython kernel install --user --name=minpot_env`

## Acknowledgements

©  Commonwealth  of  Australia  (Geoscience  Australia)  2023

Geoscience Australia eCat:  [134193](https://pid.geoscience.gov.au/dataset/ga/XXXX)

Alison Kirkby, Lachlan Adams, Rakib Hassan and Fei Zhang are thanked for development on an earlier version of this repository.

Produced as part of the Exploring for the Future program.

Geoscience Australia’s Exploring for the Future program provides precompetitive information to inform decision-making by government, community and industry on the sustainable development of Australia's mineral, energy and groundwater resources. By gathering, analysing and interpreting new and existing precompetitive geoscience data and knowledge, we are building a national picture of Australia’s geology and resource potential. This leads to a strong economy, resilient society and sustainable environment for the benefit of all Australians. This includes supporting Australia’s transition to net zero emissions, strong, sustainable resources and agriculture sectors, and economic opportunities and social benefits for Australia’s regional and remote communities. The Exploring for the Future program, which commenced in 2016, is an eight year, $225m investment by the Australian Government.

## Contact
Marcus <dot> Haynes <at> ga <dot> gov <dot> au
