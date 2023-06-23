# minpot-toolkit
The mineral potential toolkit (aka minpot-toolkit) provides tools to facilitate mineral potential analysis, from spatial associations to feature engineering and fully integrated mineral potential mapping.

## Setup

Currently, the environment list to run the minpot-toolkit is difficult to setup. To assist, the repository has been packaged with a preferred environment which can be activated using Conda:

`conda env create -f minpot_env.yaml`

`conda activate minpot_env`

Ensure to add the minpot-toolkit directory into the PYTHONPATH user environment variable. For instance, from Windows call (updating "<path_to_minpot-toolkit>"):

`set PYTHONPATH=%PYTHONPATH%;C:\<path_to_minpot-toolkit>`

or for Unix call (updating "<path_to_minpot-toolkit>"):

`export PYTHONPATH=$PYTHONPATH:<path_to_minpot-toolkit>`

To use the environment within the Jupyter Notebook tutorials it will need to be installed as a kernel. With the environment activated, run the following command:

`ipython kernel install --user --name=minpot_env`

After this you can navigate to the examples directory and launch Jupyter Notebooks (nb: ensure to select the "minpot_env" kernel within the Notebook to access the correct environment settings):

`cd <path_to_minpot-toolkit>/examples`

`jupyter notebook`

## Citation

Users of the minpot-toolkit are asked to please cite: 

Haynes, M.W., Adams, L. and Ford, A. (2023) Mineral Prospectivity as an Inference Question: Towards the Unification of ‘Knowledge’- and ‘Data’-driven Approaches. Proceedings of the World Mining Congress 2023, 26-29 June, Brisbane, Australia.

## Acknowledgements

©  Commonwealth  of  Australia  (Geoscience  Australia)  2023

Geoscience Australia eCat:  [147993](https://pid.geoscience.gov.au/dataset/ga/147993)

Alison Kirkby, Lachlan Adams, Rakib Hassan and Fei Zhang are thanked for development on earlier versions of this repository. This repository is a product of the Exploring for the Future program. Geoscience Australia’s Exploring for the Future program provides precompetitive information to inform decision-making by government, community and industry on the sustainable development of Australia's mineral, energy and groundwater resources. By gathering, analysing and interpreting new and existing precompetitive geoscience data and knowledge, we are building a national picture of Australia’s geology and resource potential. This leads to a strong economy, resilient society and sustainable environment for the benefit of all Australians. This includes supporting Australia’s transition to net zero emissions, strong, sustainable resources and agriculture sectors, and economic opportunities and social benefits for Australia’s regional and remote communities. The Exploring for the Future program, which commenced in 2016, is an eight year, $225m investment by the Australian Government.

## Contact
Marcus.Haynes@ga.gov.au
