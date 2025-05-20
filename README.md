# Python HPR Detection Toolkit
Toolkit to process grasslands and evaluating their biological value

## About
This package gives the necessary handles to automate the detection of unknown hpr(+) grasslands in Flanders. The packages holds different modules that can process landplots and availble data (like raster maps of height) to evaluate the presence of certain structures within the grassland that are indicative for an old/valuable grassland. The package does select which grasslands need to be inspected but only processes the given landplot. However, in the examples, more elaborate analyses workflows are included.

## Features
At the moment only one type of hpr grassland is included:
* Grasslands holding ditches

## Installing the package
To make sure all dependancies for the package are correct, we will work in a Virtual Environment: It's highly recommended to use a virtual environment to avoid conflicts with other projects. Here we give the instruction to do so using (mini)conda.

1. **Open a (anaconda) terminal**

2. **Create and activate new conda environment**  
Check if the python version is declared correctly (see `pyproject.toml`) and change if needed. Make sure to include pip in the creation of the environment! You can change the name of the new environment from `INBO_microrelief`to something of your liking.
```bash
conda create -n INBO_microrelief python>=3.10 pip
conda activate INBO_microrelief
```

3. **Navigate to the package directory**
```bash
cd /path/to/package/directory
```

4. **Checkout the Git repository**  
Clone the remote git repository into the package directory you navigated to. This will create a local mirror of the repository. Depending on whether you have set-up ssh or not you can use https or ssh.
```bash
git clone https://github.com/inbo/Python_Toolkit_HPR_Detection.git
```
OR
```bash
git clone git@github.com:inbo/Python_Toolkit_HPR_Detection.git
```
After cloning it, go into the repository
```bash
cd Python_Toolkit_HPR_Detection
```

5. **Install the package and dependancies**  
Pip stands for Package Installer for Python. It allows you to easily install, upgrade and manage Python packages. When installing this package, pip will search automatically for dependancies declared in the `pyproject.toml` file and install the necessary packages. Make sure to run pip with the `-e` flag. This makes the project editable and ensure the default config files are found correctly.
```bash
pip install -e .
```

6. **Verify Installation** (Optional)  
You can check if the package and its dependencies were installed correctly
```bash
pip list
```
You should see `hpr_detection_toolkit` listed, pointing to your local source directory. 

## Usage
In the folder `example_notebooks` some examples of the functionality of this package are found. They use large data files such as the BWK, agriculture landusage ... Make sure to download these data sets and point the examples to the right location using absolute paths. Some notebooks also use data that is shipped with this package. 

Do not overwrite the example files or upload new data or processed data files!

To open the notebooks, run the following command in the directory `/path/to/package/directory/Python_Toolkit_HPR_Detection`
```bash
jupyter notebook
```