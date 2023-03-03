# readme.md
# Cancer Immunotherapy Data Science Grand Challenge - Challenge 1

~~~
Author: Brody Langille
Date: 2023-01-24
~~~

The instructions below are provided to configure the environment, run the code, and inspect the results. The (2) results .csv files will be written to the `challenge_1/solution/` directory per the online competition instructions:

* `challenge_1/solution/test_output.csv` and,
* `challenge_1/solution/validation_output.csv`

### Locally

To configure the environment and run the code locally follow the steps outlined below.

__Assumptions:__

* Python version >= 3.8 and pip (the Python package manager) are installed locally
* Commands are executed in a Unix-like shell:
    * This was tested in:
        * Terminal
        * Running on MacOS Monterey 12.6.2 (21G320)

__Steps:__

* Unzip and navigate to the provided `challenge_1/code/` folder
* Install the required (Python package) dependencies
    * Run the shell command `pip3 install --no-cache-dir -r requirements.txt`
* The provided data file `sc_training.h5ad` __was too big to upload with our submission__
    * Place that file into the `code/` folder i.e. `code/sc_training.h5ad` __before running any code__
* Run the code
    * Run the shell command `python3 Main.py`
* Inspect the output files written to the `challenge_1/solution/` directory:
    * `challenge_1/solution/test_output.csv` and,
    * `challenge_1/solution/validation_output.csv`
