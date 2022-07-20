# ALDEx2-python

 Run a popular [ALDEx2 tool](https://github.com/ggloor/ALDEx_bioc) in Python.

### Dependencies

Create a new conda environment and install required dependencies:
```sh
conda create --name pyaldex2 -c conda-forge -c bioconda bioconductor-aldex2 rpy2 pandas
```
Activate the environment:
```sh
conda activate pyaldex2
```

### Run it
Clone this repository to a convinient location to store the script:
```sh
git clone https://github.com/trickovicmatija/ALDEx2-python.git
```
Import the module and set the path
```python
import pandas as pd
import sys, os

sys.path.append("./ALDEx2-python") # Append the location of cloned repository
import pyaldex2

path = f'{os.path.dirname(pyaldex2.__file__)}/run_aldex2.R' # Set the path of the "run_aldex2.R" R-script. Default: same directory as the Python module.
```

Import test data:

```python
counts = pd.read_csv("test_data/counts.tsv,sep='\t')
metadata = pd.read_csv("test_data/metadata.tsv",sep='\t')
```
Run the script:
```python
result_df = pyaldex2.run_aldex2(counts,metadata,'t',r_script_path=path)
```

For the description of ALDEx2 output see [here](https://bioconductor.org/packages/release/bioc/vignettes/ALDEx2/inst/doc/ALDEx2_vignette.html#5_ALDEx2_outputs).