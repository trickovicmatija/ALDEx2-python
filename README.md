# ALDEx2-python

 Run a popular [ALDEx2 tool](https://bioconductor.org/packages/release/bioc/html/ALDEx2.html) in Python via [rpy2](https://github.com/rpy2/rpy2).

 This tool is aware of the compositionality of sequencing data, and treats them accordingly.

### Dependencies

Create a new conda environment and install required dependencies:
```sh
conda create --name pyaldex2 -c conda-forge -c bioconda bioconductor-aldex2 rpy2 pandas
```
Activate the environment:
```sh
conda activate pyaldex2
```

## Run it
Clone this repository to a convinient location to store the script:
```sh
git clone https://github.com/trickovicmatija/ALDEx2-python.git
```
Import the module and set the path:
```python
import pandas as pd
import sys, os

sys.path.append("./ALDEx2-python") # Append the location of cloned repository
import pyaldex2

path = f'{os.path.dirname(pyaldex2.__file__)}/run_aldex2.R' # Set the path of the "run_aldex2.R" R-script. Default: same directory as the Python module.
```

Import test data:

```python
counts = pd.read_csv("test_data/raw_counts.tsv",sep='\t',index_col=0) # It will automatically orient the dataframe
metadata = pd.read_csv("test_data/metadata.tsv",sep='\t',index_col=0)
```
Run the script:
```python
result_df = pyaldex2.run_aldex2(counts,metadata,'t',r_script_path=path)
```
For the description of ALDEx2 output see [here](https://bioconductor.org/packages/release/bioc/vignettes/ALDEx2/inst/doc/ALDEx2_vignette.html#5_ALDEx2_outputs).

## Plot

Before plotting, make sure neccessary libraries (numpy, seaborn, matplotlib.pyplot) are installed.

If not, install them in order to continue.

You can create a **MA-plot** using ```pyaldex2.MA_plot()```:

```python
pyaldex2.MA_plot(result_df,effect_threshold=1)
```
Check [help of the function](https://github.com/trickovicmatija/ALDEx2-python/blob/d7dcc5a0a342f61e641f1e89c20c5ead3ffcbaea/pyaldex2.py#L57) for details.

Also, you can create a **vulcano-plot** using ```pyaldex2.vulcano_plot()```:

```python
pyaldex2.vulcano_plot(result_df)
```
Check [help of the function](https://github.com/trickovicmatija/ALDEx2-python/blob/b00bcf32d92a984d7dc279e9a85c24151cc00ff5/pyaldex2.py#L92) for details.

# Disclaimer

This script is only running the original [ALDEx2 tool](https://bioconductor.org/packages/release/bioc/html/ALDEx2.html) in Python.
I don't take any responsibility or credits from it.
For any questions about the tool itself, please write to them.