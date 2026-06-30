# MovieLens 1M Data Application


## Directory Structure

- **Code (`code/movielens/`):** Contains all R and Python scripts required to execute the models and generate the results.
- **Data (`data/movielens/raw/`):** Contains the raw MovieLens datasets. 
- **Output (`output/movielens/`):** Contains the model outputs (fit files).

## Workflow
- Requirements:
    - R packages: `IMR`, `SoftImpute`, `tidyverse`, `magrittr`, `Matrix`, `kableExtra`, `RSSthemes`, `ggh4x`

To fully reproduce the MovieLens results, execute the scripts in the following order.

### 1. Model Fitting
- **`1_fit_MovieLens.R`**: Fits the proposed IMR methods (IMR-I and IMR-IXZ) as well as the `SoftImpute` and `MCAI` models. All outputs are saved to `output/movielens/model_fits/`.

- **Configuration:** `run_cross_validation` in `config_default.R` (default `FALSE`)
controls hyperparameter selection. `FALSE` uses the fixed hyperparameters from the
paper; `TRUE` re-runs cross-validation. Leave it `FALSE`.


### 2. Deep Learning Method (Glocal-K)
#### 2.1 Data preparation
- **`2_1_prepare_python_data.R`**: Prepares the data for the GLocal-K model.
#### 2.2 Model Fitting 
There are two implementations of the GLocal-K model on the MovieLens dataset. Run one of them.
**Option 1 (PyTorch Implementation. Easy to Run):**
- **`2_2_GlocalK_torch.py`:** This is a recent unofficial PyTorch implementation (https://github.com/fleanend/TorchGlocalK). As noted by the authors, this version performs worse than the original implementation below. However, it does not require the creation of a virtual environment.
- Requirements:
    - Python >3.10
    - torch
    - numpy, scipy, pandas

- **Configuration:** Both scripts read epoch counts from an optional `config.py`
(`epoch_p`, `epoch_f`, `iter_p`, `iter_f` for TensorFlow; `max_epoch_p`, `max_epoch_f`
for PyTorch); without it they use the authors' defaults. Each writes its results to
`output/movielens/model_fits/glocalk_results.csv`, read by `3_generate_results_table.R`.

**Option 2 (Original Implementation. Matches reported results):**
- **`2_2_GlocalK_tensorflow.py`:** This is the original implementation by the authors (https://github.com/usydnlp/Glocal_K/tree/main). However, it requires TensorFlow 1.15 and Python 3.7, and hence must be run inside a Docker container. We provide instructions for creating the Docker container and running the script inside it. It has been tested on macOS and Linux. For Windows, we recommend using WSL: https://learn.microsoft.com/en-us/windows/wsl/install.

- Requirements:
    - Python 3.7
    - tensorflow==1.15.5
    - numpy, scipy, pandas
- Creating the container:
    - Dockerfile: `Dockerfile`
    - ```bash
        docker build --platform linux/amd64 -t tf1-glocalk -f code/movielens/Dockerfile .
      ```
- Running the script:
    - ```bash
        docker run --rm --platform linux/amd64 -v "$(pwd)":/project tf1-glocalk python code/movielens/2_2_GlocalK_tensorflow.py
      ```
    

### 3. Result Generation
- **`3_generate_results_table.R`**: Once all model fits have been saved to disk, run this script. It reads the files from `output/movielens/model_fits/` and generates the results tables and plot shown in the manuscript. All final artifacts are written to `output/movielens/`.

## Running the analysis:
Make sure you are in the project root directory and that all requirements are installed.
```bash
Rscript code/movielens/1_fit_MovieLens.R
Rscript code/movielens/2_1_prepare_python_data.R
python code/movielens/2_2_GlocalK_torch.py  
#or
# docker build --platform linux/amd64 -t tf1-glocalk -f code/movielens/Dockerfile .
# docker run --rm --platform linux/amd64 -v "$(pwd)":/project tf1-glocalk python code/movielens/2_2_GlocalK_tensorflow.py
Rscript code/movielens/3_generate_results_table.R
```