This repo is part of the [SOMISANA](https://somisana.ac.za/) initiative, and used for our development around WavewatchIII.

## Installation

The python tools in `ww3_tools` require a suitable environment built off the [wavespectra](https://github.com/wavespectra/wavespectra) library. 

To use the tools, clone this repo to your local machine (do this wherever you would like the code):

```sh
git clone git@github.com:SAEON/somisana-ww3.git
```

Alternatively, if you're not anticipating developing and pushing changes to this repo, you can also use HTTPS protocol (only do one of these clone commands!):

```sh
git clone https://github.com/SAEON/somisana-ww3.git
```

Then navigate to the root directory of the repo:
`cd somisana-ww3`

You can create the `wavespectra` environment by running:

```sh
mamba env create -f environment.yml
```
or use `conda` instead of `mamba` if you haven't moved to mamba yet

Then activate the environment
```sh
conda activate wavespectra
```

and install the python library so it's available in your environment:
```sh
pip install --no-deps -e .
```

If new dependencies get added to `environment.yml`, then you can update your environment by running:

```sh
mamba env update -f environment.yml --prune
```

