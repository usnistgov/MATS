Basic `MATS` installation instructions.  Note that development is alpha.  In the future, full pypi/conda installs will be made available.

# Create a conda environment with MATS dependencies

Download the [environment](environment.yml) file.  Note that without setting the environment name below,
by default an evironment `MATS-env` will be created.


``` shell
conda env create -n {optional-name-of-environment} -f environment.yml
```

Alternatively, we provide a conda metapackage `mats-dependencies`, which includes all dependencies for `MATS`` (excluding python).
This can be installed using

``` shell
conda install -n {optional-name-of-environment} -c wpk-nist mats-dependencies
```

# From source

After cloning the repo, you can do the following.

``` shell
pip install .
```

To install an editable version, use option `-e`.  To exclude dependencies, use option `--no-deps`

# With pip from github

This requires git also be installed.  Downside is that the whole repo (including all examples) are clones.

``` shell
pip install git+https://github.com/wpk-nist-gov/MATS.git@feature/master-reformat
```

# With pip from github using wheel
Note, this is experimental.  Do the following


``` shell
pip install https://raw.githubusercontent.com/wpk-nist-gov/MATS/feature/master-reformat/wheel/MATS-3-py3-none-any.whl
```

(Note that the actual version will not be 3)

Alternatively, download the wheel and run

``` shell
pip install path-to-wheel.whl
```


# Install everything with pip

(Still testing)

Just run the pip install command to install everything (sans python) with pip.
