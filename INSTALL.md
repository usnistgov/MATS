Basic `MATS` installation instructions.

# Create a conda environment with everything you need

Download the [environment](environment.yml) file.  Note that without setting the environment name below,
by default an evironment `MATS-env` will be created.


``` shell
conda create -n {optional-name-of-environment} -f environment.yml
```

# From source

After cloning the repo, you can do the following.

``` shell
pip install .
```

To install an editable version, use option `-e`.  To exclude dependencies, use option `--no-deps`

# With pip from github

This requires git also be installed.

``` shell
pip install git+https://github.com/wpk-nist-gov/MATS.git@feature/master-reformat
```
