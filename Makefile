.PHONY: clean clean-test clean-pyc clean-build help
.DEFAULT_GOAL := help

define BROWSER_PYSCRIPT
import os, webbrowser, sys

from urllib.request import pathname2url

webbrowser.open("file://" + pathname2url(os.path.abspath(sys.argv[1])))
endef
export BROWSER_PYSCRIPT

define PRINT_HELP_PYSCRIPT
import re, sys

for line in sys.stdin:
	match = re.match(r'^([a-zA-Z_-]+):.*?## (.*)$$', line)
	if match:
		target, help = match.groups()
		print("%-20s %s" % (target, help))
endef
export PRINT_HELP_PYSCRIPT

BROWSER := python -c "$$BROWSER_PYSCRIPT"

help:
	@python -c "$$PRINT_HELP_PYSCRIPT" < $(MAKEFILE_LIST)

clean: clean-build clean-pyc clean-test ## remove all build, test, coverage and Python artifacts

clean-build: ## remove build artifacts
	rm -fr build/
	rm -fr dist/
	rm -fr .eggs/
	find . -name '*.egg-info' -exec rm -fr {} +
	find . -name '*.egg' -exec rm -f {} +

clean-pyc: ## remove Python file artifacts
	find . -name '*.pyc' -exec rm -f {} +
	find . -name '*.pyo' -exec rm -f {} +
	find . -name '*~' -exec rm -f {} +
	find . -name '__pycache__' -exec rm -fr {} +

clean-test: ## remove test and coverage artifacts
	rm -fr .tox/
	rm -f .coverage
	rm -fr htmlcov/
	rm -fr .pytest_cache




################################################################################
# utilities
################################################################################
.PHONY: lint pre-commit-init pre-commit-run pre-commit-run-all init

lint: ## check style with flake8
	flake8 MATS tests

pre-commit-init: ## install pre-commit
	pre-commit install

pre-commit-run: ## run pre-commit
	pre-commit run

pre-commit-run-all: ## run pre-commit on all files
	pre-commit run --all-files

.git: ## init git
	git init

init: .git pre-commit-init ## run git-init pre-commit


################################################################################
# virtual env
################################################################################
.PHONY: conda-env conda-dev conda-all mamba-env mamba-dev mamba-all activate
conda-env: ## conda create base env
	conda env create -f environment.yml

conda-dev: ## conda update development dependencies
	conda env update -n MATS-env -f environment-dev.yml

conda-all: conda-env conda-dev ## conda create development env

mamba-env: ## mamba create base env
	mamba env create -f environment.yml

mamba-dev: ## mamba update development dependencies
	mamba env update -n MATS-env -f environment-dev.yml

mamba-all: mamba-env mamba-dev ## mamba create development env

activate: ## activate base env
	conda activate MATS-env


################################################################################
# my convenience functions
################################################################################
.PHONY: user-venv user-autoenv-zsh user-all
user-venv: ## create .venv file with name of conda env
	echo MATS-env > .venv

user-autoenv-zsh: ## create .autoenv.zsh files
	echo conda activate MATS-env > .autoenv.zsh
	echo conda deactivate > .autoenv_leave.zsh

user-all: user-venv user-autoenv-zsh ## runs user scripts


################################################################################
# Testing
################################################################################
.PHONY: test test-all coverage
test: ## run tests quickly with the default Python
	pytest -x -v

test-all: ## run tests on every Python version with tox
	tox -- -x -v

coverage: ## check code coverage quickly with the default Python
	coverage run --source MATS -m pytest
	coverage report -m
	coverage html
	$(BROWSER) htmlcov/index.html


version: ## check version of package
	python setup.py --version


################################################################################
# Docs
################################################################################
.PHONY: docs serverdocs
docs: ## generate Sphinx HTML documentation, including API docs
	rm -f docs/MATS.rst
	rm -f docs/modules.rst
	rm -fr docs/generated
	$(MAKE) -C docs clean
	$(MAKE) -C docs html
	$(BROWSER) docs/_build/html/index.html

servedocs: docs ## compile the docs watching for changes
	watchmedo shell-command -p '*.rst' -c '$(MAKE) -C docs html' -R -D .


################################################################################
# distribution
################################################################################
dist: ## builds source and wheel package (run clean?)
	python setup.py sdist
	python setup.py bdist_wheel
	ls -l dist

.PHONY: release release-test conda-dist
release: dist ## package and upload a release
	twine upload dist/*

release-test: dist ## package and upload to test
	twine upload --repository testpypi dist/*

conda-dist: ## build conda dist (run dist and clean?)
	mkdir conda_dist; \
	cd cond_dist; \
	grayskull pypi MATS ; \
	conda-build .; \
	echo 'upload now'


################################################################################
# installation
################################################################################
.PHONY: install install-dev
install: ## install the package to the active Python's site-packages (run clean?)
	python setup.py install

install-dev: ## install development version (run clean?)
	pip install -e . --no-deps
