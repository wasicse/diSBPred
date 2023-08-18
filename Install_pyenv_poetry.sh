#! /bin/bash

echo "Installing Dependencies"
export PYTHON_KEYRING_BACKEND=keyring.backends.null.Keyring
pyenv install miniconda3-3.9-4.10.3 -y
pyenv local miniconda3-3.9-4.10.3

# Create local poetry environment
rm -rf .venv
rm -rf poetry.lock
python3 -m venv .venv
./.venv/bin/pip install -U pip setuptools
./.venv/bin/pip install poetry
POETRY_VIRTUALENVS_IN_PROJECT="true"

# Install Poetry Dependencies
./.venv/bin/poetry

#Test Installation.venv/bin/poetry run which python
.venv/bin/poetry run python --version