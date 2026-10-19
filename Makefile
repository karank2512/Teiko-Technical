# Loblaw Bio immune cell analysis
#
#   make setup      Create .venv and install the dependencies into it.
#   make pipeline   Build the SQLite database and run the analysis.
#
# `make pipeline` runs setup itself if .venv does not exist yet; run
# `make setup` again if it was interrupted. Everything runs with the venv's
# interpreter; pass PYTHON=... to use another one.

VENV   ?= .venv
PYTHON ?= $(VENV)/bin/python

.PHONY: setup pipeline clean

setup:
	python3 -m venv $(VENV)
	$(VENV)/bin/python -m pip install --upgrade pip
	$(VENV)/bin/python -m pip install -r requirements.txt

# Lets pipeline depend on the venv as a file, but only when the scripts will
# actually run on it.
$(VENV)/bin/python:
	$(MAKE) setup

ifeq ($(PYTHON),$(VENV)/bin/python)
pipeline: | $(VENV)/bin/python
endif

pipeline:
	$(PYTHON) load_data.py
	$(PYTHON) run_analysis.py

clean:
	rm -f cell_counts.db cell_counts.*.tmp
	rm -rf outputs __pycache__ analysis/__pycache__
