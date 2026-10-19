# Loblaw Bio immune cell analysis
#
#   make setup      Create .venv and install the dependencies into it.
#   make pipeline   Build the SQLite database and run Parts 2 to 4.
#   make dashboard  Serve the interactive Streamlit dashboard.
#
# `make pipeline` and `make dashboard` run setup themselves if .venv does not
# exist yet; run `make setup` again if it was interrupted. Everything runs with
# the venv's interpreter; pass PYTHON=... to use another one.

VENV   ?= .venv
PYTHON ?= $(VENV)/bin/python
PORT   ?= 8501

.PHONY: setup pipeline dashboard clean

setup:
	python3 -m venv $(VENV)
	$(VENV)/bin/python -m pip install --upgrade pip
	$(VENV)/bin/python -m pip install -r requirements.txt

# Lets pipeline and dashboard depend on the venv as a file, but only when the
# scripts will actually run on it.
$(VENV)/bin/python:
	$(MAKE) setup

ifeq ($(PYTHON),$(VENV)/bin/python)
pipeline dashboard: | $(VENV)/bin/python
endif

pipeline:
	$(PYTHON) load_data.py
	$(PYTHON) run_analysis.py

dashboard:
	$(PYTHON) -m streamlit run dashboard.py --server.port $(PORT) --server.headless true

clean:
	rm -f cell_counts.db cell_counts.*.tmp
	rm -rf outputs __pycache__ analysis/__pycache__
