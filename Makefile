# Reproduce everything:   make all
# Individual steps:       make test | data | analysis | figures | notebook | paper

PY ?= python

.PHONY: all install test data analysis figures notebook paper clean

all: test data analysis notebook paper

install:
	$(PY) -m pip install -e ".[dev]"

test:
	$(PY) -m pytest

data:
	$(PY) scripts/fetch_data.py
	$(PY) scripts/build_dataset.py

analysis:
	$(PY) scripts/run_analysis.py

figures:            ## rebuild figures and LaTeX macros from existing results (seconds)
	$(PY) scripts/run_analysis.py --reuse

notebook:
	jupyter nbconvert --to notebook --execute --inplace notebooks/walkthrough.ipynb

paper:              ## needs a LaTeX installation
	cd paper && pdflatex -interaction=nonstopmode preview && bibtex preview && \
	pdflatex -interaction=nonstopmode preview && pdflatex -interaction=nonstopmode preview

clean:
	rm -f paper/*.aux paper/*.bbl paper/*.blg paper/*.log paper/*.out
