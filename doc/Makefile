# Documentation build for multifloats (Doxygen XML -> Sphinx/Breathe HTML).
#
# One-time setup:
#   uv venv .venv
#   uv pip install --python .venv -r requirements.txt
#   # plus the `doxygen` system package on PATH
#
# Then:  make html

VENV    := .venv
SPHINX  := $(VENV)/bin/sphinx-build
SRC     := .
BUILD   := _build

.PHONY: html doxygen clean

html: doxygen
	$(SPHINX) -b html $(SRC) $(BUILD)/html
	@echo ">> $(BUILD)/html/index.html"

doxygen:
	doxygen Doxyfile

clean:
	rm -rf $(BUILD) _doxygen
