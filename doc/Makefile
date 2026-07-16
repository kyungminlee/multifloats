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
VERSION_FILE := ../VERSION
MF_VERSION   := $(shell tr -d '[:space:]' < $(VERSION_FILE))

.PHONY: html doxygen clean

# Substitute the version from the root VERSION file into the .in templates.
conf.py: conf.py.in $(VERSION_FILE)
	sed 's/@MULTIFLOATS_VERSION@/$(MF_VERSION)/g' $< > $@
Doxyfile: Doxyfile.in $(VERSION_FILE)
	sed 's/@MULTIFLOATS_VERSION@/$(MF_VERSION)/g' $< > $@

html: conf.py doxygen
	$(SPHINX) -b html $(SRC) $(BUILD)/html
	@echo ">> $(BUILD)/html/index.html"

doxygen: Doxyfile
	doxygen Doxyfile

clean:
	rm -rf $(BUILD) _doxygen conf.py Doxyfile
