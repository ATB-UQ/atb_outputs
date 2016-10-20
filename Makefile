PYTHON_BIN_DIR = /usr/local/python35/bin

PYTHON_EXEC = PYTHONPATH=$(PYTHONPATH) $(PYTHON_BIN_DIR)/python3

errors:
	PYTHONPATH=$(PYTHONPATH) $(PYTHON_BIN_DIR)/pylint -E $$(find . -name '*.py')
.PHONY: errors

mypy: $(PYTHON_BIN_DIR)/mypy
	MYPYPATH=$(PYTHONPATH) $(PYTHON_BIN_DIR)/mypy format.py
.PHONY: mypy
