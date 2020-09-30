# Inspired by https://pipenv.pypa.io/en/latest/advanced/#travis-ci.

.PHONY: all init test

all: init test

init:
        # Use --site-packages + PIP_IGNORE_INSTALLED "to make sure that any
        # pip-installable dependencies actually are installed into the virtual
        # environment, with the system site-packages then being used only to
        # resolve imports that aren't covered by the Python level dependency
        # graph at all" (in our case, that's PyGTK). See
        # https://github.com/pypa/pipenv/issues/748#issue-261233053
	pipenv --two --site-packages && PIP_IGNORE_INSTALLED=1 pipenv install

test:
	pipenv run run_tests.py
