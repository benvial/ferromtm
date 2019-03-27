SHELL := /bin/bash

.PHONY: clean lint req data

#################################################################################
# GLOBALS                                                                       #
#################################################################################

PROJECT_DIR := $(shell dirname $(realpath $(lastword $(MAKEFILE_LIST))))
PROJECT_NAME = ferromtm
PYTHON_INTERPRETER = python3
HOSTING = github

ifeq (,$(shell which conda))
HAS_CONDA=False
else
HAS_CONDA=True
endif

#################################################################################
# COMMANDS                                                                      #
#################################################################################



## Set up python interpreter environment
env:
ifeq (True,$(HAS_CONDA))
		@echo -e ">>> Detected conda, creating conda environment."
ifeq (3,$(findstring 3,$(PYTHON_INTERPRETER)))
	conda create --name $(PROJECT_NAME) python=3
else
	conda create --name $(PROJECT_NAME) python=2.7
endif
		@echo -e ">>> New conda env created. Activate with:\nsource activate $(PROJECT_NAME)"
else
	$(PYTHON_INTERPRETER) -m pip install -q virtualenv virtualenvwrapper
	@echo -e ">>> Installing virtualenvwrapper if not already intalled.\nMake sure the following lines are in shell startup file\n\
	export WORKON_HOME=$$HOME/.virtualenvs\nexport PROJECT_HOME=$$HOME/Devel\nsource /usr/local/bin/virtualenvwrapper.sh\n"
	@bash -c "source `which virtualenvwrapper.sh`;mkvirtualenv $(PROJECT_NAME) --python=$(PYTHON_INTERPRETER)"
	@echo -e ">>> New virtualenv created. Activate with:\nworkon $(PROJECT_NAME)"
endif

## Test if python environment is setup correctly
testenv:
	source activate $(PROJECT_NAME); \
	$(PYTHON_INTERPRETER) testenv.py



## Install Python Dependencies
req: testenv
	source activate $(PROJECT_NAME); \
	$(PYTHON_INTERPRETER) -m pip install -U pip setuptools wheel; \
	$(PYTHON_INTERPRETER) -m pip install -r requirements.txt



## Create a github/gitlab repository
repo:
	git init
ifeq (github,$(HOSTING))
	hub create $(PROJECT_NAME)
else
	python gitlab_create.py
	git init
	git remote add origin git@gitlab.com:qmaem/$(PROJECT_NAME).git
endif


## Initialize the project: create environment, install dependencies, and create repository on github
init: env req repo

## Delete all compiled Python files
clean:
	find . -type f -name "*.py[co]" -delete
	find . -type d -name "__pycache__" -delete

## Lint using flake8
lint:
	flake8 ferromtm

## Reformat code
style:
	@echo "Styling..."
	black ferromtm


## Push to github
gh:
	@echo "Pushing to github..."
	git add -A
	@read -p "Enter commit message: " MSG; \
	git commit -a -m "$$MSG"
	git push origin master

## Clean, reformat and push to github
save: clean style gh

## Make doc css
less:
	cd docs/_custom/static/css/less;\
	lessc theme.less  ../theme.css;\
	lessc custom_styles.less  ../custom_styles.css;\
	lessc custom_gallery.less  ../custom_gallery.css;\
	lessc custom_pygments.less  ../custom_pygments.css;\

## Make html doc with Sphinx
doc:
	cd docs && make clean && sphinx-apidoc -f -o . ../$(PROJECT_NAME) && make html

## Style and build html doc
styledoc: less doc


## run the test suite
test:
	pytest ./tests -s --cov=./$(PROJECT_NAME)


## Run the codes for fitting BST measurements to a model
bst:
	python ./ferromtm/models/bst.py

## Run the codes for periodic rods
circ:
	python ./ferromtm/models/run_circ.py

## Run the codes for periodic rods convergence
conv:
	python ./ferromtm/models/run_convergence.py

## Run the codes for random rods
rand:
	python ./ferromtm/models/run_rand.py

## Run the codes and make the datasets
data: bst circ conv rand


#################################################################################
# PROJECT RULES                                                                 #
#################################################################################



#################################################################################
# Self Documenting Commands                                                     #
#################################################################################

.DEFAULT_GOAL := help

# Inspired by <http://marmelab.com/blog/2016/02/29/auto-documented-makefile.html>
# sed script explained:
# /^##/:
# 	* save line in hold space
# 	* purge line
# 	* Loop:
# 		* append newline + line to hold space
# 		* go to next line
# 		* if line starts with doc comment, strip comment character off and loop
# 	* remove target prerequisites
# 	* append hold space (+ newline) to line
# 	* replace newline plus comments by `---`
# 	* print line
# Separate expressions are necessary because labels cannot be delimited by
# semicolon; see <http://stackoverflow.com/a/11799865/1968>
.PHONY: help
help:
	@echo -e "$$(tput bold)Available rules:$$(tput sgr0)"
	@echo -e
	@sed -n -e "/^## / { \
		h; \
		s/.*//; \
		:doc" \
		-e "H; \
		n; \
		s/^## //; \
		t doc" \
		-e "s/:.*//; \
		G; \
		s/\\n## /---/; \
		s/\\n/ /g; \
		p; \
	}" ${MAKEFILE_LIST} \
	| LC_ALL='C' sort --ignore-case \
	| awk -F '---' \
		-v ncol=$$(tput cols) \
		-v indent=19 \
		-v col_on="$$(tput setaf 6)" \
		-v col_off="$$(tput sgr0)" \
	'{ \
		printf "%s%*s%s ", col_on, -indent, $$1, col_off; \
		n = split($$2, words, " "); \
		line_length = ncol - indent; \
		for (i = 1; i <= n; i++) { \
			line_length -= length(words[i]) + 1; \
			if (line_length <= 0) { \
				line_length = ncol - indent - length(words[i]) - 1; \
				printf "\n%*s ", -indent, " "; \
			} \
			printf "%s ", words[i]; \
		} \
		printf "\n"; \
	}' \
	| more $(shell test $(shell uname) = Darwin && echo '--no-init --raw-control-chars')
