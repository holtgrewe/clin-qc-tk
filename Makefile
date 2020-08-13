.PHONY: default
default: black

.PHONY: black
black:
	black -l 100 .

.PHONY: test
test:
	pytest .

.PHONY: test-vv
test-vv:
	pytest -vv .

.PHONY: test-cov
test-cov:
	pytest .
	coverage html

.PHONY: pep8
pep8:
	pycodestyle --exclude=versioneer.py --max-line-length=100 --ignore=E203,E266,E501,W503,E231

.PHONY: flake8
flake8:
	flake8
