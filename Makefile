all: lint test

lint:
	@python -m compileall -q pipelign bin

release:
	@rm -rf build dist
	@python -m build
	@twine upload dist/*

test:
	@pytest

.PHONY: docs test
