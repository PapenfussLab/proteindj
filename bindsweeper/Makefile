.PHONY: help test test-fast test-cov lint format clean install dev-install

help:  ## Show this help message
	@echo "Available commands:"
	@grep -E '^[a-zA-Z_-]+:.*?## .*$$' $(MAKEFILE_LIST) | sort | awk 'BEGIN {FS = ":.*?## "}; {printf "  \033[36m%-15s\033[0m %s\n", $$1, $$2}'

install:  ## Install package for production use
	uv sync --no-dev

dev-install:  ## Install package with development dependencies
	uv sync --extra test

test:  ## Run all tests with coverage
	uv run pytest

test-fast:  ## Run tests without coverage (faster)
	uv run pytest --no-cov

test-cov:  ## Run tests and open coverage report
	uv run pytest
	@echo "Opening coverage report..."
	@python -c "import webbrowser; webbrowser.open('htmlcov/index.html')"

test-unit:  ## Run only unit tests
	uv run pytest -m "not integration"

test-integration:  ## Run only integration tests  
	uv run pytest -m integration

lint:  ## Check code style and formatting
	uv run ruff check bindsweeper/
	uv run ruff check tests/

format:  ## Format code automatically
	uv run ruff format bindsweeper/
	uv run ruff format tests/

clean:  ## Clean up generated files
	rm -rf htmlcov/
	rm -rf .pytest_cache/
	rm -rf .coverage
	rm -rf **/__pycache__/
	rm -rf *.egg-info/
	rm -rf build/
	rm -rf dist/

check-all: lint test  ## Run all checks (lint + test)

.DEFAULT_GOAL := help