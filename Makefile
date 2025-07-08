# eccToolkit Makefile for development and maintenance tasks

.PHONY: help install install-dev format lint type-check test test-cov clean docs setup-pre-commit style-check style-fix all-checks

# Default target
help:
	@echo "eccToolkit Development Commands:"
	@echo ""
	@echo "Setup Commands:"
	@echo "  install         Install package in development mode"
	@echo "  install-dev     Install package with development dependencies"
	@echo "  setup-pre-commit Install and setup pre-commit hooks"
	@echo ""
	@echo "Code Quality Commands:"
	@echo "  format          Format code with black and isort"
	@echo "  lint            Run flake8 linting"
	@echo "  type-check      Run mypy type checking"
	@echo "  style-check     Run all style checks (format, lint, type-check)"
	@echo "  style-fix       Auto-fix all style issues"
	@echo ""
	@echo "Testing Commands:"
	@echo "  test            Run pytest tests"
	@echo "  test-cov        Run tests with coverage report"
	@echo ""
	@echo "Maintenance Commands:"
	@echo "  clean           Remove build artifacts and cache files"
	@echo "  docs            Build documentation"
	@echo "  all-checks      Run all quality checks"

# Installation commands
install:
	pip install -e .

install-dev:
	pip install -e ".[dev,plotting]"

setup-pre-commit: install-dev
	pre-commit install
	pre-commit install --hook-type commit-msg
	@echo "Pre-commit hooks installed successfully!"

# Code formatting
format:
	@echo "Running isort..."
	isort code/ tests/ --profile black
	@echo "Running black..."
	black code/ tests/
	@echo "Code formatting complete!"

# Linting
lint:
	@echo "Running flake8..."
	flake8 code/ tests/
	@echo "Linting complete!"

# Type checking
type-check:
	@echo "Running mypy..."
	mypy code/ --ignore-missing-imports
	@echo "Type checking complete!"

# Combined style checking
style-check: lint type-check
	@echo "Running black check..."
	black --check code/ tests/
	@echo "Running isort check..."
	isort --check-only code/ tests/ --profile black
	@echo "All style checks passed!"

# Auto-fix style issues
style-fix: format
	@echo "Running autoflake..."
	autoflake --in-place --remove-all-unused-imports --remove-unused-variables --recursive code/ tests/
	@echo "Style fixes applied!"

# Testing
test:
	@echo "Running pytest..."
	pytest tests/ -v
	@echo "Tests complete!"

test-cov:
	@echo "Running pytest with coverage..."
	pytest tests/ -v --cov=code --cov-report=html --cov-report=term-missing
	@echo "Coverage report generated in htmlcov/"

# Documentation
docs:
	@echo "Building documentation..."
	# Add documentation build commands here when available
	@echo "Documentation build complete!"

# Cleanup
clean:
	@echo "Cleaning build artifacts..."
	rm -rf build/
	rm -rf dist/
	rm -rf *.egg-info/
	rm -rf htmlcov/
	rm -rf .coverage
	rm -rf .pytest_cache/
	rm -rf .mypy_cache/
	find . -type d -name __pycache__ -exec rm -rf {} +
	find . -type f -name "*.pyc" -delete
	find . -type f -name "*.pyo" -delete
	@echo "Cleanup complete!"

# Run all quality checks
all-checks: style-check test
	@echo "All quality checks passed!"

# Development workflow commands
dev-setup: install-dev setup-pre-commit
	@echo "Development environment setup complete!"

# Quick development check
quick-check:
	@echo "Running quick development checks..."
	black --check code/ --quiet
	isort --check-only code/ --profile black --quiet
	flake8 code/ --quiet
	@echo "Quick checks passed!"

# CI/CD simulation
ci-check: style-check test-cov
	@echo "CI/CD checks simulation complete!"

# Format specific directories or files
format-code:
	black code/
	isort code/ --profile black

format-tests:
	black tests/
	isort tests/ --profile black

# Security checking
security-check:
	@echo "Running bandit security check..."
	bandit -r code/ -f json -o bandit-report.json || true
	@echo "Security check complete! Check bandit-report.json for results."

# Dependency management
update-deps:
	@echo "Updating dependencies..."
	pip install --upgrade pip
	pip install --upgrade -e ".[dev,plotting]"
	@echo "Dependencies updated!"

# Pre-commit management
pre-commit-run:
	@echo "Running pre-commit on all files..."
	pre-commit run --all-files

pre-commit-update:
	@echo "Updating pre-commit hooks..."
	pre-commit autoupdate

# Performance profiling helpers
profile-script:
	@echo "Usage: make profile-script SCRIPT=code/script_name.py ARGS='--arg1 value1'"
	@echo "Example: make profile-script SCRIPT=code/degs_enrichment_analysis.py ARGS='-d test.txt -e test.csv -o output.csv'"
	python -m cProfile -o profile_output.prof $(SCRIPT) $(ARGS)
	@echo "Profile saved to profile_output.prof"

# Code complexity analysis
complexity:
	@echo "Analyzing code complexity..."
	flake8 code/ --select=C901 --max-complexity=10
	@echo "Complexity analysis complete!"

# Line count statistics
stats:
	@echo "Code statistics:"
	@echo "Python files:"
	find code/ -name "*.py" | wc -l
	@echo "Total lines of Python code:"
	find code/ -name "*.py" -exec wc -l {} + | tail -1
	@echo "Test files:"
	find tests/ -name "*.py" 2>/dev/null | wc -l || echo "0"

# Git hooks
install-git-hooks: setup-pre-commit
	@echo "Git hooks installed!"

# Environment validation
validate-env:
	@echo "Validating development environment..."
	python --version
	pip --version
	black --version
	isort --version
	flake8 --version
	mypy --version
	pytest --version
	@echo "Environment validation complete!"

# Documentation generation for scripts
generate-script-docs:
	@echo "Generating documentation for all scripts..."
	python -c "
import os
import sys
sys.path.append('code')

for script in os.listdir('code'):
    if script.endswith('.py') and not script.startswith('__'):
        print(f'=== {script} ===')
        try:
            with open(f'code/{script}', 'r') as f:
                lines = f.readlines()
                # Find docstring
                in_docstring = False
                for line in lines[:20]:  # Check first 20 lines
                    if '\"\"\"' in line:
                        if not in_docstring:
                            in_docstring = True
                            print(line.strip())
                        else:
                            print(line.strip())
                            break
                    elif in_docstring:
                        print(line.strip())
        except:
            pass
        print()
"

# Backup current code
backup:
	@echo "Creating backup of current code..."
	cp -r code/ "code_backup_$(shell date +%Y%m%d_%H%M%S)/"
	@echo "Backup created!"

# Show current style violations
show-violations:
	@echo "Current style violations:"
	@echo "=== Black formatting issues ==="
	black --check code/ --diff || true
	@echo ""
	@echo "=== Import sorting issues ==="
	isort --check-only code/ --diff --profile black || true
	@echo ""
	@echo "=== Flake8 issues ==="
	flake8 code/ || true
	@echo ""
	@echo "=== MyPy issues ==="
	mypy code/ --ignore-missing-imports || true