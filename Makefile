.PHONY: test test-unit test-integration test-all test-cov test-fast clean-test

# Run all tests
test:
	pytest tests/ -v

# Run only unit tests
test-unit:
	pytest tests/ -v -m "unit or not integration"

# Run only integration tests  
test-integration:
	pytest tests/ -v -m "integration"

# Run all tests with coverage
test-cov:
	pytest tests/ -v --cov=src --cov-report=html --cov-report=term-missing

# Run tests in parallel (fast)
test-fast:
	pytest tests/ -v -n auto

# Run specific test file
test-file:
	pytest tests/test_splice_extractor.py -v

# Run tests with specific marker
test-slow:
	pytest tests/ -v -m "slow"

# Skip slow tests
test-quick:
	pytest tests/ -v -m "not slow"

# Clean test artifacts
clean-test:
	rm -rf .pytest_cache/
	rm -rf htmlcov/
	rm -rf .coverage
	find . -type d -name __pycache__ -exec rm -rf {} +
	find . -type f -name "*.pyc" -delete

# Code quality checks
lint:
	flake8 src/ tests/
	black --check src/ tests/
	isort --check-only src/ tests/
	mypy src/

# Format code
format:
	black src/ tests/
	isort src/ tests/

# Run all quality checks and tests
check: lint test-cov