set -e
rm -rf docs/generated/
black .
isort .
ruff .
coverage run --source=./tests -m unittest discover -s tests/
sphinx-build docs -W -b linkcheck -d _build/doctrees _build/html
