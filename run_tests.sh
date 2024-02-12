set -e
rm -rf docs/generated/
black .
isort .
ruff .
coverage run -m unittest
sphinx-build docs -W -b linkcheck -d _build/doctrees _build/html
