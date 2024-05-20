set -e
rm -rf docs/generated/
sphinx-build docs -W -b linkcheck -d _build/doctrees _build/html
