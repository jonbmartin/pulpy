import unittest
import importlib.metadata

from pulpy import version

if __name__ == "__main__":
    unittest.main()


class TestVersion(unittest.TestCase):
    def test_version(self):
        assert version.__version__ == importlib.metadata.version("mypackage")
