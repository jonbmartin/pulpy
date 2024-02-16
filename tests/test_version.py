import unittest

from pulpy import version

if __name__ == "__main__":
    unittest.main()


class TestVersion(unittest.TestCase):
    def test_version(self):
        assert version.__version__ == "1.8.1"
