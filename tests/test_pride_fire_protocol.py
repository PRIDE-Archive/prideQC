"""FIRE is opt-in, forwarded to pridepy, never auto-enabled."""

import unittest

from prideqc.cli import parser
from prideqc.pride import DownloadOptions


class FireProtocolTests(unittest.TestCase):
    def test_fire_accepted_in_prideqc(self):
        self.assertEqual(DownloadOptions(protocol="fire").protocol, "fire")
        args = parser().parse_args(
            [
                "analyze",
                "--accession",
                "PXD041271",
                "--file",
                "a.raw",
                "--protocol",
                "fire",
                "-o",
                "out",
            ]
        )
        self.assertEqual(args.protocol, "fire")

    def test_ftp_remains_default(self):
        self.assertEqual(DownloadOptions().protocol, "ftp")


if __name__ == "__main__":
    unittest.main()
