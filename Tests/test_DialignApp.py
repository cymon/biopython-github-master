"""
Unittests for Bio.Align.Applications interface for DIALIGN2-2

This code is part of the Biopython distribution and governed by its
license.  Please see the LICENSE file that should have been included
as part of this package.
"""

import sys
import os
import unittest
from Bio import Application
from Bio import MissingExternalDependencyError
from Bio.Align.Applications import DialignCommandline

app_name = "dialign2-2"
if sys.platform=="win32" :
    try :
        path = os.environ["DIALIGN2-2_ROOT"]
    except KeyError :
        raise MissingExternalDependencyError(\
            "Alignment application DIALIGN2-2 not found.")
    if os.path.isdir(path) :
        for name in exes_wanted :
            if os.path.isfile(os.path.join(path, app_name+".exe")) :
                exes[name] = os.path.join(path, app_name+".exe")
    del path, name
else :
    import commands
    if "not found" in commands.getoutput("%s -help" % app_name):
        raise MissingExternalDependencyError(\
            "Alignment application DIALIGN2-2 not found.")

class DialignApplication(unittest.TestCase):

    def setUp(self):
        self.infile1 = "Fasta/f002" 
        #Standard output file
        self.outfile1 = "Fasta/f002.ali"
        #MSF output
        self.outfile2 = "Fasta/f002.ms"

    def tearDown(self):
        if os.path.isfile(self.outfile1):
            os.remove(self.outfile1)
        if os.path.isfile(self.outfile2):
            os.remove(self.outfile2)

    def test_Dialign_simple(self):
        """Simple round-trip through app with infile.
        """
        cmdline = DialignCommandline()
        cmdline.set_parameter("input", self.infile1)
        stdin, stdout, stderr = Application.generic_run(cmdline)
        
        self.assert_(stdin.return_code == 0)
        self.assert_(os.path.exists(self.outfile1))
        self.assert_(stdout.read() == "")
        self.assert_(stderr.read() == "")
        self.assert_(str(stdin._cl) == "dialign2-2 Fasta/f002 ")

    def test_Dialign_simple_with_options(self):
        """Simple round-trip through app with infile and options
        """
        cmdline = DialignCommandline()
        cmdline.set_parameter("input", self.infile1)
        cmdline.set_parameter("-max_link")
        cmdline.set_parameter("stars", 4)
        stdin, stdout, stderr = Application.generic_run(cmdline)
        
        self.assert_(stdin.return_code == 0)
        self.assert_(os.path.exists(self.outfile1))
        self.assert_(stdout.read() == "")
        self.assert_(stderr.read() == "")
        self.assert_(str(stdin._cl) == "dialign2-2 -max_link -stars 4 Fasta/f002 ")

    def test_Dialign_simple_with_MSF_output(self):
        """Simple round-trip through app with infile, output MSF
        """
        cmdline = DialignCommandline()
        cmdline.set_parameter("input", self.infile1)
        cmdline.set_parameter("-msf")
        stdin, stdout, stderr = Application.generic_run(cmdline)
        
        self.assert_(stdin.return_code == 0)
        self.assert_(os.path.exists(self.outfile1))
        self.assert_(os.path.exists(self.outfile2))
        self.assert_(stdout.read() == "")
        self.assert_(stderr.read() == "")
        self.assert_(str(stdin._cl) == "dialign2-2 -msf Fasta/f002 ")

    def test_Dialign_complex_command_line(self):
        """Round-trip through app with complex command line."""
        cmdline = DialignCommandline()
        cmdline.set_parameter("input", self.infile1)
        cmdline.set_parameter("-nt")
        cmdline.set_parameter("-thr", 4)
        cmdline.set_parameter("stars", 9)
        cmdline.set_parameter("-ow")
        cmdline.set_parameter("mask")
        cmdline.set_parameter("-cs")

        stdin, stdout, stderr = Application.generic_run(cmdline)
        
        self.assert_(stdin.return_code == 0)
        self.assert_(os.path.exists(self.outfile1))
        self.assert_(stdout.read().startswith(" e_len = 633"))
        self.assert_(stderr.read() == "")
        self.assert_(str(stdin._cl) == "dialign2-2 -cs -mask -nt -ow -stars 9 -thr 4 Fasta/f002 ")

if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity = 2)
    unittest.main(testRunner=runner)
