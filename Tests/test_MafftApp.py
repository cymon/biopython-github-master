"""
Unittests for Bio.Align.Applications interface for MAFFT

This code is part of the Biopython distribution and governed by its
license.  Please see the LICENSE file that should have been included
as part of this package.
"""

import sys
import os
import unittest
from Bio import Application
from Bio import MissingExternalDependencyError
from Bio.Align.Applications import MafftCommandline

app_name = "mafft"
if sys.platform=="win32" :
    try :
        path = os.environ["MAFFT_ROOT"]
    except KeyError :
        raise MissingExternalDependencyError(\
            "Alignment application MAFFT not found.")
    if os.path.isdir(path) :
        for name in exes_wanted :
            if os.path.isfile(os.path.join(path, app_name+".exe")) :
                exes[name] = os.path.join(path, app_name+".exe")
    del path, name
else :
    import commands
    if "not found" in commands.getoutput("%s -help" % app_name):
        raise MissingExternalDependencyError(\
            "Alignment application MAFFT not found.")

class MafftApplication(unittest.TestCase):

    def setUp(self):
        self.infile1  = "Fasta/f002"

    def tearDown(self):
        if os.path.isfile("Fasta/f002.tree"):
            os.remove("Fasta/f002.tree")

    def test_Mafft_simple(self):
        """Simple round-trip through app with infile.
        Result passed to stdout.
        """
        cmdline = MafftCommandline()
        cmdline.set_parameter("input", self.infile1)
        stdin, stdout, stderr = Application.generic_run(cmdline)
        stderr_string = stderr.read()
        
        self.assert_(stdin.return_code == 0)
        self.assert_(stdout.read().startswith(">gi|1348912|gb|G26680|G26680"))
        self.assert_("STEP     2 / 2 d" in stderr_string)
        self.assert_("$#=0" not in stderr_string)
        self.assert_(str(stdin._cl) == "mafft Fasta/f002 ")

    def test_Mafft_with_options(self):
        """Simple round-trip through app with infile and options.
        Result passed to stdout.
        """
        cmdline = MafftCommandline()
        cmdline.set_parameter("input", self.infile1)
        cmdline.set_parameter("maxiterate", 100)
        cmdline.set_parameter("--localpair")
        stdin, stdout, stderr = Application.generic_run(cmdline)
        
        self.assert_(stdin.return_code == 0)
        self.assert_(stdout.read().startswith(">gi|1348912|gb|G26680|G26680"))
        self.assert_("$#=0" not in stderr.read())
        self.assert_(str(stdin._cl) == "mafft --localpair --maxiterate 100 Fasta/f002 ")

    def test_Mafft_with_Clustalw_output(self):
        """Simple round-trip through app with clustal output"""
        cmdline = MafftCommandline()
        cmdline.set_parameter("input", self.infile1)
        cmdline.set_parameter("--clustalout")
        stdin, stdout, stderr = Application.generic_run(cmdline)
        
        self.assert_(stdin.return_code == 0)
        self.assert_(stdout.read().startswith("CLUSTAL format alignment by MAFFT"))
        self.assert_("$#=0" not in stderr.read())
        self.assert_(str(stdin._cl) == "mafft --clustalout Fasta/f002 ")

    def test_Mafft_with_complex_command_line(self):
        """Round-trip with complex command line."""
        cmdline = MafftCommandline()
        cmdline.set_parameter("input", self.infile1)
        cmdline.set_parameter("--localpair")
        cmdline.set_parameter("--weighti", 4.2)
        cmdline.set_parameter("retree", 5)
        cmdline.set_parameter("maxiterate", 200)
        cmdline.set_parameter("--nofft")
        cmdline.set_parameter("op", 2.04)
        cmdline.set_parameter("--ep", 0.51)
        cmdline.set_parameter("--lop", 0.233)
        cmdline.set_parameter("lep", 0.2)
        cmdline.set_parameter("--reorder")
        cmdline.set_parameter("--treeout")
        cmdline.set_parameter("nuc")

        stdin, stdout, stderr = Application.generic_run(cmdline)

        self.assert_(stdin.return_code == 0)
        self.assert_(stdout.read().startswith(">gi|1348912|gb|G26680|G26680"))
        self.assert_("$#=0" not in stderr.read())
        self.assert_(str(stdin._cl) == "mafft --localpair --weighti 4.2 --retree 5 --maxiterate 200 --nofft --op 2.04 --ep 0.51 --lop 0.233 --lep 0.2 --reorder --treeout --nuc Fasta/f002 ")

if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity = 2)
    unittest.main(testRunner=runner)
