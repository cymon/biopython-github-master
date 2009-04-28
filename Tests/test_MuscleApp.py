"""
Unittests for Bio.Align.Applications interface for MUSCLE

This code is part of the Biopython distribution and governed by its
license.  Please see the LICENSE file that should have been included
as part of this package.
"""

import sys
import os
import unittest
from Bio import Application
from Bio import MissingExternalDependencyError
from Bio.Align.Applications import MuscleCommandline

app_name = "muscle"
if sys.platform=="win32" :
    try :
        path = os.environ["MUSCLE_ROOT"]
    except KeyError :
        raise MissingExternalDependencyError(\
            "Alignment application MUSCLE not found.")
    if os.path.isdir(path) :
        for name in exes_wanted :
            if os.path.isfile(os.path.join(path, app_name+".exe")) :
                exes[name] = os.path.join(path, app_name+".exe")
    del path, name
else :
    import commands
    if "not found" in commands.getoutput("%s -help" % app_name):
        raise MissingExternalDependencyError(\
            "Alignment application MUSCLE not found.")

class MuscleApplication(unittest.TestCase):
    
    def setUp(self):
        self.infile1  = "Fasta/f002"
        self.infile2  = "Fasta/fa01"
        self.infile3  = "Fasta/f001"

        self.outfile1 = "Fasta/temp_align_out1.fa"
        self.outfile2 = "Fasta/temp_align_out2.fa"
        self.outfile3 = "Fasta/temp_align_out3.fa"
        self.outfile4 = "Fasta/temp_align_out4.fa"

    def tearDown(self):
        if os.path.isfile(self.outfile1):
            os.remove(self.outfile1)
        if os.path.isfile(self.outfile2):
            os.remove(self.outfile2)
        if os.path.isfile(self.outfile3):
            os.remove(self.outfile3)
        if os.path.isfile(self.outfile4):
            os.remove(self.outfile4)

    def test_Muscle_simple(self):
        """Simple round-trip through app just infile and outfile
        """
        cmdline = MuscleCommandline()
        cmdline.set_parameter("in", self.infile1)
        cmdline.set_parameter("out", self.outfile1)
        stdin, stdout, stderr = Application.generic_run(cmdline)

        self.assert_(stdin.return_code == 0)
        self.assert_(stdout.read() == "")
        self.assert_("*** ERROR ***" not in stderr.read())
        self.assert_(str(stdin._cl) == "muscle -in Fasta/f002 -out " + \
                     "Fasta/temp_align_out1.fa ")

    def test_Muscle_with_options(self):
        """Round-trip through app with a switch and valued option
        """
        cmdline = MuscleCommandline()
        cmdline.set_parameter("in", self.infile1)
        cmdline.set_parameter("out", self.outfile2)
        cmdline.set_parameter("objscore", "sp")
        cmdline.set_parameter("noanchors")
        stdin, stdout, stderr = Application.generic_run(cmdline)

        self.assert_(stdin.return_code == 0)
        self.assert_(stdout.read() == "")
        self.assert_("*** ERROR ***" not in stderr.read())
        self.assert_(str(stdin._cl) == "muscle -in Fasta/f002 -out " + \
                "Fasta/temp_align_out2.fa -objscore sp -noanchors ")

    def test_Muscle_profile_simple(self):
        """Simple round-trip through app doing a profile alignment
        """
        cmdline = MuscleCommandline()
        cmdline.set_parameter("out", self.outfile3)
        cmdline.set_parameter("profile")
        cmdline.set_parameter("in1", self.infile2)
        cmdline.set_parameter("in2", self.infile3)
        stdin, stdout, stderr = Application.generic_run(cmdline)
        
        self.assert_(stdin.return_code == 0)
        self.assert_(stdout.read() == "")
        self.assert_("*** ERROR ***" not in stderr.read())
        self.assert_(str(stdin._cl) == "muscle -out Fasta/temp_align_out3.fa " + \
                     "-profile -in1 Fasta/fa01 -in2 Fasta/f001 ")

    def test_Muscle_profile_with_options(self):
        """Profile alignment, and switch and valued options.
        """
        cmdline = MuscleCommandline()
        cmdline.set_parameter("out", self.outfile4)
        cmdline.set_parameter("profile")
        cmdline.set_parameter("in1", self.infile2)
        cmdline.set_parameter("in2", self.infile3)
        cmdline.set_parameter("cluster1", "neighborjoining")
        cmdline.set_parameter("stable")
        stdin, stdout, stderr = Application.generic_run(cmdline)

        self.assert_(stdin.return_code == 0)
        self.assert_(stdout.read() == "")
        self.assert_("*** ERROR ***" not in stderr.read())
        self.assert_(str(stdin._cl) == "muscle -out Fasta/temp_align_out4.fa " + \
                     "-profile -in1 Fasta/fa01 -in2 Fasta/f001 -cluster1 " + \
                     "neighborjoining -stable ")

if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity = 2)
    unittest.main(testRunner=runner)
