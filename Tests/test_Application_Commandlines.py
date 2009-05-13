# Copyright 2009 by Cymon J. Cox.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
"""
Unittests for Bio.Align.Applications command line interfaces:
    _Clustalw.py
    _Dialign.py
    _Mafft.py
    _Muscle.py
    _Prank.py
    _Probcons.py
    _TCoffee.py

These tests do not require the applications to be installed.
"""

import unittest
from Bio.Application import AbstractCommandline
from Bio.Align.Applications import ClustalwCommandline
from Bio.Align.Applications import DialignCommandline
from Bio.Align.Applications import MafftCommandline
from Bio.Align.Applications import MuscleCommandline
from Bio.Align.Applications import PrankCommandline
from Bio.Align.Applications import ProbconsCommandline
from Bio.Align.Applications import TCoffeeCommandline

class TestApplicationCommandlines(unittest.TestCase):

    def test_Clustalw(self):
        """Clustalw command line"""
        cmdline = ClustalwCommandline("clustalw2",
                                      infile="myinfile.fa",
                                      align=True,
                                      tree=True,
                                      type="dna",
                                      output="pir",
                                      outorder="aligned")
        self.assertEquals(str(cmdline), "clustalw2 -infile=myinfile.fa "
                          "-align -tree -type=dna -output=pir -outorder=aligned ")
        cmdline.align = False
        cmdline.output = "nexus"
        self.assertEquals(str(cmdline), "clustalw2 -infile=myinfile.fa -tree "
                          "-type=dna -output=nexus -outorder=aligned ")
        self.assertRaises(ValueError, cmdline.set_parameter, "output", "msf")
        self.assertRaises(TypeError, cmdline.maxseqlen, 0.99)

    def test_Dialign(self):
        """Dialign command line"""
        cmdline = DialignCommandline(input="myinfile.fa", afc=True, lmax=4, nt=True)
        self.assertEquals(str(cmdline), "dialign2-2 -afc -lmax 4 -nt "
                          "myinfile.fa ")
        cmdline.afc = False
        cmdline.smin = True
        cmdline.o = True
        cmdline.set_parameter("lmax", 3)
        self.assertEquals(str(cmdline), "dialign2-2 -lmax 3 -nt -o -smin "
                          "myinfile.fa ")
        self.assertRaises(ValueError, cmdline.set_parameter, "stars", 100)
        self.assertRaises(TypeError, cmdline.thr, 0.1)

    def test_MafftCommandline(self):
        """Mafft command line"""
        cmdline = MafftCommandline(input="myinfile.fa",
                                   amino=True,
                                   quiet=False,
                                   LEXP=0.02,
                                   maxiterate=1000)
        self.assertEquals(str(cmdline), "mafft --maxiterate 1000 --LEXP "
                          "0.02 --amino myinfile.fa ")
        cmdline.amino = False
        cmdline.quiet = True
        cmdline.lep = 0.5
        cmdline.set_parameter("op", 0.3)
        self.assertEquals(str(cmdline), "mafft --maxiterate 1000 --op 0.3 "
                          "--lep 0.5 --LEXP 0.02 --quiet myinfile.fa ")
        self.assertRaises(ValueError, cmdline.set_parameter, "retree", 0.2)
        self.assertRaises(TypeError, cmdline.bl, "000")

    def test_MuscleCommandline(self):
        """Muscle command line"""
        cmdline = MuscleCommandline("muscle3.7",
                                    input="myinfile.fa",
                                    out="myoutfile.fa",
                                    diags=True,
                                    anchorspacing=4,
                                    diaglength=10,
                                    minsmoothscore=0.3)
        self.assertEquals(str(cmdline), "muscle3.7 -in myinfile.fa -out "
                          "myoutfile.fa -diags -anchorspacing 4 -diaglength "
                          "10 -minsmoothscore 0.3 ")
        cmdline.out = "otherfile.fa"
        cmdline.diags = False
        cmdline.set_parameter("sueff", 0.22)
        self.assertEquals(str(cmdline), "muscle3.7 -in myinfile.fa -out "
                          "otherfile.fa -anchorspacing 4 -diaglength 10 "
                          "-minsmoothscore 0.3 -sueff 0.22 ")
        self.assertRaises(ValueError, cmdline.set_parameter, "objscore", "nnn")
        self.assertRaises(TypeError, cmdline.minbestcolscore, "100")

    def test_PrankCommandline(self):
        """Prank command line"""
        cmdline = PrankCommandline(d="myinfile.fa",
                                   t="mytreefile.newick",
                                   m="WAG",
                                   o="myoutput.msf",
                                   f=15,
                                   noxml=True,
                                   maxbranches=0.5)
        self.assertEquals(str(cmdline), "prank -d=myinfile.fa "
                          "-t=mytreefile.newick -m=WAG -o=myoutput.msf "
                          "-f=15 -noxml -maxbranches=0.5 ")
        cmdline.f = 16
        cmdline.noxml=False
        cmdline.skipins = True
        self.assertEquals(str(cmdline), "prank -d=myinfile.fa -t=mytreefile.newick "
                          "-m=WAG -o=myoutput.msf -f=16 -skipins -maxbranches=0.5 ")
        self.assertRaises(ValueError, cmdline.set_parameter, "pwdist", "1")
        self.assertRaises(TypeError, cmdline.dnafreqs, 0.125)

    def test_ProbconsCommandline(self):
        """Probcons command line"""
        cmdline = ProbconsCommandline("ProbCons",
                                      clustalw=True,
                                      consistency=4,
                                      ir=666,
                                      input="myinfile.fa",
                                      viterbi=True,
                                      annot="myannot.text",
                                      emissions=True)
        self.assertEquals(str(cmdline), "ProbCons -clustalw -c 4 "
                          "-ir 666 -viterbi -annot myannot.text -e myinfile.fa ")
        cmdline.consistency = 2
        cmdline.viterbi = False
        cmdline.set_parameter("pre-training", 20)
        self.assertEquals(str(cmdline), "ProbCons -clustalw -c 2 -ir 666 -pre 20 "
                          "-annot myannot.text -e myinfile.fa ")
        self.assertRaises(ValueError, cmdline.set_parameter, "pre", "22")
        self.assertRaises(TypeError, cmdline.ir, 0.125)

    def test_TCoffeeCommandline(self):
        """TCoffee command line"""
        cmdline = TCoffeeCommandline("tcoffee",
                                     output="pir_aln",
                                     infile="myinfile.fa",
                                     outfile="myoutfile.pir",
                                     type="dna",
                                     quiet=True,
                                     gapopen=-5)
        self.assertEquals(str(cmdline), "tcoffee -output pir_aln -infile "
                          "myinfile.fa -outfile myoutfile.pir -type dna "
                          "-gapopen -5 -quiet ")
        cmdline.quiet = False
        cmdline.set_parameter("gapopen", -2)
        cmdline.set_parameter("outorder", "input")
        self.assertEquals(str(cmdline), "tcoffee -output pir_aln -infile "
                          "myinfile.fa -outfile myoutfile.pir -type dna "
                          "-outorder input -gapopen -2 ")
        self.assertRaises(ValueError, cmdline.set_parameter,
                            "type", "nucleotide")
        self.assertRaises(TypeError, cmdline.gapopen, 0.125)

if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity = 2)
    unittest.main(testRunner=runner)
