# Copyright 2009 by Cymon J. Cox.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""
Bio.Application command line for the multiple alignment programme DIALIGN2-2

http://bibiserv.techfak.uni-bielefeld.de/dialign/welcome.html

Citations:

B. Morgenstern (2004). DIALIGN: Multiple DNA and Protein Sequence Alignment
at BiBiServ. Nucleic Acids Research 32, W33-W36.

Last checked against version: 2.2
"""
import os
import types
from Bio import Application
from Bio.Application import _Option
from Bio.Application import _Argument

class DialignCommandline(Application.AbstractCommandline):
    
    def __init__(self, cmd = "dialign2-2"):

        Application.AbstractCommandline.__init__(self)
        self.program_name = cmd
        self.parameters = \
            [
            _Option(["-afc", "afc"], ["input"],
                    lambda x: 0, #Does not take a value
                    0,
                    "Creates additional output file '*.afc' " + \
                    "containing data of all fragments considered " + \
                    "for alignment WARNING: this file can be HUGE !",
                    0),

            _Option(["-afc_v", "afc_v"], ["input"],
                    lambda x: 0, #Does not take a value
                    0,
                    "Like '-afc' but verbose: fragments are explicitly " + \
                    "printed. WARNING: this file can be EVEN BIGGER !",
                    0),

            _Option(["-anc", "anc"], ["input"],
                    lambda x: 0, #Does not take a value
                    0,
                    "Anchored alignment. Requires a file <seq_file>.anc " + \
                    "containing anchor points.",
                    0),

            _Option(["-cs", "cs"], ["input"],
                    lambda x: 0, #Does not take a value
                    0,
                    "If segments are translated, not only the `Watson " + \
                    "strand' but also the `Crick strand' is looked at.",
                    0),

            _Option(["-cw", "cw"], ["input"],
                    lambda x: 0, #Does not take a value
                    0,
                    "Additional output file in CLUSTAL W format.",
                    0),

            _Option(["-ds", "ds"], ["input"],
                    lambda x: 0, #Does not take a value
                    0,
                    "`dna alignment speed up' - non-translated nucleic acid " + \
                    "fragments are taken into account only if they start " + \
                    "with at least two matches. Speeds up DNA alignment at " + \
                    "the expense of sensitivity.",
                    0),

            _Option(["-fa", "fa"], ["input"],
                    lambda x: 0, #Does not take a value
                    0,
                    "Additional output file in FASTA format.",
                    0),

            _Option(["-ff", "ff"], ["input"],
                    lambda x: 0, #Does not take a value
                    0,
                    "Creates file *.frg containing information about all " + \
                    "fragments that are part of the respective optimal " + \
                    "pairwise alignmnets plus information about " + \
                    "consistency in the multiple alignment",
                    0),

            _Option(["-fn", "fn"], ["input"],
                    None,
                    0,
                    "Output files are named <out_file>.<extension>.",
                    0),

            _Option(["-fop", "fop"], ["input"],
                    lambda x: 0, #Does not take a value
                    0,
                    "Creates file *.fop containing coordinates of all " + \
                    "fragments that are part of the respective pairwise alignments.",
                    0),

            _Option(["-fsm", "fsm"], ["input"],
                    lambda x: 0, #Does not take a value
                    0,
                    "Creates file *.fsm containing coordinates of all " + \
                    "fragments that are part of the final alignment",
                    0),

            _Option(["-iw", "iw"], ["input"],
                    lambda x: 0, #Does not take a value
                    0,
                    "Overlap weights switched off (by default, overlap " + \
                    "weights are used if up to 35 sequences are aligned). " + \
                    "This option speeds up the alignment but may lead " + \
                    "to reduced alignment quality.",
                    0),

            _Option(["-lgs", "lgs"], ["input"],
                    lambda x: 0, #Does not take a value
                    0,
                    "`long genomic sequences' - combines the following " + \
                    "options: -ma, -thr 2, -lmax 30, -smin 8, -nta, -ff, " + \
                    "-fop, -ff, -cs, -ds, -pst ",
                    0),

            _Option(["-lgs_t", "lgs_t"], ["input"],
                    lambda x: 0, #Does not take a value
                    0,
                    "Like '-lgs' but with all segment pairs assessed " + \
                    "at the peptide level (rather than 'mixed alignments' " + \
                    "as with the '-lgs' option). Therefore faster than " + \
                    "-lgs but not very sensitive for non-coding regions.",
                    0),

            _Option(["-lmax", "lmax"], ["input"],
                    lambda x: isinstance(x, types.IntType),
                    0,
                    "Maximum fragment length = x  (default: x = 40 or " + \
                    "x = 120 for `translated' fragments). Shorter x " + \
                    "speeds up the program but may affect alignment quality.",
                    0),

            _Option(["-lo", "lo"], ["input"],
                    lambda x: 0, #Does not take a value
                    0,
                    "(Long Output) Additional file *.log with information " + \
                    "about fragments selected for pairwise alignment and " + \
                    "about consistency in multi-alignment proceedure.",
                    0),

            _Option(["-ma", "ma"], ["input"],
                    lambda x: 0, #Does not take a value
                    0,
                    "`mixed alignments' consisting of P-fragments and " + \
                    "N-fragments if nucleic acid sequences are aligned.",
                    0),

            _Option(["-mask", "mask"], ["input"],
                    lambda x: 0, #Does not take a value
                    0,
                    "Residues not belonging to selected fragments are " + \
                    "replaced by `*' characters in output alignment " + \
                    "(rather than being printed in lower-case characters)",
                    0),

            _Option(["-mat", "mat"], ["input"],
                    lambda x: 0, #Does not take a value
                    0,
                    "Creates file *mat with substitution counts derived " + \
                    "from the fragments that have been selected for alignment.",
                    0),

            _Option(["-mat_thr", "mat_thr"], ["input"],
                    lambda x: 0, #Does not take a value
                    0,
                    "Like '-mat' but only fragments with weight score " + \
                    "> t are considered",
                    0),

            _Option(["-max_link", "max_link"], ["input"],
                    lambda x: 0, #Does not take a value
                    0,
                    "'maximum linkage' clustering used to construct " + \
                    "sequence tree (instead of UPGMA).",
                    0),

            _Option(["-min_link", "min_link"], ["input"],
                    lambda x: 0, #Does not take a value
                    0,
                    "'minimum linkage' clustering used.",
                    0),

            _Option(["-mot", "mot"], ["input"],
                    None, 
                    0,
                    "'motif' option.",
                    0),

            _Option(["-msf", "msf"], ["input"],
                    lambda x: 0, #Does not take a value
                    0,
                    "Separate output file in MSF format.",
                    0),

            _Option(["-n", "n"], ["input"],
                    lambda x: 0, #Does not take a value
                    0,
                    "Input sequences are nucleic acid sequences. " + \
                    "No translation of fragments.",
                    0),

            _Option(["-nt", "nt"], ["input"],
                    lambda x: 0, #Does not take a value
                    0,
                    "Input sequences are nucleic acid sequences and " + \
                    "`nucleic acid segments' are translated to `peptide " + \
                    "segments'.",
                    0),

            _Option(["-nta", "nta"], ["input"],
                    lambda x: 0, #Does not take a value
                    0,
                    "`no textual alignment' - textual alignment suppressed. " + \
                    "This option makes sense if other output files are of " + \
                    "intrest -- e.g. the fragment files created with -ff, " + \
                    "-fop, -fsm or -lo.",
                    0),

            _Option(["-o", "o"], ["input"],
                    lambda x: 0, #Does not take a value
                    0,
                    "Fast version, resulting alignments may be slightly " + \
                    "different.",
                    0),

            _Option(["-ow", "ow"], ["input"],
                    lambda x: 0, #Does not take a value
                    0,
                    "Overlap weights enforced (By default, overlap weights " + \
                    "are used only if up to 35 sequences are aligned since " + \
                    "calculating overlap weights is time consuming).",
                    0),

            _Option(["-pst", "pst"], ["input"],
                    lambda x: 0, #Does not take a value
                    0,
                    "'print status'. Creates and updates a file *.sta with " + \
                    "information about the current status of the program " + \
                    "run.  This option is recommended if large data sets " + \
                    "are aligned since it allows the user to estimate the " + \
                    "remaining running time.",
                    0),

            _Option(["-smin", "smin"], ["input"],
                    lambda x: 0, #Does not take a value
                    0,
                    "Minimum similarity value for first residue pair " + \
                    "(or codon pair) in fragments. Speeds up protein " + \
                    "alignment or alignment of translated DNA fragments " + \
                    "at the expense of sensitivity.",
                    0),

            _Option(["-stars", "stars"], ["input"],
                    lambda x: x in range(0,10), #Does not take a value
                    0,
                    "Maximum number of `*' characters indicating degree " + \
                    "of local similarity among sequences. By default, no " + \
                    "stars are used but numbers between 0 and 9, instead.",
                    0),

            _Option(["-stdo", "stdo"], ["input"],
                    lambda x: 0, #Does not take a value
                    0,
                    "Results written to standard output.",
                    0),

            _Option(["-ta", "ta"], ["input"],
                    lambda x: 0, #Does not take a value
                    0,
                    "Standard textual alignment printed (overrides " + \
                    "suppression of textual alignments in special " + \
                    "options, e.g. -lgs)",
                    0),

            _Option(["-thr", "thr"], ["input"],
                    lambda x: isinstance(x, types.IntType),
                    0,
                    "Threshold T = x.",
                    0),

            _Option(["-xfr", "xfr"], ["input"],
                    lambda x: 0, #Does not take a value
                    0,
                    "'exclude fragments' - list of fragments can be " + \
                    "specified that are NOT considered for pairwise alignment",
                    0),

            _Argument(["input"], ["input", "file"], os.path.exists, 1,
                      "Input file name. Must be FASTA format")

            ]

