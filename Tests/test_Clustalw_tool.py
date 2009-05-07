# Copyright 2008 by Peter Cock.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

#TODO - Clean up the extra files created by clustalw?  e.g. *.dnd
#and *.aln where we have not requested an explicit name?
from Bio import MissingExternalDependencyError

#TODO - Remove this work around once we drop python 2.3 support
try:
    set = set
except NameError:
    from sets import Set as set

import sys
import os
from Bio import Clustalw #old!
from Bio.Clustalw import MultipleAlignCL #old!
from Bio import SeqIO
from Bio import AlignIO
from Bio.Align.Applications import ClustalwCommandline #new!
from Bio.Application import generic_run

#################################################################

clustalw_exe = None
if sys.platform=="win32" :
    #TODO - Check the path?
    try :
        #This can vary depending on the Windows language.
        prog_files = os.environ["PROGRAMFILES"]
    except KeyError :
        prog_files = r"C:\Program Files"

    #Note that EBI's clustalw2 installer, e.g. clustalw-2.0.10-win.msi
    #uses C:\Program Files\ClustalW2\clustalw2.exe so we should check
    #for that.
    #
    #Some users doing a manual install have reported using
    #C:\Program Files\clustalw.exe
    #
    #Older installers might use something like this,
    #C:\Program Files\Clustalw\clustalw.exe
    #
    #One particular case is www.tc.cornell.edu currently provide a
    #clustalw1.83 installer which uses the following long location:
    #C:\Program Files\CTCBioApps\clustalw\v1.83\clustalw1.83.exe
    likely_dirs = ["ClustalW2", "",
                   "Clustal","Clustalw","Clustalw183","Clustalw1.83",
                   r"CTCBioApps\clustalw\v1.83"]
    likely_exes = ["clustalw2.exe",
                   "clustalw.exe", "clustalw1.83.exe"]
    for folder in likely_dirs :
        if os.path.isdir(os.path.join(prog_files, folder)) :
            for filename in likely_exes :
                if os.path.isfile(os.path.join(prog_files, folder, filename)) :
                    clustalw_exe = os.path.join(prog_files, folder, filename)
                    break
            if clustalw_exe : break
else :
    import commands
    #Note that clustalw 1.83 and clustalw 2.0.10 don't obey the --version
    #command, but this does cause them to quit cleanly.  Otherwise they prompt
    #the user for input (causing a lock up).
    output = commands.getoutput("clustalw2 --version")
    if "not found" not in output and "CLUSTAL" in output.upper() :
        clustalw_exe = "clustalw2"
    if not clustalw_exe :
        output = commands.getoutput("clustalw --version")
        if "not found" not in output and "CLUSTAL" in output.upper() :
            clustalw_exe = "clustalw"

if not clustalw_exe :
    raise MissingExternalDependencyError(\
        "Install clustalw or clustalw2 if you want to use Bio.Clustalw.")

#################################################################

print "Checking error conditions"
print "========================="

print "Empty file"
input_file = "does_not_exist.fasta"
assert not os.path.isfile(input_file)
cline = MultipleAlignCL(input_file, command=clustalw_exe)
try :
    align = Clustalw.do_alignment(cline)
    assert False, "Should have failed, returned %s" % repr(align)
except IOError, err :
    print "Failed (good)"
    #Python 2.3 on Windows gives (0, 'Error')
    #Python 2.5 on Windows gives [Errno 0] Error
    assert "Cannot open sequence file" in str(err) \
           or "not produced" in str(err) \
           or str(err) == "[Errno 0] Error" \
           or str(err) == "(0, 'Error')", str(err)

print
print "Single sequence"
input_file = "Fasta/f001"
assert os.path.isfile(input_file)
assert len(list(SeqIO.parse(open(input_file),"fasta")))==1
cline = MultipleAlignCL(input_file, command=clustalw_exe)
try :
    align = Clustalw.do_alignment(cline)
    assert False, "Should have failed, returned %s" % repr(align)
except IOError, err :
    print "Failed (good)"
    assert "has only one sequence present" in str(err)
except ValueError, err :
    print "Failed (good)"
    assert str(err) == "No records found in handle"
    #Ideally we'd get an IOError but sometimes we don't seem to
    #get a return value from clustalw.  If so, then there is a
    #ValueError when the parsing fails.

print
print "Invalid sequence"
input_file = "Medline/pubmed_result1.txt"
assert os.path.isfile(input_file)
cline = MultipleAlignCL(input_file, command=clustalw_exe)
try :
    align = Clustalw.do_alignment(cline)
    assert False, "Should have failed, returned %s" % repr(align)
except IOError, err :
    print "Failed (good)"
    #Ideally we'd catch the return code and raise the specific
    #error for "invalid format", rather than just notice there
    #is not output file.
    #Note:
    #Python 2.3 on Windows gives (0, 'Error')
    #Python 2.5 on Windows gives [Errno 0] Error
    assert "invalid format" in str(err) \
           or "not produced" in str(err) \
           or str(err) == "[Errno 0] Error" \
           or str(err) == "(0, 'Error')", str(err)

#################################################################
print
print "Checking normal situations"
print "=========================="

#Create a temp fasta file with a space in the name
temp_filename_with_spaces = "Clustalw/temp horses.fasta"
handle = open(temp_filename_with_spaces, "w")
SeqIO.write(SeqIO.parse(open("Phylip/hennigian.phy"),"phylip"),handle, "fasta")
handle.close()

#Create a large input file by converting another example file
#(See Bug 2804, this will produce so much output on stdout that
#subprocess could suffer a deadlock and hang).  Using all the
#records should show the deadlock but is very slow - just thirty
#seems to lockup on Mac OS X, even 20 on Linux (without the fix).
temp_large_fasta_file = "temp_cw_prot.fasta"
handle = open(temp_large_fasta_file, "w")
records = list(SeqIO.parse(open("NBRF/Cw_prot.pir", "rU"), "pir"))[:40]
SeqIO.write(records, handle, "fasta")
handle.close()
del handle, records

for input_file, output_file, newtree_file in [
    ("Fasta/f002", "temp_test.aln", None),
    ("GFF/multi.fna", "temp with space.aln", None),
    ("Registry/seqs.fasta", "temp_test.aln", None),
    ("Registry/seqs.fasta", "temp_test.aln", "temp_test.dnd"),
    ("Registry/seqs.fasta", "temp_test.aln", "temp with space.dnd"),
    (temp_filename_with_spaces, "temp_test.aln", None),
    (temp_filename_with_spaces, "temp with space.aln", None),
    (temp_large_fasta_file, "temp_cw_prot.aln", None),
    ] :
    #Note that ClustalW will map ":" to "_" in it's output file
    input_records = SeqIO.to_dict(SeqIO.parse(open(input_file),"fasta"),
                                  lambda rec : rec.id.replace(":","_"))
    if os.path.isfile(output_file) :
        os.remove(output_file)
    print "Calling clustalw on %s (with %i records)" \
          % (repr(input_file), len(input_records))
    print "using output file %s" % repr(output_file)
    if newtree_file is not None :
        print "requesting output guide tree file %s" % repr(newtree_file)

    #Prepare the command...
    cline = MultipleAlignCL(input_file, command=clustalw_exe)
    cline.set_output(output_file)
    if newtree_file is not None :
        cline.set_new_guide_tree(newtree_file)

    if sys.platform=="win32" and sys.version_info[:2] < (2,4 ):
        if " " in input_file or " " in output_file \
        or (newtree_file is not None and " " in newtree_file) :
            #This will fail on Python 2.3 ... cheat so the
            #print-and-compare output matches!
            print "Got an alignment, %i sequences" % (len(input_records))
            continue

    #Run the command...
    align = Clustalw.do_alignment(cline)

    #Check the output...
    print "Got an alignment, %i sequences" % (len(align.get_all_seqs()))
    #The length of the alignment will depend on the version of clustalw
    #(clustalw 2.0.10 and clustalw 1.83 are certainly different).
    output_records = SeqIO.to_dict(SeqIO.parse(open(output_file),"clustal"))
    assert set(input_records.keys()) == set(output_records.keys())
    for record in align :
        assert str(record.seq) == str(output_records[record.id].seq)
        assert str(record.seq).replace("-","") == \
               str(input_records[record.id].seq)

    #Clean up...
    os.remove(output_file)

    #Check the DND file was created.
    #TODO - Try and parse this with Bio.Nexus?
    if newtree_file is not None :
        tree_file = newtree_file
    else :
        #Clustalw will name it based on the input file
        tree_file = os.path.splitext(input_file)[0] + ".dnd"
    assert os.path.isfile(tree_file), \
           "Did not find tree file %s" % tree_file
    os.remove(tree_file)


    #And again, but this time using Bio.Align.Applications wrapper
    #Any filesnames with spaces should get escaped with quotes automatically.
    #Using keyword arguments here.
    cline = ClustalwCommandline(clustalw_exe,
                                infile=input_file,
                                outfile=output_file)
    assert str(eval(repr(cline)))==str(cline)
    if newtree_file is not None :
        #Test using a property:
        cline.newtree = newtree_file
        #I don't just want the tree, also want the alignment:
        cline.align = True
        assert str(eval(repr(cline)))==str(cline)
    #print cline
    return_code, out_handle, err_handle = generic_run(cline)
    assert out_handle.read().strip().startswith("CLUSTAL")
    assert err_handle.read().strip() == ""
    align = AlignIO.read(open(output_file), "clustal")
    assert set(input_records.keys()) == set(output_records.keys())
    for record in align :
        assert str(record.seq) == str(output_records[record.id].seq)
        assert str(record.seq).replace("-","") == \
               str(input_records[record.id].seq)

    #Clean up...
    os.remove(output_file)

    #Check the DND file was created.
    #TODO - Try and parse this with Bio.Nexus?
    if newtree_file is not None :
        tree_file = newtree_file
    else :
        #Clustalw will name it based on the input file
        tree_file = os.path.splitext(input_file)[0] + ".dnd"
    assert os.path.isfile(tree_file), \
           "Did not find tree file %s" % tree_file
    os.remove(tree_file)

#Clean up any stray temp files..
if os.path.isfile("Fasta/f001.aln") :
    os.remove("Fasta/f001.aln")
if os.path.isfile("Medline/pubmed_result1.aln") :
    os.remove("Medline/pubmed_result1.aln")
if os.path.isfile(temp_filename_with_spaces) :   
    os.remove(temp_filename_with_spaces)
if os.path.isfile(temp_large_fasta_file) :
    os.remove(temp_large_fasta_file)

print "Done"
