# Copyright 2008-2009 by Peter Cock.  All rights reserved.
# Revisions copyright 2009 by Cymon J. Cox.  All rights reserved.
#
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Bio.SeqIO support for the "phd" file format.

PHD files are output by PHRED and used by PHRAP and CONSED.

You are expected to use this module via the Bio.SeqIO functions.
See also the underlying Bio.Sequencing.Phd module."""

from Bio.SeqRecord import SeqRecord
from Bio.Sequencing import Phd
from Bio.SeqIO.Interfaces import SequentialSequenceWriter
from Bio.SeqIO import QualityIO

#This is a generator function!
def PhdIterator(handle) :
    """Returns SeqRecord objects from a PHD file.

    This uses the Bio.Sequencing.Phd module to do the hard work.
    """
    phd_records = Phd.parse(handle)
    for phd_record in phd_records:
        #Convert the PHY record into a SeqRecord...
        seq_record = SeqRecord(phd_record.seq,
                               id = phd_record.file_name,
                               name = phd_record.file_name)
        #Just re-use the comments dictionary as the SeqRecord's annotations
        seq_record.annotations = phd_record.comments
        seq_record.letter_annotations["phred_quality"] = \
                [int(site[1]) for site in phd_record.sites]
        #Peak locations are optional
        try:
            seq_record.letter_annotations["peak_location"] = \
                [int(site[2]) for site in phd_record.sites]
        except IndexError:
            pass
        yield seq_record 
    #All done

class PhdWriter(SequentialSequenceWriter):
    """Class to write Phd format files"""

    def __init__(self, handle):
        SequentialSequenceWriter.__init__(self, handle)

    def write_record(self, record):
        """Write a single Phd record to the file."""
        assert record.seq, "No sequence present in SeqRecord"
        # This method returns the 'phred_quality' scores or converted 'solexa_quality'
        # scores if present, else raises a value error
        phred_qualities = QualityIO._get_phred_quality(record)
        peak_locations = record.letter_annotations.get("peak_location", None)
        assert len(record.seq) == len(phred_qualities), "Number of " + \
                "phd quality scores does not match length of sequence"
        if peak_locations:
            assert len(record.seq) == len(peak_locations), "Number " + \
                    "of peak location scores does not match length of sequence"
        self.handle.write("BEGIN_SEQUENCE %s\nBEGIN_COMMENT\n" % record.name)
        for annot in [k.lower() for k in Phd.CKEYWORDS]:
            value = None
            if annot == "trim":
                if record.annotations.get("trim", None):
                    value = "%s %s %.4f" % record.annotations["trim"]
            elif annot == "trace_peak_area_ratio":
                if record.annotations.get("trace_peak_area_ratio", None):
                    value = "%.4f" % record.annotations["trace_peak_area_ratio"]
            else:
                value = record.annotations.get(annot, None)
            if value or value == 0:
                self.handle.write("%s: %s\n" % (annot.upper(), value))

        self.handle.write("\nEND_COMMENT\nBEGIN_DNA\n")
        for i, site in enumerate(record.seq):
            if peak_locations:
                self.handle.write("%s %i %i\n" % (
                        site,
                        phred_qualities[i],
                        peak_locations[i])
                        )
            else:
                self.handle.write("%s %i\n" % (
                        site,
                        phred_qualities[i])
                        )

        self.handle.write("END_DNA\nEND_SEQUENCE\n")

if __name__ == "__main__" :
    import os
    print "Quick self test"
    handle = open("../../Tests/Phd/phd1")
    records = list(PhdIterator(handle))
    handle.close()
    outfile = "../../Tests/Phd/test.phd"
    handle = open(outfile, "w")
    pw = PhdWriter(handle)
    for record in records:
        pw.write_record(record)
    handle.close()
    os.remove(outfile)
    print "Done"
