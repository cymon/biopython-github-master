# Copyright 2008 by Peter Cock.  All rights reserved.
#
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Bio.SeqIO support for the "phd" file format.

PHD files are output by PHRED and used by PHRAP and CONSED.

You are expected to use this module via the Bio.SeqIO functions.
See also the underlying Bio.Sequencing.Phd module."""

import time
from Bio.SeqRecord import SeqRecord
from Bio.Sequencing import Phd
from Bio.SeqIO.Interfaces import SequentialSequenceWriter
from Bio.SeqIO import QualityIO

# Defaults as BioPerl - no defaults for 'trim',
# or 'trace_peak_area_ratio'
# 'chromat_file' and 'time' are set if necessary in PhdWriter.write_record()
ANNOTATION_DEFAULTS = {
                "abi_thumbprint": 0,
                "phred_version": "0.980904.e",
                "call_method": "phred",
                "quality_levels": 99,
                "trace_array_min_index": 0,
                "trace_array_max_index": "unknown",
                "chem": "unknown",
                "dye": "unknown"
                }

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
        seq_record.letter_annotations["peak_location"] = \
                [int(site[2]) for site in phd_record.sites]
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
        if not record.letter_annotations.has_key("peak_location"):
            raise ValueError("No suitable 'peak_location' found in letter_annotations "
                             "of SeqRecord (id=%s)." % record.id)
        assert len(record.seq) == \
                len(record.letter_annotations["phred_quality"]), "Number of " + \
                "phd quality scores does not match length of sequence"
        assert len(record.seq) == \
                len(record.letter_annotations["peak_location"]), "Number " + \
                "of peak location scores does not match length of sequence"
        self.handle.write("BEGIN_SEQUENCE %s\n\nBEGIN_COMMENT\n\n" % record.name)
        for annot in [k.lower() for k in Phd.CKEYWORDS]:
            if annot == "trim":
                # No default
                if record.annotations.get("trim", None):
                    value = "%s %s %.4f" % record.annotations["trim"]
                else:
                    continue
            elif annot == "trace_peak_area_ratio":
                # No default
                if record.annotations.get("trace_peak_area_ratio", None):
                    value = "%.4f" % record.annotations["trace_peak_area_ratio"]
                else:
                    continue
            elif annot == "time":
                if record.annotations.get("time", None):
                    value = record.annotations["time"]
                else:
                    value = time.strftime("%a %b %d %H:%M:%S %Y", time.gmtime())
            elif annot == "chromat_file":
                if not record.annotations.get("chromat_file", None):
                    value = record.id
                else:
                    value = record.annotations["chromat_file"]
            else:
                if not record.annotations.get(annot, None):
                    value = ANNOTATION_DEFAULTS[annot]
                else:
                    value = record.annotations[annot]
            self.handle.write("%s: %s\n" % (annot.upper(), value))

        self.handle.write("\nEND_COMMENT\n\nBEGIN_DNA\n")
        for i, site in enumerate(record.seq):
            self.handle.write("%s %i %i\n" % (
                    site,
                    phred_qualities[i],
                    record.letter_annotations["peak_location"][i])
                    )
        self.handle.write("END_DNA\n\nEND_SEQUENCE\n")

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
