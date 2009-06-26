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
        seq_record.letter_annotations["phd_qualities"] = \
                [site[1] for site in phd_record.sites]
        seq_record.letter_annotations["peak_locations"] = \
                [site[2] for site in phd_record.sites]
        yield seq_record 
    #All done

class PhdWriter(SequentialSequenceWriter):
    """Class to write Phd format files"""

    def __init__(self, handle):
        SequentialSequenceWriter.__init__(self, handle)
        # Defaults as BioPerl - no defaults for 'trim',
        # or 'trace_peak_area_ratio'
        # 'chromat_file' and 'time' are set if necessary in write_record()
        self.header_defaults = {
                    "abi_thumbprint": 0,
                    "phred_version": "0.980904.e",
                    "call_method": "phred",
                    "quality_levels": 99,
                    "trace_array_min_index": 0,
                    "trace_array_max_index": "unknown",
                    "chem": "unknown",
                    "dye": "unknown"
                    }

    def write_record(self, record):
        """Write a single Phd record to the file."""
        self.handle.write("BEGIN_SEQUENCE %s\n\nBEGIN_COMMENT\n\n" % record.name)
        for keyword in Phd.CKEYWORDS:
            if keyword == "TRIM":
                # No default
                if record.annotations.has_key('trim'):
                    value = "%s %s %.4f" % record.annotations['trim']
                else:
                    continue
            elif keyword == "TRACE_PEAK_AREA_RATIO":
                # No default
                if record.annotations.has_key("trace_peak_area_ratio"):
                    value = "%.4f" % record.annotations["trace_peak_area_ratio"]
                else:
                    continue
            elif keyword == "TIME":
                value = record.annotations.get("time", 
                            time.strftime("%a %b %d %H:%M:%S %Y", time.gmtime()))
            elif keyword == "CHROMAT_FILE":
                value = record.annotations.get("chromat_file", record.id)
            else:
                value = record.annotations.get(keyword.lower(),
                            self.header_defaults[keyword.lower()])
            self.handle.write("%s: %s\n" % (keyword, value))

        self.handle.write("\nEND_COMMENT\n\nBEGIN_DNA\n")
        assert record.seq, "No sequence present in SeqRecord"
        assert record.letter_annotations["phd_qualities"], "No Phd quality " + \
                "scores in the per letter annotations"
        assert record.letter_annotations["peak_locations"], "No peak " + \
                "location values in the per letter annotations"
        assert len(record.seq) == \
                len(record.letter_annotations["phd_qualities"]), "Number of " + \
                "phd quality scores does not match length of sequence"
        assert len(record.seq) == \
                len(record.letter_annotations["peak_locations"]), "Number " + \
                "of peak location scores does not match length of sequence"
        for i, site in enumerate(record.seq):
            self.handle.write("%s %s %s\n" % (
                    site,
                    record.letter_annotations["phd_qualities"][i],
                    record.letter_annotations["peak_locations"][i])
                    )
        self.handle.write("END_DNA\n\nEND_SEQUENCE\n")

if __name__ == "__main__" :
    import os
    print "Quick self test"
    handle = open("../../Tests/Phd/phd1")
    records = list(PhdIterator(handle))
    #for record in records:
    #    print record
    handle.close()
    outfile = "../../Tests/Phd/test.phd"
    handle = open(outfile, "w")
    pw = PhdWriter(handle)
    for record in records:
        pw.write_record(record)
    handle.close()
    os.remove(outfile)
    print "Done"
