import os
import unittest

from Bio.Sequencing import Phd
from Bio import SeqIO
from Bio.SeqIO.PhdIO import ANNOTATION_DEFAULTS

class PhdTestOne(unittest.TestCase):
    def setUp(self):
        self.handle = open("Phd/phd1")

    def tearDown(self):
        self.handle.close()

    def test_check_record_parser(self):
        """Test to check that record parser parses all records of a contig.
        """
        records = Phd.parse(self.handle)
        # Record 1
        record = records.next()
        self.assertEqual(record.file_name, "34_222_(80-A03-19).b.ab1")
        self.assertEqual(record.comments['abi_thumbprint'], 0)
        self.assertEqual(record.comments['call_method'], "phred")
        self.assertEqual(record.comments['chem'], "term")
        self.assertEqual(record.comments['chromat_file'], "34_222_(80-A03-19).b.ab1")
        self.assertEqual(record.comments['dye'], "big")
        self.assertEqual(record.comments['phred_version'], "0.020425.c")
        self.assertEqual(record.comments['quality_levels'], 99)
        self.assertEqual(record.comments['time'], "Fri Feb 13 09:16:11 2004")
        self.assertEqual(record.comments['trace_array_max_index'], 10867)
        self.assertEqual(record.comments['trace_array_min_index'], 0)
        self.assertAlmostEqual(record.comments['trace_peak_area_ratio'], 0.1467)
        self.assertEqual(record.comments['trim'][0], 3)
        self.assertEqual(record.comments['trim'][1], 391)
        self.assertAlmostEqual(record.comments['trim'][2], 0.05)
        center = len(record.sites)/2
        self.assertEqual(record.sites[0], ('c', '9', '6'))
        self.assertEqual(record.sites[1], ('t', '9', '18'))
        self.assertEqual(record.sites[2], ('c', '10', '26'))
        self.assertEqual(record.sites[3], ('c', '19', '38'))
        self.assertEqual(record.sites[4], ('g', '22', '49'))
        self.assertEqual(record.sites[5], ('t', '37', '65'))
        self.assertEqual(record.sites[6], ('c', '28', '76'))
        self.assertEqual(record.sites[7], ('g', '28', '87'))
        self.assertEqual(record.sites[8], ('g', '24', '100'))
        self.assertEqual(record.sites[9], ('a', '22', '108'))
        self.assertEqual(record.sites[center-5], ('c', '11', '5259'))
        self.assertEqual(record.sites[center-4], ('c', '11', '5273'))
        self.assertEqual(record.sites[center-3], ('t', '9', '5286'))
        self.assertEqual(record.sites[center-2], ('g', '10', '5300'))
        self.assertEqual(record.sites[center-1], ('a', '10', '5316'))
        self.assertEqual(record.sites[center], ('t', '8', '5323'))
        self.assertEqual(record.sites[center+1], ('c', '8', '5343'))
        self.assertEqual(record.sites[center+2], ('g', '8', '5352'))
        self.assertEqual(record.sites[center+3], ('c', '8', '5366'))
        self.assertEqual(record.sites[center+4], ('c', '8', '5378'))
        self.assertEqual(record.sites[-10], ('c', '8', '10756'))
        self.assertEqual(record.sites[-9], ('c', '8', '10764'))
        self.assertEqual(record.sites[-8], ('a', '8', '10769'))
        self.assertEqual(record.sites[-7], ('a', '8', '10788'))
        self.assertEqual(record.sites[-6], ('a', '8', '10803'))
        self.assertEqual(record.sites[-5], ('g', '10', '10816'))
        self.assertEqual(record.sites[-4], ('c', '11', '10826'))
        self.assertEqual(record.sites[-3], ('g', '11', '10840'))
        self.assertEqual(record.sites[-2], ('t', '11', '10855'))
        self.assertEqual(record.sites[-1], ('g', '11', '10864'))
        self.assertEqual(record.seq.tostring()[:10], 'ctccgtcgga')
        self.assertEqual(record.seq.tostring()[-10:], 'ccaaagcgtg')
        self.assertEqual(record.seq_trimmed.tostring()[:10], 'cgtcggaaca')
        self.assertEqual(record.seq_trimmed.tostring()[-10:], 'tatttcggag')
        # Record 2
        record = records.next()
        center = len(record.sites)/2
        self.assertEqual(record.file_name, "425_103_(81-A03-19).g.ab1")
        self.assertEqual(record.comments['abi_thumbprint'], 0)
        self.assertEqual(record.comments['call_method'], 'phred')
        self.assertEqual(record.comments['chem'], 'term')
        self.assertEqual(record.comments['chromat_file'], '425_103_(81-A03-19).g.ab1')
        self.assertEqual(record.comments['dye'], 'big')
        self.assertEqual(record.comments['phred_version'], '0.020425.c')
        self.assertEqual(record.comments['quality_levels'], 99)
        self.assertEqual(record.comments['time'], 'Tue Feb 17 10:31:15 2004')
        self.assertEqual(record.comments['trace_array_max_index'], 10606)
        self.assertEqual(record.comments['trace_array_min_index'], 0)
        self.assertAlmostEqual(record.comments['trace_peak_area_ratio'], 0.0226)
        self.assertEqual(record.comments['trim'][0], 10)
        self.assertEqual(record.comments['trim'][1], 432)
        self.assertAlmostEqual(record.comments['trim'][2], 0.05)
        self.assertEqual(record.sites[0], ('c', '14', '3'))
        self.assertEqual(record.sites[1], ('g', '17', '11'))
        self.assertEqual(record.sites[2], ('g', '22', '23'))
        self.assertEqual(record.sites[3], ('g', '10', '35'))
        self.assertEqual(record.sites[4], ('a', '10', '53'))
        self.assertEqual(record.sites[5], ('t', '10', '68'))
        self.assertEqual(record.sites[6], ('c', '15', '75'))
        self.assertEqual(record.sites[7], ('c', '8', '85'))
        self.assertEqual(record.sites[8], ('c', '8', '94'))
        self.assertEqual(record.sites[9], ('a', '9', '115'))
        self.assertEqual(record.sites[center-5], ('c', '33', '5140'))
        self.assertEqual(record.sites[center-4], ('c', '28', '5156'))
        self.assertEqual(record.sites[center-3], ('g', '25', '5167'))
        self.assertEqual(record.sites[center-2], ('c', '28', '5178'))
        self.assertEqual(record.sites[center-1], ('c', '18', '5193'))
        self.assertEqual(record.sites[center], ('a', '16', '5204'))
        self.assertEqual(record.sites[center+1], ('a', '15', '5213'))
        self.assertEqual(record.sites[center+2], ('a', '10', '5230'))
        self.assertEqual(record.sites[center+3], ('a', '10', '5242'))
        self.assertEqual(record.sites[center+4], ('t', '8', '5249'))
        self.assertEqual(record.sites[-10], ('c', '8', '10489'))
        self.assertEqual(record.sites[-9], ('c', '8', '10503'))
        self.assertEqual(record.sites[-8], ('c', '8', '10514'))
        self.assertEqual(record.sites[-7], ('a', '8', '10516'))
        self.assertEqual(record.sites[-6], ('g', '8', '10530'))
        self.assertEqual(record.sites[-5], ('c', '8', '10550'))
        self.assertEqual(record.sites[-4], ('c', '10', '10566'))
        self.assertEqual(record.sites[-3], ('a', '8', '10574'))
        self.assertEqual(record.sites[-2], ('a', '7', '10584'))
        self.assertEqual(record.sites[-1], ('g', '7', '10599'))
        self.assertEqual(record.seq.tostring()[:10], 'cgggatccca')
        self.assertEqual(record.seq.tostring()[-10:], 'cccagccaag')
        self.assertEqual(record.seq_trimmed.tostring()[:10], 'cctgatccga')
        self.assertEqual(record.seq_trimmed.tostring()[-10:], 'ggggccgcca')
        # Record 3
        record = records.next()
        center = len(record.sites)/2
        self.assertEqual(record.file_name, '425_7_(71-A03-19).b.ab1')
        self.assertEqual(record.comments['abi_thumbprint'], 0)
        self.assertEqual(record.comments['call_method'], 'phred')
        self.assertEqual(record.comments['chem'], 'term')
        self.assertEqual(record.comments['chromat_file'], '425_7_(71-A03-19).b.ab1')
        self.assertEqual(record.comments['dye'], 'big')
        self.assertEqual(record.comments['phred_version'], '0.020425.c')
        self.assertEqual(record.comments['quality_levels'], 99)
        self.assertEqual(record.comments['time'], 'Thu Jan 29 11:46:14 2004')
        self.assertEqual(record.comments['trace_array_max_index'], 9513)
        self.assertEqual(record.comments['trace_array_min_index'], 0)
        self.assertAlmostEqual(record.comments['trace_peak_area_ratio'], 100.0)
        self.assertEqual(record.comments['trim'][0], -1)
        self.assertEqual(record.comments['trim'][1], -1)
        self.assertEqual(record.comments['trim'][2], 0.05)
        self.assertEqual(record.sites[0], ('a', '10', '7'))
        self.assertEqual(record.sites[1], ('c', '10', '13'))
        self.assertEqual(record.sites[2], ('a', '10', '21'))
        self.assertEqual(record.sites[3], ('t', '10', '28'))
        self.assertEqual(record.sites[4], ('a', '8', '33'))
        self.assertEqual(record.sites[5], ('a', '8', '40'))
        self.assertEqual(record.sites[6], ('a', '6', '50'))
        self.assertEqual(record.sites[7], ('t', '6', '53'))
        self.assertEqual(record.sites[8], ('c', '6', '66'))
        self.assertEqual(record.sites[9], ('a', '6', '68'))
        self.assertEqual(record.sites[center-5], ('a', '6', '4728'))
        self.assertEqual(record.sites[center-4], ('t', '10', '4737'))
        self.assertEqual(record.sites[center-3], ('a', '10', '4746'))
        self.assertEqual(record.sites[center-2], ('a', '8', '4756'))
        self.assertEqual(record.sites[center-1], ('t', '8', '4759'))
        self.assertEqual(record.sites[center], ('t', '8', '4768'))
        self.assertEqual(record.sites[center+1], ('a', '8', '4775'))
        self.assertEqual(record.sites[center+2], ('g', '10', '4783'))
        self.assertEqual(record.sites[center+3], ('t', '8', '4788'))
        self.assertEqual(record.sites[center+4], ('g', '8', '4794'))
        self.assertEqual(record.sites[-10], ('a', '8', '9445'))
        self.assertEqual(record.sites[-9], ('t', '6', '9453'))
        self.assertEqual(record.sites[-8], ('c', '6', '9462'))
        self.assertEqual(record.sites[-7], ('t', '6', '9465'))
        self.assertEqual(record.sites[-6], ('g', '6', '9478'))
        self.assertEqual(record.sites[-5], ('c', '6', '9483'))
        self.assertEqual(record.sites[-4], ('t', '6', '9485'))
        self.assertEqual(record.sites[-3], ('t', '8', '9495'))
        self.assertEqual(record.sites[-2], ('t', '3', '9504'))
        self.assertEqual(record.sites[-1], ('n', '0', '9511'))
        self.assertEqual(record.seq.tostring()[:10], 'acataaatca')
        self.assertEqual(record.seq.tostring()[-10:], 'atctgctttn')
        # Make sure that no further records are found
        self.assertRaises(StopIteration, records.next)

class PhdTestTwo(unittest.TestCase):
    def setUp(self):
        handle1 = open("Phd/phd1", "r")
        self.records = list(SeqIO.parse(handle1, "phd"))
        handle1.close()
        self.temp_filename = "Phd/phd_temp"
        handle2 = open(self.temp_filename, "w")
        SeqIO.write(self.records, handle2, "phd")
        handle2.close()
        handle3 = open(self.temp_filename, "r")
        self.new_records = list(SeqIO.parse(handle3, "phd"))
        handle3.close()

    def tearDown(self):
        os.remove(self.temp_filename)

    def test_check_record_writer(self):
        """Check that written records when parsed are sane
        """
        for i, record in enumerate(self.records):
            self.assertEqual(record.name, self.new_records[i].name)
            for annot in record.annotations:
                self.assertEqual(record.annotations[annot], 
                            self.new_records[i].annotations[annot])
            for j, site in enumerate(record.seq):
                self.assertEqual(record.seq[j], self.new_records[i].seq[j])
                self.assertEqual(record.letter_annotations['phred_quality'][j],
                        self.new_records[i].letter_annotations['phred_quality'][j])
                self.assertEqual(record.letter_annotations['peak_location'][j],
                        self.new_records[i].letter_annotations['peak_location'][j])

class PhdTestThree(unittest.TestCase):
    """The phd file Phd/Phd1_no_comments has no annotations in the COMMENTS
    block. By default, some of these should be written with default values.
    """

    def setUp(self):
        handle1 = open("Phd/phd1_no_comments", "rU")
        self.records = list(SeqIO.parse(handle1, "phd"))
        handle1.close()
        self.temp_filename = "Phd/phd_temp_no_comments"
        handle2 = open(self.temp_filename, "w")
        SeqIO.write(self.records, handle2, "phd")
        handle2.close()
        handle3 = open(self.temp_filename, "rU")
        self.new_records = list(SeqIO.parse(handle3, "phd"))
        handle3.close()

    def tearDown(self):
        os.remove(self.temp_filename)

    def test_check_annotation_default_values(self):
        """Check that annotation defaults are written"""
        for i, record in enumerate(self.records):
            for annot in record.annotations:
                if annot in ["trim", "trace_peak_area_ratio", "chromat_file",
                             "time"]:
                    # No default values
                    continue
                if ANNOTATION_DEFAULTS[annot] or ANNOTATION_DEFAULTS[annot] == 0:
                    self.assertEqual(ANNOTATION_DEFAULTS[annot], 
                        self.new_records[i].annotations[annot])
                else:
                    self.fail("Unknown annotation: %s" % annot)

if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity = 2)
    unittest.main(testRunner=runner)

