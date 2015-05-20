
import pybedtools

class TestSpliceAnnotator(object):

    def test_init(self, miso_ids_splice_type):
        from rnaseek import miso
        miso_ids, splice_type = miso_ids_splice_type
        genome_fasta = pybedtools.example_filename('test.fa')
        sa = miso.SpliceAnnotator(miso_ids, splice_type,
                                  'test', genome_fasta)