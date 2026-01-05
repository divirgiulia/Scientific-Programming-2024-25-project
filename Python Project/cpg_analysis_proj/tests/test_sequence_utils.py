import unittest
from sequence_utils import read_sequence, is_valid_dna, global_cpg, sliding_window_cpg, detect_cpg_islands

class TestSequenceUtils(unittest.TestCase):
    """Unit tests for sequence_utils.py functions."""

    def test_read_sequence_raw(self):
        """Test reading a raw DNA sequence string."""
        seq = 'ATCGATCG'
        self.assertEqual(read_sequence(seq), "ATCGATCG")
    
    def test_read_sequence_fasta(self):
        """Test reading a sequence from FASTA formatted input."""
        fasta = ">header\nATCGATCG"
        self.assertEqual(read_sequence(fasta), 'ATCGATCG')

    def test_is_valid_dna(self):
        """Test DNA validation function."""
        self.assertTrue(is_valid_dna("ATCG"))
        self.assertFalse(is_valid_dna("ATXBG"))
    
    def test_global_cpg(self):
        """Test global CpG content calculation."""
        seq = 'CGCGATAT'
        result = global_cpg(seq)
        self.assertAlmostEqual(result, 28.57)
    
    def test_sliding_window_cpg(self):
        """Test sliding window CpG content calculation."""
        seq = 'CGCGATATCGCG'
        results = sliding_window_cpg(seq, window=4, step=2)
        self.assertIsInstance(results, list)
        self.assertTrue(all('position' in r and 'cpg_percent' in r for r in results))

    def test_detect_cpg_islands(self):
        """Test detection of CpG islands."""
        seq = "CGCGCGCGCGCGATATATATATATAT"
        islands = detect_cpg_islands(seq, window=6, step=2, threshold=50)
        self.assertIsInstance(islands, list)

if __name__ == '__main__':
    unittest.main()

