import unittest
from analysis import run_analysis

class TestAnalysis(unittest.TestCase):
    """Unit tests for analysis.py's run_analysis function."""

    def setUp(self):
        self.seq = "CGCGATATCGCG"
        
    def test_global_task(self):
        """Test global CpG content analysis."""
        result = run_analysis('global', self.seq)
        self.assertIn('cpg_percentage', result)

    def test_sliding_task(self):
        """Test sliding window analysis."""
        result = run_analysis('sliding', self.seq, window=4, step=2, threshold=60)
        self.assertIn('sliding_results', result)
        self.assertEqual(result['window'], 4)
        self.assertEqual(result['step'], 2)

    def test_island_cpg(self):
        """Test CpG island detection."""
        result = run_analysis('island', self.seq, window=4, step=2, threshold=60)
        self.assertIn('cpg_islands', result)
        self.assertEqual(result['threshold'], 60)

    def test_missing_params(self):
        """Test that missing parameters raise ValueError."""
        with self.assertRaises(ValueError):
            run_analysis('sliding', self.seq)
        with self.assertRaises(ValueError):
            run_analysis('island', self.seq, window=4, step=2)

    def test_unknown_task(self):
        """Test that an unknown task raises ValueError."""
        with self.assertRaises(ValueError):
            run_analysis('unknown', self.seq)

    if __name__ == '__main__':
        unittest.main()


