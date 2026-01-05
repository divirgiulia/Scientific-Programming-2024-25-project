from sequence_utils import read_sequence, global_cpg,  is_valid_dna
import unittest

class TestSequence(unittest.TestCase):
    def test_long_sequence(self):
        seq = read_sequence("test_data/long_seq.txt")
        result = global_cpg(seq)
        assert isinstance(result, float), "Result should be a number"
        assert 0 <= result <= 100, "CpG percentage should be between 0 and 100"

    def test_short_sequence(self):
        seq = read_sequence("test_data/short_seq.txt")
        result = global_cpg(seq)
        assert isinstance(result, float), "Result should be a number"
        assert 0 <= result <= 100, "CpG percentage should be between 0 and 100"

    def test_invalid_sequence(self):
        seq = read_sequence("test_data/invalid_seq.txt")
        assert not is_valid_dna(seq), "Invalid sequence should be rejected"

if __name__ == '__main__':
    unittest.main()
