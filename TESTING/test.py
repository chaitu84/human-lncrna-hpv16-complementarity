import unittest
from pathlib import Path

from config_code.complements.shuffled import shuffled

class ShuffleTestCase(unittest.TestCase):
    def test_shuffle(self):
        file_path = Path(r"C:\Users\chait\Downloads\BINF\BINF_2111\src\config_code\hpv_step2\targets\HPV16_E6starI_E7_spliced_region.rna.fasta")
        
        # Verify the file exists before running the test (optional, but good practice)
        self.assertTrue(file_path.exists(), "FASTA file path does not exist.")
        
        # Call your function (pass the path if `shuffled` expects it)
        result = shuffled()
        self.assertTrue(result)
        print(result)

if __name__ == "__main__":
    unittest.main()