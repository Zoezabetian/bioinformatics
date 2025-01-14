## DNA Data Encoding and Decoding System

### Overview
This fictional project explores using DNA as long-term data storage. The program provides a system to encode binary data into DNA sequences with embedded tags for segmentation and error detection. It also decodes the tagged DNA sequence back into the original binary file. The implementation is designed for robust data integrity and efficient sequence parsing. This project was made with three other team members. 

---

### Features
1. **Encoding**:
   - Converts binary data into a DNA sequence using a nucleotide substitution cipher.
   - Inserts location-specific tags every 256 nucleotides for segmentation and easier error tracking.
   - Adds non-coding filler nucleotides to maintain consistent segment lengths.

2. **Decoding**:
   - Removes location-specific tags and filler nucleotides from DNA sequences.
   - Accurately reconstructs the original binary file using the substitution cipher.

3. **Error Handling**:
   - Validates alignment of location tags.
   - Ensures sequence integrity by removing extraneous nucleotides.

---

### Usage
1. **Encoding**:
   - Initialize the `Encoder` class with the binary file to be encoded.
   - Call the `addTags` method to encode the file into a tagged DNA sequence.

2. **Decoding**:
   - Initialize the `Decoder` class with the tagged DNA sequence.
   - Use the `convertToBin` method to decode the sequence back into the original binary format.
