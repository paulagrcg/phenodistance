from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML

# Define the RNA sequences to search
sequences = {
    "FR456139": "GCAGAGCTGCAGCGATGCGGA",
    "FR422481": "AGTCAAACGGCAAAACGG",
    "FR132646|DQ703740": "CTATCCAACCCAGATGGG",
    "FR412251": "AAGGAGCATCGTCCACCGGTG",
    "FR446582": "CTCTTTAGCTCAGTGGAGAGCACT",
}

# Perform BLAST search for each sequence
for identifier, sequence in sequences.items():
    #print(f"Running BLAST for {identifier}...")
    result_handle = NCBIWWW.qblast("blastn", "nt", sequence)
    
    # Save the BLAST result as an XML file
    file_name = f"{identifier}_blast.xml"
    with open(file_name, "w") as save_file:
        save_file.write(result_handle.read())
    print(f"Results saved to {file_name}")
    with open(file_name) as result_handle:
        blast_record = NCBIXML.read(result_handle)

    #with open(file_name) as result_handle:
    #    blast_record = NCBIXML.read(result_handle)

    with open(file_name+'_analysis.txt', 'w') as file:
        for alignment in blast_record.alignments:
            for hsp in alignment.hsps:
                file.write(f"Sequence: {alignment.title}\n")
                file.write(f"Length: {alignment.length}\n")
                file.write(f"Score: {hsp.score}\n")
                file.write(f"E-value: {hsp.expect}\n")
                file.write(f"Alignment:\n{hsp.query}\n{hsp.match}\n{hsp.sbjct}\n\n")
