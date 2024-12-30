from Bio.Blast import NCBIXML

# Parse a single BLAST XML result
file_name = "FR446582_blast.xml"
with open(file_name) as result_handle:
    blast_record = NCBIXML.read(result_handle)

# Print top alignments
print(f"Results for {file_name}:")
for alignment in blast_record.alignments:
    for hsp in alignment.hsps:
        print(f"Sequence: {alignment.title}")
        print(f"Length: {alignment.length}")
        print(f"Score: {hsp.score}")
        print(f"E-value: {hsp.expect}")
        print(f"Alignment:\n{hsp.query}\n{hsp.match}\n{hsp.sbjct}\n")