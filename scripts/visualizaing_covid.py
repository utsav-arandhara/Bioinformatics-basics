"""Module to analyze COVID-19 strains"""
 # pylint: disable=line-too-long

# Importing Biopythop
from Bio import SeqIO
import matplotlib.pyplot as plt

# Reading in COVID-19 strain sequences
fasta = list(SeqIO.parse('D:/Private/Dihang/Bioinformatics/Bioinformatics-basics/Bioinformatics-basics/datasets/seqs_aligned.fasta', format = "fasta"))

# Creating a disctionary for easier parsal
strain_seq = {}
for entry in fasta:
    strain_seq[entry.id] = entry

# All the strains of COVID-19 in the aligned fasta file
# for entry in fasta:
#     print(entry.id)

# Checking for the length of all the sequences in each strain
# for entry in fasta:
#     print(len(entry.seq))

# seq: aligned sequence; pos: position of spiked protein in un-aligned version.
def gapped_pos(seq, pos):
    """Function to find the spike protein site in the aligned sequences."""
    non_gaps = 0
    gaps = 0
    for nt in seq:
        if nt != '-':
            non_gaps += 1
        else:
            gaps += 1
        if non_gaps == pos:
            return pos + gaps

# Using one strain to find the spiked protein position
# print(gapped_pos(strain_seq['Wuhan_strain'].seq, 21563))
# print(gapped_pos(strain_seq['Wuhan_strain'].seq, 25384))

spikes = {}
for seq in fasta:
    spikes[seq.id] = seq.seq[21563-1:25393]

# for spike in spikes:
#     print(spikes[spike][0:10])

def get_mutations(initial, variant):
    """Module to get all the mutations in the strains"""
    zip_seqs = list(zip(initial, variant))
    for pos, nt in enumerate(zip_seqs):
        if nt[0] != nt[1]:
            print(nt[0].upper() + str(pos) + nt[1].upper())

# get_mutations(spikes['Wuhan_strain'], spikes['B.1.1.529|Omicron'])
# print(spikes)
# Next step is to work with the proteins directly, but we need to get rid of all the gaps in the aligned file before we translate

with open('datasets/spikes2.fasta', 'w', encoding="utf-8") as f:
    for strain, spike in spikes.items():
        out = spike.replace('-','').translate()
        f.write('>' + strain + '\n')
        f.write(str(out).upper() + '\n')

# Aligning done through EMBL-EBI MAFFT tool.        
spikes_align = list(SeqIO.parse('D:/Private/Dihang/Bioinformatics/Bioinformatics-basics/Bioinformatics-basics/datasets/spikes_aligned.fasta', format = "fasta"))
# print(spikes_align)

sp_align = {}
for entry in spikes_align:
    sp_align[entry.id] = entry.seq

print(sp_align['Wuhan_strain'])

# Re-using the mutations functions to find AA changes in different strains
def get_spikepro_mutations(initial, variant):
    """This function is to find the mutations in the protein changes for two strains"""
    out_list = []
    zip_seqs = list(zip(initial, variant))
    for pos, aa in enumerate(zip_seqs):
        if aa[0] != aa[1]:
            out_list.append(aa[0].upper() + str(pos) + aa[1].upper())
    return out_list

print(len(get_spikepro_mutations(sp_align['Wuhan_strain'], sp_align['B.1.1.529|Omicron'])))
# print(len(get_spikepro_mutations(sp_align['Wuhan_strain'], sp_align['B.1.1.7|Alpha'])))

# Loop to find the differences in spike proteins in each strain.
print("Comparisons of AA changes wrt original Wuhan Strain")
for item, seqs in sp_align.items():
    print(item + ':' + str(len(get_spikepro_mutations(sp_align['Wuhan_strain'], seqs))))

# Plotting the changes

for y, item in enumerate(sp_align):
    plt.plot((0, len(sp_align['Wuhan_strain'])), (y,y), color = 'lightgray')
plt.xlim(-150, len(sp_align['Wuhan_strain']) + 100)
plt.ylim(-.75, 5.75)
plt.show()
