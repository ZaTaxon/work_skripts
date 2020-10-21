#https://pypi.org/project/ProtParCon/

from ProtParCon import imc

sequence = '/under_protein_coding_translation.fasta'
tree = 'under_aligment_consensus_tree.newick'
muscle = 'muscle'
codeml = 'codeml'
evolver = '/paml4.9g/bin/evolver'

imc(sequence, tree, aligner=muscle, ancestor=codeml, simulator=evolver)

