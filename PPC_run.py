#https://pypi.org/project/ProtParCon/

from ProtParCon import imc

sequence = '/under_protein_coding_translation.fasta'
#Привести к единообразным названиям все в дереве и белковой фасте
#Лишнее в скобках удаляем выраженем \(.*\) в выравнивании

tree = 'under_aligment_consensus_tree.newick'
#Все узловые значения (типа 100.0 и другие) заменить на NODE с цифрой
#Например (E_luteus_1816:0.06526060298578323,L_lagurus_5135:0.06278815662600405)NODE1:0.0251830230723456

muscle = 'muscle'
codeml = 'codeml'
evolver = '/paml4.9g/bin/evolver'

imc(sequence, tree, aligner=muscle, ancestor=codeml, simulator=evolver)

