"""Считывает таблицу с ортологичными группами и создает отдельные файлы для каждой из них (по каждому ряду создается
файл с порядковым номером этого ряда, в котором соответствующие сиквенсы всех видов
Важно! Первое слово в названии фасты до знака _ должно совпадать с тем, что есть в списке видов!!"""

import os
import Bio.SeqIO as SeqIO
import re

ortho_file = 'Arvicolinae_RNA.proteinortho.clean_with_NCBI_nucl.txt'
#Сюда нужно вставить ссылку на таблицу с очищенными ортологическими рядами. Сейчас разделитель запятая. можно менять

species_list = ['A.fortis',	'A.lemminus', 'A.amphibius',	'C.griseus',	'C.nivalis',
'C.glareolus',	'D.torquatus',	'E.lutescens',	'L.brandtii',	'L.gregalis',
'L.mandarinus',	'L.raddei',	'L.sibiricus',	'M.ochrogaster', 'M.pennsylvanicus',	'M.schisticolor',
'O.zibethicus',	'P.schaposchnikovi', 'T.subterraneus']
#Здесь нужно указать порядок видов, как он идет в таблице с ортологическими рядами

species = [{}, {}, {}, {}, {}, {},{}, {}, {}, {}, {}, {},{}, {}, {}, {}, {}, {}, {}]
#В листе сделать нужное количество словарей, соответствующее количеству изучаемых видов

directory = 'all scpecies/fasta/'
#Ссылка на папку, в которой лежат cds фасты (нуклеотидные) всех изучаемых видов

output_dir = 'all scpecies/clusters'
#Ссылка на папку, куда будут сохраняться результаты

def write_dictionaries(line, count):
    'Заполняет файл со словарями по порядку, присваивая каждому гену номер ортологической группы'
    org = 0
    for genes in line:
        if org != len(line)-1:
            name = species[org]
            name[genes] = count
        else:
            name = species[org]
            genes = genes[:-1]
            name[genes] = count
        org += 1

def create_dictionaries():
    "Читает таблицу с ортологическими рядами и по очереди отдает их функции write_dictionaries"
    with open(ortho_file, 'r') as ortologues:
        count = 1
        for line in ortologues:
            line = line.split(',')
            write_dictionaries(line, count)
            count += 1

def write_clusters(name, sequence, ortho_group, org_name):
    #записывает все в файлы. Номер файла соответствует номеру ортогруппы
    os.chdir(output_dir)
    cluster_name = str(ortho_group) + '.fa'
    with open(cluster_name, 'a') as answer:
        answer.write('>' + org_name)
        answer.write('\n')
        answer.write(sequence)
        answer.write('\n')

def read_fasta(file, org_number, org_name):
    "Читает фасты изучаемых видов и ищет нужные транскрипты в созданных словарях"
    full_way = directory + file
    fasta_file = SeqIO.parse(open(full_way), 'fasta')
    for fasta in fasta_file:
        name, sequence = fasta.id, str(fasta.seq)
        #При работе только с protheinortho 2 последних знака можно будет не убирать
        my_d = species[org_number]
        if name in my_d.keys():
            ortho_group = my_d[name]
            write_clusters(name, sequence, ortho_group, org_name)

create_dictionaries()
# for elem in species:
#     print(elem)

files = os.listdir(directory)
print(files)

for file in files:
    org_name = re.split(r'_', file)[0]
    print(org_name)
    org_number = species_list.index(org_name)
    read_fasta(file, org_number,org_name)



