"Берет на вход белковую таблицу ортологов и создает таблицу с ID транскриптов и количеством копий для белка"
# Пример таблицы в итоге. Столбцы: Prot- количество rna - rna- transcript
#['XP_003497743.1', ['1'], ['XM_003497695.5'], ['XM_003497695.5']]

import re
import itertools

#proteins = {}
#transcripts = {}
#order = {}

def open_function(protein_file):
    #записываем все белки из списка в словарь. Поскольу в конце есть перенос строки, то его убираем
    with open(protein_file, 'r') as ortologues:
        count = 1
        for line in ortologues:
            proteins[line[:-1]] = []
            order[line[:-1]] = count
            count += 1
    #print(len(proteins))

#создает словарь соответствия rna - transcript ID и protein ID - Parent
def distributor(rna, ids,value):
    #print(rna, ids, value)
    if value == 'protein' and ids in proteins.keys():
        print('Find it into ortologues')
        proteins[ids] += [rna]
    if value == 'transcript' and rna not in transcripts.keys():
        print('Add it to transcripts')
        transcripts[rna] = ids


#из gff файла выцеляет соответствие  белковое ID - Parent RNA
def retriever(gff_file):
    with open(gff_file, 'r') as gff:
        # for line in itertools.islice(gff, 8, None):
        for line in gff:
            #print(line.split(';'))
            # print(line)
            match_rna = re.compile(r'Parent=rna-([A-Za-z\_0-9\.]+)')
            if match_rna.findall(line):
                print(match_rna.findall(line)[0])
                rna = match_rna.findall(line)[0]
                match_protein = re.compile(r'protein_id=([A-Za-z\_0-9\.]+)')
                # print(line)
                if match_protein.findall(line):
                    ids = match_protein.findall(line)[0]
                    #print(ids)
                    distributor(rna, ids, 'protein')
                match_transcript = re.compile(r'transcript_id=([A-Za-z\_0-9\.]+)')
                if match_transcript.findall(line):
                    ids = match_transcript.findall(line)[0]
                    distributor(rna, ids, 'transcript')
                    #print(ids)
        #print(len(transcripts))

#Создает таблицу Prot- количество rna - rna- transcript для каждого вида

def table_divide(protein,i,length, new_value, table):
    table[i][0] = protein
    table[i][1] = list(str(length))
    table[i][2] = new_value
    for value in new_value:
        if value in transcripts.keys():
            table[i][3].append(transcripts[value])
        else:
            table[i][3] = 'Non transcripts found'

def table_creator(table_name):
    table = [[[] for i in range(4)] for _ in range(len(order))]
    for key, value in proteins.items():
        i = order[key] -1
        protein = key
        new_value = list(set(value))
        length = len(new_value)
        table_divide(protein,i,length, new_value, table)
    with open(table_name, 'w') as answer:
        for elem in table:
            answer.write(str(elem))
            answer.write('\n')
            print(elem)

def common_function(gff_file, protein_file, table_name):
    global proteins, transcripts, order
    proteins = {}
    transcripts = {}
    order = {}
    open_function(protein_file)
    retriever(gff_file)
    table_creator(table_name)

#На вход дается gff, файл с белковыми ID, соответствие которым нужно установить и название таблички на выход

common_function('/C_griseus/GCF_000223135.1_CriGri_1.0_genomic.gff',
    'Cricetulus.ID.txt',
    'Cricetulus.NCBI_prot-nucl.txt')

# Пример файла Cricetulus ID
# XP_003500610.1
# XP_003497743.1
# XP_003503636.1
# XP_003497721.1
# XP_003506454.1
# XP_003506257.1
# XP_003495843.1
# XP_007650601.1
# XP_007639342.1

