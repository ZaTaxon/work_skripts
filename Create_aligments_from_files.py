'Берет на вход номера ортогрупп, сливает воедино и добавляет номер группы в название сиквенса'

import os
import Bio.SeqIO as SeqIO
import re

directory = '/home/olga/Documents/Rodents/Under_transcropts_story/all scpecies/Lasyo/aligned/'
#Ссылка на папку, в которой лежат фасты (нуклеотидные) всех изучаемых видов

file_distribution = '/home/olga/Documents/Rodents/Under_transcropts_story/all scpecies/Lasyo/GO files/Distribution_groups'
#Файл типа: 262	First_common, чтобы было понятно, какой файл к какой группе относится

output_dir = '/home/olga/Documents/Rodents/Under_transcropts_story/all scpecies/Lasyo/GO files/'
#Папка, куда записать все созданные файлы

ortogroups = {}
#Словарь, где будут храниться записи о группах из файла file_distribution

def write_fasta(name, sequence, org_name, file_name):
    #Записываем сиквенс в файл согласно file_distribution
    name = org_name + '_' + name
    #print(name)
    answer_file = output_dir + file_name
    with open(answer_file, 'a') as answer:
        answer.write('>' + name + '\n' + sequence + '\n')

def read_fasta(file, org_name):
    if org_name in ortogroups.keys():
        #print('In keys')
        file_name = ortogroups[org_name]
        full_way = directory + file
        fasta_file = SeqIO.parse(open(full_way), 'fasta')
        for fasta in fasta_file:
            #берем из фасты первый сиквенс, отправляем его на запись и прерываем цикл. Другие сиквенсы из файла не нужны
            #Заодно удаляем из него гэпы
            name, sequence = fasta.id, str(fasta.seq.ungap('-'))
            write_fasta(name, sequence, org_name, file_name)
            break

#делаем словарь. который содержит все записи из файла file_distribution
with open(file_distribution, 'r') as file:
    for line in file:
        line = line.split()
        #print(line)
        ortogroups[line[0]] = line[1]
#print(ortogroups)

files = os.listdir(directory)
print(files)

#по очереди считываем фасты из папки и отправляем на анализ
for file in files:
    org_name = re.split(r'\.', file)[0]
    #print(org_name)
    read_fasta(file, org_name)
