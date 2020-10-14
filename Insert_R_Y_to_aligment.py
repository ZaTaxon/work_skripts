'Скрипт заменяет в выравнивании 3 позиции кодона на R или Y в зависимости от нуклеотида'

import Bio.SeqIO as SeqIO

aligment = '/home/olga/Documents/Rodents/mt_genome_story/full_MT/underground_mitogenomes/codon usage test/62_part1.fasta'
#Выравнивание, в котором нужно заменить позиции. Важно, чтобы выравнивание не сбивало рамку считывания и начиналось с 1 позиции

edit_aligment = '/home/olga/Documents/Rodents/mt_genome_story/full_MT/underground_mitogenomes/codon usage test/62_edit.fa'
#Выравнивание с исправленными 3 позициями кодона

nucleotides = { 'A':'R', 'G':'R', 'C':'Y', 'T':'Y'}
#словарь с обозначением, какой нуклеотид заменять какой буквой

#ШАГ ВСПОМОГАТЕЛЬНЫЙ
def write_in_file(text):
    #функция, которая все записывает в файл
    with open(edit_aligment, 'a') as answer:
        answer.write(text)

#ШАГ 2
def replace_nucleotide(nucl):
    if nucl in nucleotides.keys():
    #проверяем, есть ли эти нуклеотиды в словаре и нужно ли их заменять. Если нужно - отправляем на запись замененный
        write_in_file(nucleotides[nucl])
    else:
        #если это не нуклеотиды, а N или -, то отправляем на запись без изменений
        write_in_file(nucl)

#ШАГ 1
def count_nucleotides(sequence):
    #делаем счетчик по нуклеотидам. Проверяем позиции на то, что они 3.
    for i in range(0, len(sequence)):
        if i in range(2,len(sequence),3):
            #если позиция оказывается в списке третьих позиций, отправляем на замену нуклеотида
            nucl = sequence[i]
            replace_nucleotide(nucl)
        else:
            #если это 1 или 2 позиция, то просто пишем обратно в файл без изменений
            write_in_file(sequence[i])
    write_in_file('\n')

#ШАГ 0
fasta_file = SeqIO.parse(aligment, 'fasta')
for fasta in fasta_file:
    name, sequence = fasta.id, str(fasta.seq)
    # считываем по отдельности фасты для каждого вида из выравнивания
    # print(name)
    write_in_file('>')
    write_in_file(name)
    write_in_file('\n')
    #сразу записываем название сиквенса в файл, чтобы не забылось нигде
    count_nucleotides(sequence)
    #саму последовательность отправляем дальше