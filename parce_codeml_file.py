'Считывает файл codeml и достает оттуда информацию о велечине dN/dS на филогенетическом дереве'
#таблица создается формата ген - организм - значение

import os
import re

directory = 'results codeml/'
#Папка, из котрой будут считываться результаты

table_name = 'Codeml_free_model_results.txt'
#Таблица, в которую будут записываться все результаты

label = False

def parce_line(line, ortho_group):
    line = line.split(',')
    print(line)
    for elem in line:
        elem = elem.split(' ')
        # print(elem)
        # print(elem[1])
        # print(elem[2][1:])
    # dn = re.findall(r'\d+.\d+', line)
    # print(dn)
    #нахдим все цифры в строке. Это и есть значения
        with open(table_name, 'a') as answer:
            #записываем их в таблицу
            answer.write(ortho_group + '\t')
            answer.write(elem[1] + '\t' + elem[2][1:])
            answer.write('\n')


def find_line(file, ortho_group):
    #Функция читает файл до тех пор, пока не найдет нужную строку с значениями. При этом меняется label
    global label
    way = directory + file
    with open(way, 'r') as f:
        for line in f:
            if label:
                print(line)
                line = ' ' + line
                #добавляем пробел в начале строки, чтобы было удобно дальше парсить
                label = False
                parce_line(line, ortho_group)
                #меняем значение label, а саму строку передаем на парсинг
            if line == 'w ratios as labels for TreeView:\n':
                #мы нашли нужную строку, следующая будет со значениями dN/dS
                label = True

files = os.listdir(directory)
#создаем список файлов, которые будем открывать

for file in files:
    ortho_group = re.split(r'_', file)[2]
    find_line(file, ortho_group)
    #циклично считываем каждый файл, храним название гена для таблицы и передаем дальше