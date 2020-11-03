from ete3 import EvolTree
import re
import os

tree_file = '/home/olga/Documents/Rodents/Under_transcropts_story/all scpecies/trees/NJ_transcript.nwk'
#Ссылка на дерево в формате newik

directory = '/home/olga/Documents/Rodents/Under_transcropts_story/all scpecies/aligned/'
#Ссылка на папку, в которой лежат выровненные нуклеотидные фасты

#align_file = '/home/olga/Documents/Rodents/CytB+COI/cytB_62_spes/easyCodeML/Lasy/test_ali'
#если хочется обсчитать отдельное выравнивание, можно вставить путь до него сюда

temp = '/home/olga/Documents/Rodents/Under_transcropts_story/all scpecies/temp'
#временный файл для хранения результатов

final_table = '/home/olga/Documents/Rodents/Under_transcropts_story/all scpecies/free_branch_table'
#итоговая таблица с результатами. Формат записи результата: ген - организм - омега

problem_files = '/home/olga/Documents/Rodents/Under_transcropts_story/all scpecies/Lasyo/free_branch_problems'
#тут будет список файлов, с которыми возникли ошибки при подсчете отбора (из-за некратности 3)

def write_in_table(gene_name):
    with open(temp, 'r') as temp_lines:
        for line in temp_lines:
            line = line.split()
            # print(line[0])
            # print(line)
            if len(line) > 5:
                if line[10] == 'ROOT' or line[10] == 'EDGE':
                    continue
                else:
                    with open(final_table, 'a') as answer:
                        answer.write(gene_name + '\t')
                        answer.write(line[10] + '\t' + line[4] + '\n')
                    print(line)
                    # print(line[10], line[4])

def count_omega(align_file, gene_name):
    print(gene_name)
    tree = EvolTree(tree_file)
    tree.link_to_alignment(align_file)
    #
    # #free branch ratio count
    tree.run_model('fb')
    fb_results = tree.get_evol_model('fb')
    print(fb_results)
    with open(temp, 'w') as temp_file:
        temp_file.write(str(fb_results))
    write_in_table(gene_name)


files = os.listdir(directory)
print(files)

for file in files:
    gene_name = re.split(r'\.', file)[0]
    file_aligment = directory + file
    # count_omega(file_aligment, gene_name)
    try:
        count_omega(file_aligment, gene_name)
    except IndexError:
        print('Have problem with',file)
        with open(problem_files, 'a') as prob_file:
            prob_file.write(file + '\n')
        continue