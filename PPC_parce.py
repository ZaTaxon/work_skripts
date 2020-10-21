#Считывает таблицу очищенную в R (однострочную) и делает сводную таблицу по заменам
#Пример строки: 1	C_nivalis_2737	YF

#Важное:
#1. Считает с 0, а не с 1 позиции
#2. Делает анализ по вычещенной фасте без гэпов, так что могут быть проблемы с соотношением позиций, если конкатенат
#3. Он не берет в анализ очень вариабельные сайты

#Формат аутпута details:
#P - произошла одна и та же замена
#D - в одной позиции независимые замены
#R1-R2 - по каждой ветке информация. Важное! Дается пара аминокислот в режиме Было-Стало (т.е. 2)

TreeSAAP_file = 'PPC_one_string_clear'
#очищенный файл с информацией о заменах

TreeSAAP_result = 'PPC_summary_table.txt'
#файл, куда записывается полученная таблица

#Словари для всех интересующих нас вариантов: подземные по отдельности и наземные в целом
storage = {'terrestrial':{}, 'E_lutescens_4905':{}, 'E_fuscocapillus_4907':{}, 'P_schaposchnikovi_5136':{},
           'E_talpinus_1924':{}, 'L_mandarinus_NC_025283':{}, 'H_fertilis_4909':{}, 'T__subterraneus_4875':{}}

positions = []
#список со всеми позициями, где найдены замены у любых видов

#Виды, кого мы считаем подземными и кого хотем посмотреть по отдельности
undergrounds = {'E_lutescens_4905':'', 'E_fuscocapillus_4907':'', 'P_schaposchnikovi_5136':'',
           'E_talpinus_1924':'', 'L_mandarinus_NC_025283':'', 'H_fertilis_4909':'', 'T__subterraneus_4875':''}


#для расчета плотности замен на каждой позиции

def table_filling():
    #заполняем таблицу со всеми заменами. Наземные идут после названий и заполняются отдельно
    order = 2
    for species in undergrounds.keys():
        #заполнение подземных видов, для каждого вида свой столбец
        substitututions[0][order] = species
        for key, value in storage[species].items():
            j = positions.index(int(key)) + 1
            substitututions[j][order] = value
        order += 1
    substitututions[0][1] = 'terrestrial'
    for key, value in storage['terrestrial'].items():
        #заполняем наземных, удаляя дубликаты аминокислот в позиции и делая разделителем "/"
        j = positions.index(int(key)) + 1
        value = list(set(value))
        value = ('/').join(value)
        substitututions[j][1] = value

def insert_terrestrial(line, species):
    #записываем аминокислоты наземных видов в лист, т.к. их много на 1 позицию
    pos = line[0]
    AA = line[2][:-1]
    #print(species, pos, AA)
    positions.append(int(pos))
    if pos in storage['terrestrial'].keys():
        temp = storage['terrestrial'][pos]
        temp.append(AA)
        storage['terrestrial'][pos] = temp
    else:
        temp = []
        temp.append(AA)
        storage['terrestrial'][pos] = temp

def insert_to_dict(line, species):
    #каждую замену для подземных сохраняем в соответствующий словарь
    pos = line[0]
    AA = line[2][:-1]
    positions.append(int(pos))
    #print(species, pos, AA)
    storage[species][pos] = AA

with open(TreeSAAP_file, 'r') as file:
    #считываем файл построчно
    for line in file:
        line = line.split('\t')
        species = line[1]
        if species in undergrounds.keys():
            #проверяем, является ли вид подземным
            #print(line)
            insert_to_dict(line, species)
        else:
            #это наземные виды, их мы записываем отдельно
            #print(line)
            insert_terrestrial(line, species)

positions = list(set(positions))
positions.sort()

substitututions = [['*' for _ in range(len(storage)+1)] for j in range(len(positions)+1)]
substitututions[0][0] = 'Positions'

for elem in positions:
    substitututions[positions.index(elem) + 1][0] = str(elem + 1)
    #print(elem, positions.index(elem))

table_filling()

# for elem in substitututions:
#     print(elem)
#print(positions)

for elem in substitututions:
    print(elem)
    elem = ('\t').join(elem)
    with open(TreeSAAP_result, 'a') as answer:
        answer.write(elem)
        answer.write('\n')


# print(storage['E_lutescens_4905'])
# print(positions)
# print(storage['terrestrial'])