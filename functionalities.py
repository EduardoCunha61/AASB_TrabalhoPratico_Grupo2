from classes import Seq
import re
from os import path

# Funcao que guarda novas propriedades para uma dada proteina
def edit_protein(protein):
    aux = []
    keys = []

    for prop in protein.propriedades:
        if(prop!='Size'):
            keys.append(prop)
            print('\t' + prop + ': ',end='')
            x = input()
            aux.append(x)
    
    updated_props = create_properties(protein.sequencia,keys,aux)
    res = add_protein(protein.sequencia,updated_props)

    return res

# Cria o dicionario com as propriedades inseridas
def create_properties(seq,keys,props):
    res = {}

    if(len(keys)==len(props)):
        for index in range(len(keys)):
            res[keys[index]] = props[index]
    
    res['Size'] = len(seq)

    return res

# Cria um objeto novo e valida a sequencia deste
def add_protein(seq,properties):
    obj = Seq(seq,properties)
    if(Seq.valida(obj)):
        return obj
    else:
        return None

# Le o ficheiro filename. Path.exists verifica se este ja existe, caso seja verdade retorna a informacao la dentro
def read_file(filename):
    if path.exists(filename):
        f = open(filename, "r")
        current = f.readlines()
        f.close()
    else:
        return 0

    return current

# Importa o ficheiro com as sequencias que foram inseridas na aplicacao
def import_db(filename):
    info = read_file(filename)
    if info:
        res = {}
        for line in info:
            obj = {}
            key = line.split(',')[0]
            seq = line.split(',')[1]
            props_str = re.search("{.+}",line).group()
            props = eval(props_str)
            obj = Seq(str(seq),props)
            if(Seq.valida(obj)):
                res[key] = obj
            else:
                return None

        return res

    else:
        print("File not found!")
    
# Exporta todas as sequencias inseridas na aplicacao para o ficheiro com nome output
def export(database,output):
    if read_file(output)!=0:
        c = read_file(output)
    else:
        c = ""

    f = open(output,"w")
    for index,protein in database.items():
        f.write(str(index) + ',' + str(protein.sequencia) + ',' + str(protein.propriedades) + '\n')
    f.close()

    l = read_file(output)

    if c==l: 
        return 0
    else:
        return 1

def read_fasta(filename):
    from re import sub
    fh = open(filename)
    lines = fh.readlines()
    sequence = "" 
    for line in lines[1:]:
        sequence += sub('\s','',line)
    fh.close()
    return sequence

def read_txt(filename):
    from re import sub
    fh = open(filename)
    lines = fh.readlines()
    sequence = "" 
    for line in lines:
        sequence += sub('\s','',line)
    fh.close()
    return sequence












    
    


