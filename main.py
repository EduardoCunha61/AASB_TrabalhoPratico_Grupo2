import os
import re
import functionalities
from classes import Seq, MatrixNum, SubstMatrix, MyAlign, AlignSeq, MultipleAlign, BinaryTree, ClustHier, UPGMA, Filogenia, AlinhamentoMultiplo, SearchDomain, Blast, SearchPTM, SaveSequence

#Classe seq cria o objeto para cada proteina
#O script functionalities e onde esta o codigo todo importante
#Esta classe apenas tem I/O e insercao da bd

inp = "sequences.txt"       #Ficheiro que contem todas as sequencias inseridas
output = "sequences.txt"    #Ficheiro onde vai ser feito o export (neste momento estao iguais para ler e escrever para o mesmo sitio)

clear_shell = os.system('cls' if os.name == 'nt' else 'clear') #Comando para limpar a shell

def create_new_protein(bd,inserting_id):
    print("Sequence: ",end='')
    sequence = input()
    s = Seq(sequence)
    if s.valida():
        print("Properties:")
        list_of_props = ["Id NCBI","Id Uniprot","Function","Genes","Taxonomic lineage"]
        aux = []
        for elem in list_of_props:
            print('\t' + elem + ': ',end='')
            x = input()
            aux.append(x)
        
        props_dic = functionalities.create_properties(sequence,list_of_props,aux)
        obj = functionalities.add_protein(sequence,props_dic)
        flag = 1
        if obj:
            for key in bd.keys():
                if obj.sequencia == bd[key].sequencia:
                    print("Sequence already in database! \n")
                    flag = 0
                    break
            if flag:
                bd[inserting_id] = obj
                if(obj in bd.values()):
                    clear_shell
                    print("Sequence saved!\n\n")
    else:
        print("Invalid Sequence!\n")
    
def main():
    bd = {}
    if(functionalities.read_file(inp)):
        if(functionalities.import_db(inp)!=None):
            bd = functionalities.import_db(inp)
        else:
            print("There are invalid sequences in " + str(inp))
    else:
        print(str(inp) + " is empty!")

    done = 1

    while(done):
        id_bd=len(bd)
        
        clear_shell
        print("\t\t\tProtein sequence manager")
        print("1. Add manually sequence")
        print("2. Read and validate a sequence from a file")
        print("3. Consult all sequences")
        print("4. Consult sequence by id")
        print("5. Pattern search")
        print("6. Export all sequences to " + output)
        print("7. Subsequences of size k search")
        print("8. Search the most similar sequence in database")
        print("9. Filogenetic Tree from a list of sequences")
        print("10. Multiple Alignment and Phylogenetic Tree")
        print("11. Search possible domains for an input sequence")
        print("12. Import sequences from NCBI protein using NCBI id")
        print("13. Make BLAST against an external database")
        print("14. Multiple Alignment using external programs")
        print("15. Phylogenetic Tree using external programs")
        print("16. Detect post transcriptional modifications using external programs")
        print("0. Exit\n")
        print("Insert option: ", end='')
        x = input()

        if(x=='0'):
            x = functionalities.export(bd,output)
            if x:
                print("Information exported to " + str(output) + '!\n' )
            else:
                print("Nothing added to file!\n")

            done=0

        elif(x=='1'):
            create_new_protein(bd,id_bd)

        elif(x=='2'):
            file_type = input('Is your file a text file, or a FASTA file? [txt/FASTA] ')
            filename = input('File name: ')
            print()
            if (file_type=='txt'):
                seq = functionalities.read_txt(filename)
                sequence = Seq(seq)
                if sequence.valida():
                    print('Your input sequence is a valid protein sequence!')
                else:
                    print('Your input sequence is not a valid protein sequence!')
                    
            elif (file_type=='FASTA'):
                seq = functionalities.read_fasta(filename)
                sequence = Seq(seq)
                if sequence.valida():
                    print('Your input sequence is a valid protein sequence!')
                else:
                    print('Your input sequence is not a valid protein sequence!')
            
            else:
                print('Invalid option!')
            
            pass

        elif(x=='3'):
            number_of_prints = 0
            print("Saved sequences:")
            for key,obj in bd.items():
                print("\t" + str(key) + ': ' + str(obj.sequencia))
                for prop_key in obj.propriedades:
                    print("\t\t" + str(prop_key) + ': ' + str(obj.propriedades[prop_key]))
                
                number_of_prints += 1                                       #Estas linhas para garantir que imprime 4 de cada vez, para ser mais facil de ler
                if number_of_prints>=4 and (len(bd)-number_of_prints)!= 0:
                    input("\nPress Enter to continue... " + str((len(bd)-number_of_prints)) + " more to show!\n")
                    number_of_prints = 0
                        
        
        elif(x=='4'):
            protein_id = input("Insert protein id: ")
            if(protein_id in bd.keys()): 
                print("\t" + str(protein_id) + ': ' + str(bd[protein_id].sequencia))
                for prop_key in bd[protein_id].propriedades:
                    print("\t\t" + str(prop_key) + ': ' + str(bd[protein_id].propriedades[prop_key]))
                
                option = input("Edit protein [y/n]?")
                if(option=='y' or option=='Y'):
                    res = functionalities.edit_protein(bd[protein_id])
                    if(res):
                        bd[protein_id] = res
            else:
                print("Invalid Id")
        
        elif(x=='5'):
            re_pattern = input("Insert regex pattern: ")
            pattern = re.compile(re_pattern)

            ans = input("All db or one sequence only [db/id]? ") #Escola entre pesquisar padrao na base de dados toda ou apenas numa dada sequencia
            
            if(ans == 'db'):
                for key in bd.keys():
                    number_of_prints = 0
                    all_patterns = pattern.finditer(bd[key].sequencia)
                    if(all_patterns):
                        print("Protein id: " + key)
                        for pat in all_patterns:
                            print("\t" + str(pat.group()) + " found in between positions " + str(pat.span()))
                            number_of_prints += 1
                        print("\nFound pattern " + str(number_of_prints) + " times in protein " + str(key) + '!\n')

            elif(ans.isdigit()):
                if ans in bd.keys():
                    number_of_prints = 0                    
                    all_patterns = pattern.finditer(bd[ans].sequencia)
                    if(all_patterns):
                        print("Protein id: " + ans)
                        for pat in all_patterns:
                            print("\t" + str(pat.group()) + " found in between positions " + str(pat.span()))
                            number_of_prints += 1
                        print("\nFound pattern " + str(number_of_prints) + " times in protein " + str(ans) + '!\n')
                    else:
                        print("Pattern not found in protein " + str(ans) + '!')
                else:
                    print("Invalid id!")

        elif(x=='6'):
            x = functionalities.export(bd,output)
            if x:
                print("Information exported to " + str(output) + '!\n' )
            else:
                print("Nothing added to file!\n")
            

        elif(x=='7'):
            ans = input('All db or one sequence only [db/id]? ')
            k = int(input('Insert an integer, k: '))
            mapa = {}
            if (ans=='db'):
                for ID in bd:
                    seq = bd[ID]
                    for i in range(len(seq)-k+1):
                        subseq = seq[i:i+k]
                        if subseq in mapa:
                            mapa[subseq] +=1
                        else:
                            mapa[subseq] = 1
            elif (ans=='id'):
                ID = str(input('Insert an id: '))
                if ID in bd:
                    seq = bd[ID]
                    for i in range(len(seq)-k+1):
                        subseq = seq[i:i+k]
                        if subseq in mapa:
                            mapa[subseq] +=1
                        else:
                            mapa[subseq] = 1
                else:
                    print('Invalid id!')
            else:
                print('Invalid option!')
            
            for pattern in mapa.keys():
                print('Pattern %s occurs %d times. \n' %(pattern,mapa[pattern]))
            print()
            
        elif(x=='8'):
            option = input('Input a sequence or input a file name that contains the sequence [seq/file]? ')
            
            if (option=='seq'):
                s = input('Sequence: ')
                seq1 = s.upper()
                sm = input('Do you want to use BLOSUM62 as substitution matrix? [yes/no] ')
                gap = int(input('What is the value for penalty gap? '))
                if (sm=='yes'):
                    scores = []
                    id_list = []
                    smatrix = SubstMatrix()
                    s_matrix = smatrix.loadFromFile('blosum62.mat','\t')
                    for ID in bd:
                        seq2 = bd[ID]
                        alseq = AlignSeq(s_matrix, gap)
                        alseq.needlemanWunsch(seq1, seq2)
                        align = alseq.recoverAlignment()
                        score = 0
                        for i in range(len(align[0])):
                            score += int(alseq.scorePos(align[0][i], align[1][i]))
                        scores.append(score)
                        id_list.append(ID)
                        
                    max_score = max(scores)
                    index = scores.index(max_score)
                    id_seq = id_list[index]
                    sequence = bd[id_seq]
                    
                    print('The most similar sequence in database is %s.\n' %sequence)
                        
                elif (sm=='no'):
                    match = int(input('Score para match: '))
                    mismatch = int(input('Score para mismatch: '))
                    scores = []
                    id_list = []
                    smatrix = SubstMatrix()
                    s_matrix = smatrix.createFromMatchPars(match, mismatch)
                    for ID in bd:
                        seq2 = bd[ID]
                        alseq = AlignSeq(s_matrix, gap)
                        alseq.needlemanWunsch(seq1, seq2)
                        align = alseq.recoverAlignment()
                        score = 0
                        for i in range(len(align[0])):
                            score += int(alseq.scorePos(align[0][i], align[1][i]))
                        scores.append(score)
                        id_list.append(ID)
                        
                    max_score = max(scores)
                    index = scores.index(max_score)
                    id_seq = id_list[index]
                    sequence = bd[id_seq]
                        
                    print('The most similar sequence in database is %s.\n' %sequence)
                        
                else:
                    print('Invalid option!')
                    
            elif (option=='file'):
                file_type = input('Is it a FASTA file or a simple txt file only with the sequence? [FASTA/txt]')
                filename = input('File name: ')
                if (file_type=='FASTA'):
                    s = functionalities.read_fasta(filename)
                    seq1 = s.upper()
                    sm = input('Do you want to use BLOSUM62 as substitution matrix? [yes/no] ')
                    gap = int(input('What is the value for penalty gap? '))
                    if (sm=='yes'):
                        scores = []
                        id_list = []
                        smatrix = SubstMatrix()
                        s_matrix = smatrix.loadFromFile('blosum62.mat','\t')
                        for ID in bd:
                            seq2 = bd[ID]
                            alseq = AlignSeq(s_matrix, gap)
                            alseq.needlemanWunsch(seq1, seq2)
                            align = alseq.recoverAlignment()
                            score = 0
                            for i in range(len(align[0])):
                                score += int(alseq.scorePos(align[0][i], align[1][i]))
                            scores.append(score)
                            id_list.append(ID)
                            
                        max_score = max(scores)
                        index = scores.index(max_score)
                        id_seq = id_list[index]
                        sequence = bd[id_seq]
                        
                        print('The most similar sequence in database is %s. \n' %sequence)
                        
                    elif (sm=='no'):
                        match = int(input('Score para match: '))
                        mismatch = int(input('Score para mismatch: '))
                        scores = []
                        id_list = []
                        smatrix = SubstMatrix()
                        s_matrix = smatrix.createFromMatchPars(match, mismatch)
                        for ID in bd:
                            seq2 = bd[ID]
                            alseq = AlignSeq(s_matrix, gap)
                            alseq.needlemanWunsch(seq1, seq2)
                            align = alseq.recoverAlignment()
                            score = 0
                            for i in range(len(align[0])):
                                score += int(alseq.scorePos(align[0][i], align[1][i]))
                            scores.append(score)
                            id_list.append(ID)
                            
                        max_score = max(scores)
                        index = scores.index(max_score)
                        id_seq = id_list[index]
                        sequence = bd[id_seq]
                            
                        print('The most similar sequence in database is %s. \n' %sequence)
                    
                elif (file_type=='txt'):
                    s = functionalities.read_txt(filename)
                    seq1 = s.upper()
                    sm = input('Do you want to use BLOSUM62 as substitution matrix? [yes/no] ')
                    gap = int(input('What is the value for penalty gap? '))
                    if (sm=='yes'):
                        scores = []
                        id_list = []
                        smatrix = SubstMatrix()
                        s_matrix = smatrix.loadFromFile('blosum62.mat','\t')
                        for ID in bd:
                            seq2 = bd[ID]
                            alseq = AlignSeq(s_matrix, gap)
                            alseq.needlemanWunsch(seq1, seq2)
                            align = alseq.recoverAlignment()
                            score = 0
                            for i in range(len(align[0])):
                                score += int(alseq.scorePos(align[0][i], align[1][i]))
                            scores.append(score)
                            id_list.append(ID)
                            
                        max_score = max(scores)
                        index = scores.index(max_score)
                        id_seq = id_list[index]
                        sequence = bd[id_seq]
                        
                        print('The most similar sequence in database is %s. \n' %sequence)
                        
                    elif (sm=='no'):
                        match = int(input('Score para match: '))
                        mismatch = int(input('Score para mismatch: '))
                        scores = []
                        id_list = []
                        smatrix = SubstMatrix()
                        s_matrix = smatrix.createFromMatchPars(match, mismatch)
                        for ID in bd:
                            seq2 = bd[ID]
                            alseq = AlignSeq(s_matrix, gap)
                            alseq.needlemanWunsch(seq1, seq2)
                            align = alseq.recoverAlignment()
                            score = 0
                            for i in range(len(align[0])):
                                score += int(alseq.scorePos(align[0][i], align[1][i]))
                            scores.append(score)
                            id_list.append(ID)
                            
                        max_score = max(scores)
                        index = scores.index(max_score)
                        id_seq = id_list[index]
                        sequence = bd[id_seq]
                            
                        print('The most similar sequence in database is %s. \n' %sequence)
                    
                else:
                    print('The program only supports this 2 types. Please convert your sequence in one of this.')
                    
        
        elif(x=='9'):
            n_sequences = int(input('How many sequences do you want to use to create the Phylogenetic Tree? '))
            if n_sequences > len(bd):
                print('Number of sequences exceeded database size.')
            else:
                gap = int(input('What is the value for penalty gap? '))
                seq_list = []
                i=0
                while i < n_sequences:
                    ID = input('Input sequence  ID: ')
                    if ID in bd.keys():
                        seq_list.append(bd[ID])
                        i+=1
                    else: 
                        print('Invalid ID. Try again.')
                
                sm = input('Do you want to use BLOSUM62 as substitution matrix? [yes/no] ')
                if (sm=='yes'):
                    print()
                    smatrix = SubstMatrix()
                    s_matrix = smatrix.loadFromFile('blosum62.mat','\t')
                    alseq = AlignSeq(s_matrix, gap)
                    up  = UPGMA(seq_list, alseq)
                    arv = up.run()
                    arv.printtree()
                    
                elif (sm=='no'):
                    match = int(input('Score para match: '))
                    mismatch = int(input('Score para mismatch: '))
                    print()
                    smatrix = SubstMatrix()
                    s_matrix = smatrix.createFromMatchPars(match, mismatch)
                    alseq = AlignSeq(s_matrix, gap)
                    up  = UPGMA(seq_list, alseq)
                    arv = up.run()
                    arv.printtree()
                    

        elif(x=='10'):
            n_sequences = int(input('How many sequences do you want to use to create the Phylogenetic Tree? '))
            gap = int(input('What is the value for penalty gap? '))
            seq_list = []
            i=0
            while i < n_sequences:
                ID = input('Input sequence  ID: ')
                if ID in bd.keys():
                    seq_list.append(bd[ID])
                    i+=1
                else: 
                    print('Invalid ID. Try again.')
            
            sm = input('Do you want to use BLOSUM62 as substitution matrix? [yes/no] ')
            if (sm=='yes'):
                print()
                smatrix = SubstMatrix()
                s_matrix = smatrix.loadFromFile('blosum62.mat','\t')
                alseq = AlignSeq(s_matrix, gap)
                up  = UPGMA(seq_list, alseq)
                arv = up.run()
                arv.printtree()
                seqs_ord = arv.getCluster()
                nova_lista = []
                for i in seqs_ord:
                    nova_lista.append(seq_list[i])
                ma = MultipleAlign(seq_list, alseq)
                alinm = ma.alignConsensus()
                print()
                print('Multiple alignment score: %d' %(ma.scoreSP(alinm)))
                print(alinm)

            
            elif (sm=='no'):
                match = int(input('Score para match: '))
                mismatch = int(input('Score para mismatch: '))
                print()
                smatrix = SubstMatrix()
                s_matrix = smatrix.createFromMatchPars(match, mismatch)
                alseq = AlignSeq(s_matrix, gap)
                up  = UPGMA(seq_list, alseq)
                arv = up.run()
                arv.printtree()
                seqs_ord = arv.getCluster()
                nova_lista = []
                for i in seqs_ord:
                    nova_lista.append(seq_list[i])
                ma = MultipleAlign(nova_lista, alseq)
                alinm = ma.alignConsensus()
                print()
                print('Multiple alignment score: %d \n' %(ma.scoreSP(alinm)))
                print(alinm)


        elif(x=='11'):
            s = input('Input protein sequence to search: ')
            seq = s.upper()
            sequence = Seq(seq)
            if sequence.valida():
                handle = SearchDomain(sequence)
                handle.Prosite_Domain()

        elif(x=='12'):
            id_seq = input('Insert the sequence id:')
            form_type = input('Select the type of format (gb or fasta):')
            s=SaveSequence(id_seq)
            if form_type.lower() == 'gb':
                file = s.protein_GenBank()
            elif form_type.lower() == 'fasta':
                file = s.protein_fasta()
            else:
                print('Invalid option!')
            print(file)
        
        elif(x=='13'):
            e = float(input('Insert the trheshold value: '))
            b=Blast(e)
            print('Select the parameters for the blast, if the parameter is not necessary, write None.')
            print(b.make_blast())
        
        elif(x=='14'):
            ans =  input('This functionality requir an input file with sequences. Do you want to continue? [yes/no] ')
            file = input('File name: ')
            if(ans=='yes'):
                try:
                    seqs = AlinhamentoMultiplo(file)
                    seqs.alinhamento_multiplo()
                except:
                    print('It was not possible to answer your request.')
        
        elif(x=='15'):
            ans =  input('This functionality requir an input file with sequences. Do you want to continue? [yes/no] ')
            file = input('File name: ')
            if(ans=='yes'):
                try:
                    align = Filogenia(file)
                    align.ArvoreFilo()
                except:
                    print('It was not possible to answer your request.')
        
        elif(x=='16'):
            ans =  input('This functionality requir an UniProt id. Do you want to continue? [yes/no] ')
            if(ans=='yes'):
                ID = input('Input an UniProt id: ')
                ptm=SearchPTM(ID)
                ptm.Uniprot_records()
                ptm.prot_basic_info()
                more = input('\nDo you want to see more detailed information, like location? [yes/no] ')
                if(more=='yes'):
                    ptm.prot_deep_info()
            print()
            
            

if __name__ == "__main__":
    main()
    