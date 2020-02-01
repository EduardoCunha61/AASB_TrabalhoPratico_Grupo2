#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 29 11:52:45 2020

@author: anabarbosa
"""
#from abc import ABC

# =============================================================================
#   Classe Seq
# =============================================================================

class Seq:

    def __init__(self, sequencia, propriedades=None):
        self.sequencia = sequencia.upper()
        self.propriedades = propriedades

    def __eq__(self,obj):
        if isinstance(obj,Seq):
            return self.sequencia == obj.sequencia
        else:
            return False
    
    def __len__(self):
        return len(self.sequencia)
    
    def __getitem__(self, n):
        return self.sequencia[n]

    def __getslice__(self, i, j):
        return self.sequencia[i:j]

    def __str__(self):
        return self.sequencia
    
    def add_prop(self,key,value):
        self.propriedades[key] = value
    
    def get_prop(self,key):
        if key not in self.propriedades.keys():
            raise IndexError
        return self.propriedades[key]

    def alfabeto(self): 
        return "ACDEFGHIKLMNPQRSTVWY"
    
    def valida(self):
        alf = self.alfabeto()
        res = True
        i = 0
        while i < len(self.sequencia) and res:
            if self.sequencia[i] not in alf: 
                res = False
            else: i += 1
        return res


# =============================================================================
#  Classe MatrixNum
# =============================================================================

class MatrixNum:
    def __init__(self, rows, cols): 
        self.mat = []
        for i in range(rows):
            self.mat.append([]) 
            for j in range(cols):
                self.mat[i].append(0)
        
    def __getitem__ (self,ij):
        i,j = ij
        return self.mat[i][j]
    
    def __setitem__(self,ij,value):
        i,j = ij
        self.mat[i][j] = value
    
    def numRows (self): 
        return len(self.mat)
    
    def numCols (self): 
        return len(self.mat[0])
    
    def getValue (self, i, j): 
        return self.mat[i][j]
    
    def setValue (self, i, j, value): 
        self.mat[i][j] = value
    
    def printmat (self):
        for r in self.mat: 
            print(r) 
        print()

    def addRow(self, newrow):
        self.mat.append(newrow)
        
    def addCol(self, newcol):
        for r in range(self.numRows()):
            self.mat[r].append(newcol[r])
    
    def removeRow(self, ind):
        del self.mat[ind]
        
    def removeCol(self, ind):
        for r in range(self.numRows()):
            del self.mat[r][ind]
            
    def copy(self):
        newm = MatrixNum(self.numRows(), self.numCols())
        for i in range(self.numRows()):
            for j in range(self.numCols()):
                newm.mat[i][j] = self.mat[i][j]
        return newm
        
    def minDistIndexes (self):
        m = self.mat[1][0]
        res= (1,0)
        for i in range(self.numCols()):
            for j in range(i+1, self.numRows()):
                if self.mat[j][i] < m:
                    m = self.mat[j][i]
                    res = (j, i)
        return res

# =============================================================================
#   Classe SubstMatrix
# =============================================================================

class SubstMatrix:

    def __init__(self):
        self.alphabet = ""
        self.sm = {}
        
    def loadFromFile(self, filename, sep):
        f = open(filename, "r")
        line = f.readline()
        tokens = line.split(sep)
        ns = len(tokens)
        self.alphabet = ""
        for i in range(0, ns): 
            self.alphabet += tokens[i][0]
        for i in range(0,ns):
            line = f.readline();
            tokens = line.split(sep);
            for j in range(0, len(tokens)):
                k = self.alphabet[i]+self.alphabet[j]
                self.sm[k] = int(tokens[j])
        f.close()
        return self.sm
                
    def createFromMatchPars(self, match, mismatch):
        alphabet = 'ACDEFGHIKLMNPQRSTVWY'
        for c1 in alphabet:
            for c2 in alphabet:
                if (c1 == c2):
                    self.sm[c1+c2] = match
                else:
                    self.sm[c1+c2] = mismatch
        return self.sm
    
    def scorePair(self, c1, c2):
        if c1 not in self.alphabet or c2 not in self.alphabet:
            return None
        return self.sm[c1+c2]
    
    def __getitem__(self, ij):
        i, j = ij
        return self.scorePair(i, j)
    
# =============================================================================
#   Classe MyAlign
# =============================================================================
    
class MyAlign:

    def __init__(self, lseqs, tipo="protein"):
        self.listseqs = lseqs
        self.tipo = tipo
    
    def __len__(self):
        return len(self.listseqs[0])
    
    def __getitem__(self, n):
        if type(n) is tuple and len(n) ==2: 
            i, j = n
            return self.listseqs[i][j]
        elif type(n) is int: return self.listseqs[n]
        return None
    
    def __str__(self):
        res = ""
        for seq in self.listseqs:
            res += "\n" + seq 
        return res
    
    def numSeqs(self):
        return len(self.listseqs)
    
    def column(self, indice):
        res = []
        for k in range(len(self.listseqs)):
            res.append(self.listseqs[k][indice])
        return res

    def consensus(self):
        cons = ""

        for i in range(len(self)):
            cont = {}
            for k in range(len(self.listseqs)):
                c = self.listseqs[k][i]
                if c in cont:
                    cont[c] = cont[c] + 1
                else:
                    cont[c] = 1
            maximum = 0
            cmax = None
            for ke in cont.keys():
                if ke != "-" and cont[ke] > maximum:
                    maximum = cont[ke]
                    cmax = ke
            cons = cons + cmax
        return cons

# =============================================================================
#   Classe AlignSeq
# =============================================================================

class AlignSeq:

    def __init__(self, sm, g):
        self.g = g
        self.sm = sm
        self.S = None
        self.T = None
        self.seq1 = None
        self.seq2 = None
        
    def scorePos (self, c1, c2):
        if c1 == "-" or c2=="-":
            return self.g
        else:
            return self.sm[c1+c2]
        
    def scoreAlin (self, alin):
        res = 0;
        for i in range(len(alin)):
            res += self.scorePos (alin[0][i], alin[1][i])
        return res
    
    def needlemanWunsch (self, seq1, seq2):
        self.S = [[0]]
        self.T = [[0]]
        self.seq1 = seq1
        self.seq2 = seq2
        for j in range(1, len(seq2)+1):
            self.S[0].append(self.g * j)
            self.T[0].append(3)
        for i in range(1, len(seq1)+1):
            self.S.append([self.g * i])
            self.T.append([2])
        for i in range(0, len(seq1)):
            for j in range(len(seq2)):
                s1 = self.S[i][j] + self.scorePos (seq1[i], seq2[j])
                s2 = self.S[i][j+1] + self.g
                s3 = self.S[i+1][j] + self.g
                self.S[i+1].append(max(s1, s2, s3))
                self.T[i+1].append(max3t(s1, s2, s3))
        return self.S[len(seq1)][len(seq2)]
    
    def recoverAlignment (self):
        res = ["", ""]
        i = len(self.seq1)
        j = len(self.seq2)
        while i>0 or j>0:
            if self.T[i][j]==1:
                res[0] = self.seq1[i-1] + res[0]
                res[1] = self.seq2[j-1] + res[1]
                i -= 1
                j -= 1
            elif self.T[i][j] == 3:
                res[0] = "-" + res[0]
                res[1] = self.seq2[j-1] + res[1] 
                j -= 1
            else:
                res[0] = self.seq1[i-1] + res[0]
                res[1] = "-" + res[1]
                i -= 1
        return MyAlign(res)

def max3t (v1, v2, v3):
    if v1 > v2:
        if v1 > v3: return 1
        else: return 3
    else:
        if v2 > v3: return 2
        else: return 3

def printMat (mat):
    for i in range(0, len(mat)):
        print(mat[i])

# =============================================================================
#   Classe MultipleAlign
# =============================================================================
    
class MultipleAlign(object):

    def __init__(self, seqs, alignseq):
        self.seqs = seqs
        self.alignpars = alignseq
    
    def numSeqs(self):
        return len(self.seqs)

    def addSeqAlignment(self, alignment, seq):
        res = []
        for i in range(len(alignment.listseqs) + 1):
            res.append("")
        cons = Seq(alignment.consensus())
        self.alignpars.needlemanWunsch(cons, seq)
        align2 = self.alignpars.recoverAlignment()
        orig = 0
        for i in range(len(align2)):  # cada coluna
            if align2[0, i] == '-':
                for k in range(len(alignment.listseqs)):
                    res[k] +=    "-"
            else:
                for k in range(len(alignment.listseqs)):
                    res[k] += alignment[k, orig]
                orig += 1
        res[len(alignment.listseqs)] = align2.listseqs[1]
        return MyAlign(res)

    def alignConsensus(self):
        self.alignpars.needlemanWunsch(self.seqs[0], self.seqs[1])
        res = self.alignpars.recoverAlignment()
        for i in range(2, len(self.seqs)):
            res = self.addSeqAlignment(res, self.seqs[i])
        return res

    def scoreColumn (self, charsCol):
        sc = 0;
        for i in range(len(charsCol)):
            for j in range(i+1, len(charsCol)):
                if charsCol[i]!='-' or charsCol[j]!='-':
                    sc += self.alignpars.scorePos(charsCol[i], charsCol[j])
        return sc
    
    def scoreSP (self, alinhamento):
        sp = 0;
        ncols = len(alinhamento[0])
        for j in range(ncols):
            charsCol = alinhamento.column(j)
            scoreCol = self.scoreColumn(charsCol)
            sp += scoreCol
        return sp

    def printMat (self,mat):
        for i in range(0, len(mat)):
            print(mat[i])
    
# =============================================================================
#   Classe BinaryTree
# =============================================================================

class BinaryTree:

    def __init__(self, val, dist=0, left = None, right = None):
        self.value = val
        self.distance = dist
        self.left = left
        self.right = right
    
    def getCluster(self):
        res = []
        if self.value >= 0:
            res.append(self.value)
        else:
            if (self.left != None):
                res.extend(self.left.getCluster()) 
            if (self.right != None):
                res.extend(self.right.getCluster())
        return res

    def printtree(self):
        self.printtreerec(0, "Root")
        
    
    def printtreerec (self, level, side):
        import sys
        for i in range(level): 
            sys.stdout.write("\t")
        al = self.getCluster()
        sys.stdout.write(side + ":" + str(al)+ " Dist.: " + str(self.distance) + "\n")
        if self.value < 0:
            if (self.left != None): 
                self.left.printtreerec(level+1, "Left")
            else: 
                sys.stdout.write("Null")
            if (self.right != None): 
                self.right.printtreerec(level+1, "Right")
            else: 
                sys.stdout.write("Null\n")
        return al

# =============================================================================
#   Classe ClustHier
# =============================================================================

class ClustHier:

    def __init__(self, matdists):
        self.matdists = matdists
        
    def distance (self, tree1, tree2):
        c1 = tree1.getCluster()
        c2 = tree2.getCluster()
        sd = 0.0
        for i in range(len(c1)):
            for j in range(len(c2)):
                sd += self.matdists.getValue(c1[i],c2[j])
        return sd/(len(c1)*len(c2))
    
    def executeClustering(self):
        trees = []
        tableDist = self.matdists.copy()
        for i in range(self.matdists.numRows()):
            t = BinaryTree(i)
            trees.append(t)
        for k in range(self.matdists.numRows(), 1, -1):
            mins = tableDist.minDistIndexes()
            i = mins[0]
            j = mins[1]
            n = BinaryTree(-1, tableDist.getValue(i, j)/2.0, trees[i], trees[j])
            if k>2:
                trees.pop(i)
                trees.pop(j)
                tableDist.removeRow(i)
                tableDist.removeRow(j)
                tableDist.removeCol(i)
                tableDist.removeCol(j)
                dists = []
                for x in range(len(trees)):
                    dists.append(self.distance(n, trees[x]))
                tableDist.addRow(dists)
                cdists = []
                for y in range(len(dists)):
                    cdists.append(dists[y])
                cdists.append(0.0)
                tableDist.addCol(cdists)
                trees.append(n)
            else: 
                return n

# =============================================================================
#   Classe UPGMA
# =============================================================================
  
class UPGMA:

    def __init__(self, seqs, alseq):
        self.seqs = seqs
        self.alseq = alseq
        self.matdist = MatrixNum(len(seqs), len(seqs))
        self.criaMatDists()
        
    def criaMatDists(self):
        for i in range(len(self.seqs)): 
            self.matdist.setValue(i, i, 0.0)
        
        for i in range(len(self.seqs)):
            for j in range(i, len(self.seqs)):
                s1 = self.seqs[i]
                s2 = self.seqs[j]
                self.alseq.needlemanWunsch(s1, s2) 
                alin = self.alseq.recoverAlignment() 
                ncd = 0
                for k in range(len(alin)):
                    col = alin.column(k)
                    if (col[0] != col[1]): 
                        ncd += 1
                self.matdist.setValue(i, j, ncd) 
                self.matdist.setValue(j, i, ncd)
                
    def run(self):
        ch = ClustHier(self.matdist)
        arv= ch.executeClustering()
        return arv

# =============================================================================
#   Classe SearchDomain
# =============================================================================


class SearchDomain:
    
    def __init__(self,seq_input):
        self.__seq_input = seq_input
        
    def Prosite_Domain(self):
        from Bio import ExPASy
        from Bio.ExPASy import Prosite, ScanProsite
        try:
            handle = ScanProsite.scan(seq=self.__seq_input)
            result = ScanProsite.read(handle)
            if len(result) != 0:
                for res in range(len(result)):
                    prosite_acession = result[res]['signature_ac']
                    r = ExPASy.get_prosite_raw(prosite_acession)
                    html = Prosite.read(r)
                    r.close()
                    print('Foi encontrado um dominio %s.' %(html.name))            
            else:
                print('Não foram encontradas correspondências.')
        except:
            print('A sequência fornecida não é uma sequência proteica.')

# =============================================================================
#   Classe SaveSequence
# =============================================================================


class SaveSequence:
    
    def __init__(self, id_seq=None):
        self.__id_seq = id_seq
        
    def protein_fasta(self):
        
        from Bio import Entrez, SeqIO
        Entrez.email = 'example@gmail.com'
        
        protein_fasta = Entrez.efetch(db="protein", id=self.__id_seq, retmode='text', rettype='fasta')
        protein_read_fasta = SeqIO.read(protein_fasta,'fasta')
        protein_fasta.close()
        print(protein_read_fasta.seq)
        SeqIO.write(protein_read_fasta,'protein_fasta.txt','fasta')
        print('The fasta file was successfully saved: sequence.txt')
        
    def protein_GenBank(self):
        
        from Bio import Entrez, SeqIO
        Entrez.email = 'example@gmail.com'

        gene_gb = Entrez.efetch(db="protein", id=self.__id_seq, retmode='text', rettype='gb')
        gene_read_gb = SeqIO.read(gene_gb,'gb')
        gene_gb.close()
        SeqIO.write(gene_read_gb,'protein_GenBank.txt','gb')
        print('The Genbank file was successfully saved: sequence_GenBank.txt')


# =============================================================================
#   Classe SearchPTM
# =============================================================================

class SearchPTM:
    
    def __init__(self,uniprot_id):
        self.__uniprot_id = uniprot_id
    
    def Uniprot_records(self):
        
        import requests
        from Bio import Entrez, ExPASy, Seq, SeqIO, SeqRecord, SwissProt
        from urllib.request import urlopen   
        
        handle = ExPASy.get_sprot_raw(self.__uniprot_id)#ID do NCBI, para tirar ficheiro xml da Uniprot
        url = handle.url 
        url = url.replace('txt','xml') 
        response = requests.get(url) 
        with open('Uniprot' + self.__uniprot_id + '.xml','wb') as file: #b para escrever em modo binário
            file.write(response.content)
        
    def prot_basic_info(self):
        
        import requests
        from Bio import Entrez, ExPASy, Seq, SeqIO, SeqRecord, SwissProt
        from urllib.request import urlopen   

        handle_uniprot = SeqIO.read('Uniprot' + self.__uniprot_id + '.xml','uniprot-xml') #Objeto SeqRecord com informações do xml da Uniprot
        try:
            dic = {}
            for feat in handle_uniprot.features:
                if  feat.type in dic.keys():
                    dic[feat.type] += 1
                else:
                    dic[feat.type] = 1
                    
            for key in dic.keys():
                print('Existem %d %s.' %(dic[key],key))

        except:
            print('Sem informação correspondente a PTM')

    def prot_deep_info(self):

        from Bio import Seq, SeqIO, SeqRecord
        
        handle_uniprot = SeqIO.read('Uniprot' + self.__uniprot_id + '.xml','uniprot-xml') #Objeto SeqRecord com informações do xml da Uniprot
        try:
            for feat in handle_uniprot.features:
                print('\nFeature type: ' + feat.type)
                print('Location: ' + str(feat.location))
        except:
            print('Sem informação correspondente a PTM')


# =============================================================================
#   Classe AlinhamentoMultiplo
# =============================================================================


class AlinhamentoMultiplo:

    def __init__(self, seq_input):
        self.seq_input = seq_input
        self.diretoria = r'C:\Program Files (x86)\ClustalW2\clustalw2.exe'

    def alinhamento_multiplo(self):
        from Bio.Align.Applications import ClustalwCommandline
        clustalw_cline = ClustalwCommandline(self.diretoria, infile = self.seq_input)
        clustalw_cline()

    def get_diretoria(self): 
        return self.diretoria
    
    def set_diretoria(self,nova_diretoria):
        self.diretoria = nova_diretoria

# =============================================================================
#   Classe Filogenia
# =============================================================================

class Filogenia:
    
    def __init__(self, input_file):
        self.__input_file = input_file
        
    def ArvoreFilo(self):
        from Bio import Phylo
        tree = Phylo.read(self.__input_file, 'newick')
        Phylo.draw_ascii(tree)

# =============================================================================
#   Classe Blast
# =============================================================================


class Blast:

    def __init__(self, e_value, database="nr"):
        self.__database = database
        self.__dic_blast = {}
        self.__e_value = e_value

    def make_blast(self):
        from Bio.Blast import NCBIWWW, NCBIXML
        self.__dic_blast['program']=str(input('Blast progam (blastp or tblastn):'))
        self.__dic_blast['sequence']=str(input('Insert the sequence:'))
        self.__dic_blast['entrez_query']=str(input('This restricts the search(ex:Homo sapiens [organism]):'))
        self.__dic_blast['hitlist_size']=int(input('Number of databases sequences to keep:'))
        result_handle = NCBIWWW.qblast(program = self.__dic_blast['program'], database = self.__database, sequence = self.__dic_blast['sequence'], entrez_query=self.__dic_blast['entrez_query'], hitlist_size=self.__dic_blast['hitlist_size'])
        blast_records = NCBIXML.parse(result_handle)
        seqs_result = []
        for blast_record in blast_records:
            for alignment in blast_record.alignments:
                for hsp in alignment.hsps:
                    if hsp.expect < self.__e_value:
                        seqs_result.append([alignment.title, hsp.expect, hsp.sbjct.replace('-', '')])
        return seqs_result





