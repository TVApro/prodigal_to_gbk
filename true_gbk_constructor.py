import os
import sys
import re
from re import sub
import textwrap as tw

print('@@@ START GENERATING OF GBK FILE @@@')
def killsymbol_last(X, string):
    string_lenth = len(string)
    X_lenth = len(X)
    if string[string_lenth-X_lenth:string_lenth] == X:
        string = string[:string_lenth-X_lenth]
    return string

def addsymbol_last(X, string):
    string_lenth = len(string)
    X_lenth = len(X)
    if string[string_lenth-X_lenth:string_lenth] != X:
        string = string+X
    return string

def killsymbol_first(X, string):
    X_lenth = len(X)
    if string[0:X_lenth] == X:
        string = string[X_lenth:]
    return string

def addsymbol_first(X, string):
    X_lenth = len(X)
    if string[0:X_lenth] != X:
        string = X+string
    return string

def printlist(name, x):
    global string_X
    for i in x:
        string_X += str(i)
        string_X += '\n'
    with open(name, 'w') as myfile5:
        print(string_X, file=myfile5)

def fasta_obtain(text):
    def fasta_add_spaces(fasta_seq):
        def convert_tuple(x):
            a = str('')
            for i in x:
                a += i+' '
            return a
        result = tw.wrap(fasta_seq, 60)
        number = 1
        true_fasta = str("")
        for i in result:
            p = tw.wrap(i, 10)
            space = str((9 - len(str(number)))*' ')
            fa = convert_tuple(p)
            true_fasta += space + str(number) + ' ' + fa + '\n'
            number += len(i)
        return true_fasta
    
    fasta_true_list = {}
    fasta_list = text.split(">")
    for i in fasta_list:
            if i:
                end = i.find("\n")
                fasta_name = '>'+sub('\n', '', i[:end])
                fasta_sequence = sub('\n', '', i[end:])
                X = fasta_add_spaces(fasta_sequence)
                fasta_true_list[fasta_name]=X
    return fasta_true_list

def faa_to_gbk(gbk, faa):
    def faa_obtain(text):
        fasta_faa_list = {}
        fasta_list = text.split(">")
        for i in fasta_list:
            if i:
                end = i.find("\n")
                fasta_name = '>'+sub('\n', '', i[:end])
                fasta_sequence = sub('\n', '', i[end:])
                fasta_faa_list[fasta_name]=fasta_sequence
        return fasta_faa_list

    def faa_work_with_name(name):
        faa_rest = re.sub('>', '', name)
        faa_microlist = faa_rest.split('#')
        for i in faa_microlist:
            if 'NODE_' in i:
                contig = sub(' ', '', i)
            if 'ID=' in i:
                protein_ID = sub(' ', '', i)
        name_info = [contig, protein_ID]
        return name_info

    faa_list = faa_obtain(faa)
    gbk_list = gbk.split('//')
    for g, i in enumerate(gbk_list):
        ind_DEFINITION = i.find("DEFINITION")
        ind_FEATURES = i.find("FEATURES")
        ind_ORIGIN = i.find("ORIGIN")
        definition = i[ind_DEFINITION:ind_FEATURES]
        contig_find = re.findall(r';seqhdr="NODE.+?";', definition)
        contig_orig = str(contig_find).replace('"', '').replace("[';seqhdr=", '').replace(";']", '')
        features = i[ind_FEATURES:ind_ORIGIN]
        CDSs = features.split('     CDS             ')
        for index, CDS in enumerate(CDSs):
            if not 'FEATURES' in CDS:
                if CDS:
                    newCDS = str('')
                    CDS = '     CDS             '+CDS
                    CDS = addsymbol_last('\n', CDS)
                    for faa in faa_list:
                        name_info = faa_work_with_name(faa)
                        contig = name_info[0]
                        protein_ID = name_info[1]
                        if contig_orig in contig:
                            if protein_ID in CDS:
                                newCDS += CDS+'                     /translation="'+faa_list[faa]+'"'+'\n'
                                CDSs[index] = newCDS
        CDS_string = str('')
        for j in CDSs:
            CDS_string += j
        i_new = i.replace(features, CDS_string)
        linecount_i = len(i_new)
        i_new = killsymbol_last('\n', i_new)
        gbk_list[g] = i_new
    return gbk_list

path = os.getcwd()
fasta_name = sys.argv[1]
faa_name = sys.argv[2]
gbk_name = sys.argv[3]

fasta_file = os.path.join(path, fasta_name)
faa_file = os.path.join(path, faa_name)
gbk_file = os.path.join(path, gbk_name)

with open(fasta_file, 'r') as myfile0:
    fasta = myfile0.read()
with open(gbk_file, 'r') as myfile:
    gbk = myfile.read()
with open(faa_file, 'r') as myfile3:
    faa = myfile3.read()
fasta_true_list = fasta_obtain(fasta)

gbk_list = faa_to_gbk(gbk, faa)
for i in gbk_list:
    if not i:
         gbk_list.remove(i)
for index, item in enumerate(gbk_list):
    for i in fasta_true_list:
        without = sub('>', '', i)
        if without in item:
            inew = "LOCUS       "+without+addsymbol_first('\n', item)+'ORIGIN'+'\n'+fasta_true_list[i]+'//'
            gbk_list[index] = inew
string_X = ''  
printlist(fasta_name+'.gbk', gbk_list)
print('@@@ FINISH GENERATING OF GBK FILE @@@')
