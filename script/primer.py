import re
def revers(seq):
    codes={"A":"A", "C":"C", "G":"G", "T":"T", 
           "R":"[AG]", "S":"[GC]", "B":"[CGT]", "Y":"[CT]", 
           "W":"[AT]", "D":"[AGT]", "K":"[GT]", "N":"[ACGT]", 
           "H":"[ACT]", "M":"[AC]", "V":"[ACG]", "X": "[ACGT]"}
    
    d={'A':'T','G':'C','C':'G','T':'A',
      'R':'Y','S':'S','B':'V','Y':'R',
       'W':'W','D':'H','K':'M','N':'N',
       'H':'D','M':'K','V':'B','X':'X'
      }
    new=''
    for i in seq[::-1]:
        new+=d[i]
    return new

def replace_ambiguity_codes(sequence):
    codes={"A":"A", "C":"C", "G":"G", "T":"T", 
           "R":"[AG]", "S":"[GC]","B":"[CGT]", "Y":"[CT]", "W":"[AT]",
           "D":"[AGT]", "K":"[GT]", "N":"[ACGT]", "H":"[ACT]", "M":"[AC]", "V":"[ACG]", "X": "[ACGT]"}
    regexlist=[]
    for base in sequence:
        if not base in codes:
            print("Unidentified nucleotide code in primer:", base)
            sys.exit()
        else:
            regexlist.append(codes[base])

    return ''.join(regexlist)


'''
In solico PCR
def get_target_region(input_fasta,f_primer,r_primer):
    f_primer=replace_ambiguity_codes(f_primer)
    r_primer=replace_ambiguity_codes(r_primer)
    f=open(input_fasta,'r')
    fasta_data=f.read().replace('\r','').split('>')[1:]
    f.close()
    d1={} # d1: key: ID  value: ATCG
    target={}
    for each_seq in fasta_data:
        d1[each_seq.split(' ')[0]]=each_seq.split('\n',1)[1].replace('\n','')
    for each_seq in d1:
        findall = re.findall(f_primer, d1[each_seq])
        if len(findall) == 1:
            start=d1[each_seq].split(findall[0])[1]
            findall2 = re.findall(r_primer,start)
            if len(findall2)==1:
                target[each_seq]=start.split(findall2[0])[0]
    return len(target),len(d1),len(target)/len(d1),target
'''




def cov(f,r,d):
    result=[]
    for i in d:
        hit=[]
        
        for j in d[i]:
            forward=re.findall(replace_ambiguity_codes(f),j)
            if len(forward) ==1:
                reverse=re.findall(replace_ambiguity_codes(revers(r)),j.split(forward[0])[1])
                if len(reverse) ==1:
                    hit.append(j)
        result.append(str((len(hit)/len(d[i]))*100)+'%')
        #print(str((len(hit)/len(d[i]))*100)+'%')
    return result
        #print(len(d[i]))
        #print(i)
        
def cov_e1(f,r,d):
    result=[]
    for i in d:
        hit=[]
        
        for j in d[i]:
            forward=regex.findall('('+replace_ambiguity_codes(f)+'){e<=1}',j)
            if len(forward) ==1:
                reverse=regex.findall('('+replace_ambiguity_codes(revers(r))+'){e<=1}',j.split(forward[0])[1])
                if len(reverse) ==1:
                    hit.append(j)
        result.append(str((len(hit)/len(d[i]))*100)+'%')
        #print(str((len(hit)/len(d[i]))*100)+'%')
    return result
        #print(len(d[i]))
        #print(i)                    
