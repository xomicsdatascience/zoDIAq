import pandas as pd
import re
from Bio import SeqIO
import pickle

#this process takes a while (~10-15 minutes), so I just pickled it immediately after.
'''
input_file = "/Users/calebcranney/Desktop/Meyer Lab Project/Jesse_Github/DI2A/outputs/Python/peptide_target_table_unfiltered.txt"
fasta_file = "/Users/calebcranney/Desktop/Meyer Lab Project/Jesse_Github/DI2A/Data/targeted_quant/2019-03-14-td-UP000005640.fasta"
df = pd.read_csv(input_file, sep="\t").sort_values('cosine', ascending=False)

#12
df.Peptide = [re.sub("\+\d+\.\d+", "", x) for x in df.Peptide.tolist()]

#13
matching_proteins = [[] for i in df.Peptide]
count = 0
match_count = 0
for record in SeqIO.parse(fasta_file, "fasta"):
    count += 1
    if count % 100 == 0:
        print(count)
        print(match_count)
    for i in range(len(df.Peptide)):
        if df.Peptide.iloc[i] in record.seq:
            matching_proteins[i].append(record.id)
            match_count += 1

#15
df['Protein'] = matching_proteins

df.to_pickle("./protein_table.pkl")
'''
df = pd.read_pickle("./protein_table.pkl")

test = pd.read_csv("/Users/calebcranney/Desktop/Meyer Lab Project/Jesse_Github/DI2A/outputs/num_prot_matches.csv")


#code where I determined which loci are different between Jesse's R code and my python code.
'''
print(len(df.Protein))
print(len(test))

fails = []
for i in range(len(df.Protein)):
    if len(df.Protein.loc[i]) != test.x.loc[i]:
        fails.append([i,len(df.Protein.loc[i]),str(test.x.loc[i])])
        #print(str(i) + ", df= " + str(len(df.Protein.loc[i])) + ", test= "+ str(test.x.loc[i]))
fails = pd.DataFrame(fails, columns=['index','pythonScript','RScript'])
fails.to_csv("./fails.csv")
'''

#16 - Note that this includes those with length 0
df = df[ df['Protein'].map(len) < 2 ]


#18
df = df.reset_index(drop=True)
for i in range(len(df)):
    if df.Name.loc[i] == 'DECOY_null':
        df.Protein.loc[i] = ['DECOY_null'+str(i)]


#I added this myself - to be discussed with Jesse, but I'm removing all non-decoy lines that did not have a protein match
remove = []
for i in range(len(df.Protein)):
    if len(df.Protein.loc[i])==0:
        remove.append(i)
df = df.take(list(set(range(df.shape[0])) - set(remove)))

df['Protein'] = [ x[0] for x in df.Protein.tolist() ]

#19
df = df.drop_duplicates(subset="Protein")

#20 - I think there are two parts to this, the first not used, so I'm doing this based on my understanding
# Review with Jesse
df = df.reset_index(drop=True)
decoy_indices = []
for i in range(len(df)):
    if df.Name.loc[i] == 'DECOY_null':
        decoy_indices.append(i)

def FDRLastDecoyIndex(FDR, df, decoy_indices):
    for i in range(len(decoy_indices)):
        curr_FDR = (i+1)/(decoy_indices[i]+1)
        if curr_FDR > FDR:
            return i
    return -1

last_decoy = FDRLastDecoyIndex(0.025,df,decoy_indices)
if last_decoy == 0:
    df = pd.DataFrame()
elif last_decoy > 0:
    df = df.truncate(after=last_decoy-1)

#22
if last_decoy != -1:
    decoy_indices = decoy_indices[:last_decoy]
df = df.take(list(set(range(df.shape[0])) - set(decoy_indices)))


#24
df = df.sort_values(by=['CV','cosine'], ascending = False)
df = df.reset_index(drop=True)

df = df[['Peptide','Mz.1','Heavymz','z.1','Scan#','ylight_quant_frags','yheavy_quant_frags','fragment_ordinals','CV']]
new_columns = ['Peptide','prec_light_mz','prec_heavy_mz','z','scan','ylight','yheavy','yordinals','CV']
df.columns = new_columns
df.scan = list(range(1,len(df)+1))
