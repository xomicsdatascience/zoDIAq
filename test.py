import pandas as pd

df = pd.read_csv('Data/MSPLIToutput.txt', sep='\t')
df = df.drop(columns='Name')
df.to_csv('Data/MSPLIToutput.csv', index=False)
