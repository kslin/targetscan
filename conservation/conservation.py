import pandas as pd
import pomegranate as pom

model = pom.HiddenMarkovModel().from_json('parameters1000.txt')

UTRs = pd.read_csv('~/Downloads/tiny_utrs.txt', sep='\t', header=None)


