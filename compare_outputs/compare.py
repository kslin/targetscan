import sys
import pandas as pd

myfile = pd.read_csv('small_output.txt',sep='\t').astype(str).drop('Species ID',1)
myfile = myfile[['Gene ID','Mirbase ID','Site type','Site start','Site end','Threep score','Local AU score','Min dist','UTR length score','ORF length','ORF 8mers','Off6m score','TA','SPS']]
actual = pd.read_csv('Targets.BL_PCT.context_scores_final.txt',sep='\t').astype(str)

actual = actual[actual['Species ID'] == '9606'].drop([
	'Site type contribution','Species ID','sRNA1A contribution','sRNA1C contribution','sRNA1G contribution','sRNA8A contribution','sRNA8C contribution',
	'sRNA8G contribution','site8A contribution','site8C contribution','site8G contribution','context++ score','context++ score percentile',
	'AIR','weighted context++ score','weighted context++ score percentile','UTR region','UTR-miRNA pairing','Group #',
	'mature miRNA sequence','miRNA family'],1)

vikram_params = pd.read_csv('Agarwal_2015_parameters.txt',sep='\t').fillna(1)
vikram_params = vikram_params.set_index('Feature')
vikram_params = vikram_params.astype(float)
# print vikram_params

# actual.columns = ['Gene ID', 'Mirbase ID', 'Site Type', 'UTR start', 'UTR end',
#        '3' pairing contribution', 'local AU contribution',
#        'Min_dist contribution', '3'UTR length contribution',
#        'SA contribution', 'ORF length contribution',
#        'ORF 8mer contribution', 'Offset 6mer contribution',
#        'TA contribution', 'SPS contribution']

actual.columns = [
'Gene ID','Mirbase ID','Site Type','UTR start','UTR end',
'3P_score', 'Local_AU','Min_dist', 'Len_3UTR', 
'SA', 'Len_ORF', 'ORF8m', 'Off6m', 'TA_3UTR', 'SPS','PCT']


for row in vikram_params.iterrows():
	feature = row[0]
	if feature in actual.columns:
		if feature == 'ORF8m':
			params = {}
			params['8mer-1a'] = [row[1]['8mer coeff'],row[1]['8mer min'],row[1]['8mer max']]
			params['7mer-m8'] = [row[1]['7mer-m8 coeff'],row[1]['7mer-m8 min'],row[1]['7mer-m8 max']]
			params['7mer-1a'] = [row[1]['7mer-A1 coeff'],row[1]['7mer-A1 min'],row[1]['7mer-A1 max']]
			params['6mer'] = [row[1]['6mer coeff'],row[1]['6mer min'],row[1]['6mer max']]
			print params
			zipped = zip(actual[feature],actual['Site Type'])
			actual[feature] = [(float(x)/params[stype][0])*(params[stype][2] - params[stype][1]) + params[stype][1] for x,stype in zipped]

myfile = myfile.sort_values(by=['Gene ID','Mirbase ID'])
actual = actual.sort_values(by=['Gene ID','Mirbase ID'])

print list(actual['ORF8m'])
sys.exit()
# for row in myfile.iterrows():
# 	if float(row[1]['ORF 8mers']) != 0:
# 		print row

for row in actual.iterrows():
	if float(row[1]['ORF8m']) != 0:
		print row

myfile.to_csv('compare_myfile.txt',sep='\t',index=False)
actual.to_csv('compare_actual.txt',sep='\t',index=False)
