#!/usr/bin/env python
# encoding: utf-8

# POLYMORPH v.2
# Julie BOGOIN
# June 2021

import subprocess, pandas, vcf, io, numpy
import xlsxwriter
import os
from pandasql import sqldf
import datetime
from db_location import *
import functions

##### MAIN #####
################

print('\n*************************')
print('***** POLYMPORPH v2 *****')
print('*************************\n')

pathfolder = './input'

for fi in os.listdir(pathfolder):

	fullsplit = fi.split('.')
	filename = fullsplit[0]

	print('Fichier en cours de traitement: {}'.format(fi))
	print('Version du genome: hg19')
	print("ok\n")

	print("--> Nettoyage du fichier d'entrée.")

	primers_clean = functions.win2unix('./input/'+ fi)

	df = pandas.read_csv(primers_clean, sep='\t', \
		names=['name_F', 'seq_F', 'name_R', 'seq_R'])
	print('ok')

	print("\n--> Ajout des noms de gènes.")

	gene_name_list = []
	indexNames = []
	index = 0

	for x in df.itertuples():
		namef = x.name_F
		
		if (namef.startswith('#')):
			gene_name = namef.replace("#","")
			gene_name_list.append(gene_name)
			indexNames.append(index)
		
		else:
			gene_name_list.append(gene_name)
			
		index = index + 1

	df['gene'] = pandas.Series(gene_name_list, dtype=str)
	print('ok')

	print("\n--> Suppression des lignes #gènes.")

	# Remove #gene line
	df.drop(indexNames, inplace=True)
	print('ok')

	print("\n--> Vérification des caractères A,T,C,G.")

	size_F = []
	size_R = []

	for index, value in df['seq_F'].items():
	    functions.containsOnlyATGC(value)
	    size_F.append(len(value))

	for index, value in df['seq_R'].items():
	    functions.containsOnlyATGC(value)
	    size_R.append(len(value))

	df['size_F'] = pandas.Series(size_F, dtype='int')
	df['size_R'] = pandas.Series(size_R, dtype='int')

	df.reset_index(inplace=True)

	# columns order
	cols = ['gene','name_F', 'seq_F', 'size_F', 'name_R', 'seq_R', 'size_R']
	df = df[cols]

	df.to_csv(('data/primers_final.csv'), index=False)
	print('ok')

	print("\n--> Création input In-Silico PCR.")

	ispcr_input = df.copy(deep=True)
	del ispcr_input['gene']
	del ispcr_input['name_R']
	del ispcr_input['size_F']
	del ispcr_input['size_R']

	ispcr_input.to_csv('data/ispcr_input.txt', header=None, index=None, sep='\t', mode='a')
	print('ok\n')

	print("--> Dataframe noms de gènes.")

	uniq = pandas.read_csv(guniq, index_col=None, header=None, sep='\t')

	uniq.columns = ['contig', 'strand','gene']

	#Supprimer le 'chr' de la colonne contig
	uniq['contig'] = uniq['contig'].str[3:]
	print('ok')

	print("\n--> Présence des gènes dans la liste hg19?")

	genes_list = uniq['gene'].tolist()

	for x in df.itertuples():
		gene = x.gene
		
		if gene not in genes_list:
			print("Erreur: {} n'est pas dans la liste de gènes.".format(gene))
	print('ok')

	print("\n--> In-Silico PCR.")

	subprocess.call("isPcr +seq+ ./data/ispcr_input.txt \
        -out=bed ./data/ispcr_output.bed", shell="/bin/bash")
	print('ok')

	print("\n--> Creation fichier targets.")

	ispcr_output = pandas.read_csv('data/ispcr_output.bed', sep='\t',\
	    names=['contig', 'amplicon_start', 'amplicon_end', 'name_F', 'score', 'strand'])

	targets_tot = pandas.merge(df, ispcr_output, \
	    how='left', left_on='name_F', right_on='name_F')

	F_start = []
	F_end = []
	R_start = []
	R_end = []

	for x in targets_tot.itertuples():
		if x.strand == "+":
			F_start.append(x.amplicon_start)
			F_end.append(x.amplicon_start + x.size_F - 1)
			R_start.append(x.amplicon_end - x.size_R + 1)
			R_end.append(x.amplicon_end)
		
		if x.strand == "-":
			F_start.append(x.amplicon_end - x.size_F + 1)
			F_end.append(x.amplicon_end)
			R_start.append(x.amplicon_start)
			R_end.append(x.amplicon_start + x.size_R- 1)

	targets_tot['F_start'] = pandas.Series(F_start, dtype='float')
	targets_tot['F_end'] = pandas.Series(F_end, dtype='float')
	targets_tot['R_start'] = pandas.Series(R_start, dtype='float')
	targets_tot['R_end'] = pandas.Series(R_end, dtype='float')

	targets_list = []

	targets_F = targets_tot.copy(deep=True)
	del targets_F['name_R']
	del targets_F['seq_R']
	del targets_F['size_R']
	del targets_F['R_start']
	del targets_F['R_end']
	cols = ['contig', 'F_start', 'F_end', 'name_F', 'score', 'strand', 'gene',\
	    'amplicon_start', 'amplicon_end' ]
	targets_F = targets_F[cols]
	targets_F.rename(columns={'F_start': 'start'}, inplace=True)
	targets_F.rename(columns={'F_end': 'end'}, inplace=True)
	targets_F.rename(columns={'name_F': 'name'}, inplace=True)
	targets_list.append(targets_F)

	targets_R = targets_tot.copy(deep=True)
	del targets_R['name_F']
	del targets_R['seq_F']
	del targets_R['size_F']
	del targets_R['F_start']
	del targets_R['F_end']
	cols = ['contig', 'R_start', 'R_end', 'name_R', 'score', 'strand', 'gene',\
	    'amplicon_start', 'amplicon_end' ]
	targets_R = targets_R[cols]
	targets_R.rename(columns={'R_start': 'start'}, inplace=True)
	targets_R.rename(columns={'R_end': 'end'}, inplace=True)
	targets_R.rename(columns={'name_R': 'name'}, inplace=True)
	targets_list.append(targets_R)

	targets_final = pandas.concat(targets_list, ignore_index=True)

	del targets_final['gene']
	del targets_final['amplicon_start']
	del targets_final['amplicon_end']

	targets_final = targets_final.dropna()

	targets_final['start'] = targets_final['start'].astype(int)
	targets_final['end'] = targets_final['end'].astype(int)
	targets_final['score'] = targets_final['score'].astype(int)

	targets_final.to_csv('data/targets.bed', sep='\t', header=None, index=False)
	print('ok')

	print('\n--> gnomAD')
	print(('\n*** Exomes ***'))

	gnomad_exomes = vcf.Reader(filename=gnex, compressed=True)

	vcf_writer = vcf.Writer(open('gnomad_exomes_fetch.vcf', 'w'), gnomad_exomes, lineterminator='\n')
	
	for x in targets_final.itertuples():
		tig = x.contig
		chromosome = tig.replace('chr','')
		fetch = gnomad_exomes.fetch(chromosome, x.start - 1, x.end + 1)
		
		for record in fetch:
			vcf_writer.write_record(record)

	variants_exomes = functions.read_vcf('gnomad_exomes_fetch.vcf')

	info_split = variants_exomes['INFO'].str.split(pat=';', n=3, expand=True)
	variants_exomes['AC_e'] = info_split[0].str.replace('AC=','')
	variants_exomes['AN_e'] = info_split[1].str.replace('AN=','')
	variants_exomes['AF_e'] = info_split[2].str.replace('AF=','')

	del variants_exomes['INFO']
	del variants_exomes['QUAL']
	del variants_exomes['FILTER']
	print('ok')

	print(('\n*** Genomes ***'))

	gnomad_genomes = vcf.Reader(filename=gnge, compressed=True)

	vcf_writer = vcf.Writer(open('gnomad_genomes_fetch.vcf', 'w'), gnomad_genomes, lineterminator='\n')
	
	for x in targets_final.itertuples():
		tig = x.contig
		chromosome = tig.replace('chr','')
		fetch = gnomad_genomes.fetch(chromosome, x.start - 1, x.end + 1)

		for record in fetch:
			vcf_writer.write_record(record)

	variants_genomes = functions.read_vcf('gnomad_genomes_fetch.vcf')

	info_split = variants_genomes['INFO'].str.split(pat=';', n=3, expand=True)
	variants_genomes['AC_g'] = info_split[0].str.replace('AC=','')
	variants_genomes['AN_g'] = info_split[1].str.replace('AN=','')
	variants_genomes['AF_g'] = info_split[2].str.replace('AF=','')

	del variants_genomes['INFO']
	del variants_genomes['QUAL']
	del variants_genomes['FILTER']
	print('ok')

	print(('\n*** Merging ***'))

	variants_genomes['CHROM'] = variants_genomes['CHROM'].astype(str)
	variants_genomes['POS'] = variants_genomes['POS'].astype(str)
	variants_genomes['REF'] = variants_genomes['REF'].astype(str)
	variants_genomes['ALT'] = variants_genomes['ALT'].astype(str)

	variants_genomes['pos-id'] = variants_genomes['CHROM'] + '-' \
	    + variants_genomes['POS'] + '-' + variants_genomes['REF'] + '-' \
	    + variants_genomes['ALT']

	variants_exomes['CHROM'] = variants_exomes['CHROM'].astype(str)
	variants_exomes['POS'] = variants_exomes['POS'].astype(str)
	variants_exomes['REF'] = variants_exomes['REF'].astype(str)
	variants_exomes['ALT'] = variants_exomes['ALT'].astype(str)

	variants_exomes['pos-id'] = variants_exomes['CHROM'] + '-' \
	    + variants_exomes['POS'] + '-' + variants_exomes['REF'] + '-' \
	    + variants_exomes['ALT']

	gnomad_final = pandas.merge(variants_genomes, variants_exomes, \
	    on='pos-id', how='outer', sort=True)
	print('ok')

	print('\n*** Suppression des duplicats ***')

	gnomad_final.drop_duplicates(subset=None, keep='first', inplace=True, \
	    ignore_index=False)
	print('ok')

	print(('\n*** Calcul des frequences ***'))

	gnomad_final['AC_e'] = gnomad_final['AC_e'].astype(float)
	gnomad_final['AC_g'] = gnomad_final['AC_g'].astype(float)
	gnomad_final['AN_e'] = gnomad_final['AN_e'].astype(float)
	gnomad_final['AN_g'] = gnomad_final['AN_g'].astype(float)
	
	gnomad_final['Freq_gnomAD'] = (gnomad_final['AC_e'] + gnomad_final['AC_g']) / \
	    (gnomad_final['AN_e'] + gnomad_final['AN_g'])

	gnomad_final.loc[gnomad_final['AC_e'].isna(), 'Freq_gnomAD'] = \
	    gnomad_final['AF_g']

	gnomad_final.loc[gnomad_final['AC_g'].isna(), 'Freq_gnomAD'] = \
	    gnomad_final['AF_e']

	gnomad_final.loc[gnomad_final['CHROM_x'].isna(), 'CHROM_x'] = \
	    gnomad_final['CHROM_y']  

	gnomad_final.loc[gnomad_final['POS_x'].isna(), 'POS_x'] = \
	    gnomad_final['POS_y']

	gnomad_final.loc[gnomad_final['ID_x'].isna(), 'ID_x'] = \
	    gnomad_final['ID_y']

	gnomad_final.loc[gnomad_final['REF_x'].isna(), 'REF_x'] = \
	    gnomad_final['REF_y']   

	gnomad_final.loc[gnomad_final['ALT_x'].isna(), 'ALT_x'] = \
	    gnomad_final['ALT_y']

	del gnomad_final['pos-id']
	del gnomad_final['CHROM_y']
	del gnomad_final['ID_y']
	del gnomad_final['REF_y']
	del gnomad_final['ALT_y']

	gnomad_final.to_csv('./data/gnomad_final.csv', index=False)
	print('ok')

	print('\n--> Fichier final.')

	del gnomad_final['AC_e']
	del gnomad_final['AN_e']
	del gnomad_final['AF_e']
	del gnomad_final['AC_g']
	del gnomad_final['AN_g']
	del gnomad_final['AF_g']
	del gnomad_final['POS_y']

	# Make your pysqldf object:
	pysqldf = lambda q: sqldf(q, globals())

	# Write your query in SQL syntax, here you can use df as a normal SQL table
	cond_join= '''
	    SELECT * FROM targets_final
	    JOIN gnomad_final ON ((gnomad_final.[POS_x] >= targets_final.[start]) \
		AND (gnomad_final.[POS_x] <= targets_final.[end]) )
	    '''

	# Now, get your queries results as dataframe using the sqldf object that you created
	polymorph = pysqldf(cond_join)

	# Distance à l'amorce
	polymorph['end'] = polymorph['end'].astype(int)
	polymorph['POS_x'] = polymorph['POS_x'].astype(int)
	# polymorph['Freq_gnomAD'] = polymorph['Freq_gnomAD'].astype(float)

	polymorph['Freq_gnomAD'] = pandas.to_numeric(polymorph['Freq_gnomAD'],errors='coerce')

	polymorph["dist_to_3'"] = polymorph['end'] - polymorph['POS_x']

	# Evaluation
	conditions = [
	    ( (polymorph["dist_to_3'"] <= 10) & (polymorph['Freq_gnomAD'] >= 0.01) ),
	  ( (polymorph["dist_to_3'"] > 10) & (polymorph['Freq_gnomAD'] >= 0.05) ),
	    ]

	choices = ['A revoir', 'A revoir']

	polymorph['Evalutation'] = numpy.select(conditions, choices, default='OK')

	polymorph['Freq_gnomAD'] = polymorph['Freq_gnomAD'].round(decimals=6)

	cols = ['Contig', 'Primer_start', 'Primer_end', 'Primer_name', 'score', \
	    'strand', 'CHROM_x', 'dbSNP_position', 'dbSNP_rs', 'REF', 'ALT', \
	    'Freq_gnomAD', 'dist_to_3', 'Evaluation']

	polymorph.columns = cols    

	noms_genes = df[['gene', 'name_F', 'name_R']]

	final = pandas.merge(polymorph, noms_genes, \
	    left_on='Primer_name', right_on='name_F', how='inner', sort=True, \
	    suffixes=(False, False))

	del final['name_F']
	del final['name_R']

	final['CHROM_x'] = final['CHROM_x'].astype(str)
	final['dbSNP_position'] = final['dbSNP_position'].astype(str)
	final['REF'] = final['REF'].astype(str)
	final['ALT'] = final['ALT'].astype(str)

	hyperlink = pandas.DataFrame()

	hyperlink['gnomAD_hyperlink'] = 'https://gnomad.broadinstitute.org/variant/' \
	    + final['CHROM_x'] + '-' + final['dbSNP_position'] + '-' \
	    + final['REF'] + '-' + final['ALT'] + hl_version

	hl_list = hyperlink['gnomAD_hyperlink'].tolist()
	final['Freq_gnomAD'] = final['Freq_gnomAD'].astype(str)
	freq_list = final['Freq_gnomAD'].tolist()

	final['Conclusion/Derogation'] = ''
	final['Biologiste'] = ''
	final['Date'] = ''

	del final['CHROM_x']

	final = final[['gene', 'Contig', 'Primer_name', 'Primer_start', 'Primer_end', 'score', \
	    'strand', 'dbSNP_rs', 'dbSNP_position', 'REF', 'ALT', \
	    'Freq_gnomAD', 'dist_to_3', \
	    'Evaluation', 'Conclusion/Derogation', \
	    'Biologiste', 'Date']]

	final.drop_duplicates(keep = 'first', inplace=True)

	final.to_csv('./data/polymorph_final.csv', index=False)

	print('ok')

	#########################################
	print('\n--> Export Excel.')

	# Create a Pandas Excel writer using XlsxWriter as the engine.
	writer = pandas.ExcelWriter('./output/'+ filename + '_results_polymorph.xlsx', engine='xlsxwriter')

	final.to_excel(writer, sheet_name='Sheet1', startrow=9, \
	    header=True, index=False)

	# Get the xlsxwriter objects from the dataframe writer object.
	workbook  = writer.book
	# default cell format to size 10 
	workbook.formats[0].set_font_size(9)

	worksheet = writer.sheets['Sheet1']

	# Adapter la largeur des colonnes
	#Iterate through each column and set the width == the max length in that column. 
	#A padding length of 2 is also added.
	for i, col in enumerate(final.columns):
	    # find length of column i
	    column_len = final[col].astype(str).str.len().max()
	    # Setting the length if the column header is larger
	    # than the max column value length
	    column_len = max(column_len, len(col)) + 2
	    # set the column length
	    worksheet.set_column(i, i, column_len)

	# Add a header format.
	header_format = workbook.add_format({
	    'bold': True,
	    'text_wrap': True,
	    'valign': 'top',
	    'fg_color': '#E6E6FA',
	    'border': 1})
	header_format.set_font_size(9)

	# Add a title format
	title_format = workbook.add_format({'bold': True})
	title_format.set_font_size(12)

	# Write the column headers with the defined format.
	for col_num, value in enumerate(final.columns.values):
	    worksheet.write(9, col_num, value, header_format)
	    
	worksheet.write(0, 0, '** POLYMORPH v2.0 **', title_format)

	info_format = title_format = workbook.add_format({'bold': True})
	info_format.set_font_size(10)

	number = datetime.date.today()
	worksheet.write(2, 0, 'Date:', info_format)
	format2 = workbook.add_format({'num_format': 'dd/mm/yy'})
	worksheet.write(2, 2, number, format2) 

	worksheet.write(3, 0, 'Genome version:', info_format)
	worksheet.write(3, 2, 'hg19', info_format)
	
	worksheet.write(4, 0, 'dbSNP version:', info_format)
	worksheet.write(4, 2, '154', info_format)
	
	worksheet.write(5, 0, 'gnomAD genome:', info_format)
	worksheet.write(5, 2, vgAG, info_format)

	worksheet.write(6, 0, 'gnomAD exome:', info_format)
	worksheet.write(6, 2, veAG, info_format)

	worksheet.write(7, 0, 'File:', info_format)
	worksheet.write(7, 2, filename, info_format)

	# Add hyperlinks
	for i in range(len(hl_list)):
		worksheet.write_url(10+i, 11, hl_list[i], string=freq_list[i])

	# Close the Pandas Excel writer and output the Excel file.
	writer.save()
	print('ok')

print('\n**********************')
print('***** JOB DONE ! *****')
print('**********************\n')