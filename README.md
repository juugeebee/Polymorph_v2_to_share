# Polymorph_v2

Logiciel de recherche de polymorphismes connus dans dbSNP dans une liste de primers et annotation de leurs fréquences GnomAD.

Julie Bogoin, juin 2021.

1. Organisation des répertoires racines

./output/
./input/
./data/

2. Fichiers nécessaires au bon fonctionnement du script. Les chemins sont à renseigner dans le fichier db_location.py 

guniq : liste des gènes hg19 au format texte (.txt)
seq : séquence de référence au format fasta (.fa)
gnex : fichier gnomAD exomes 2.1.1 au format .sites.vcf.bgz
gnge : fichier gnomAD genomes 2.1.1 au format .sites.vcf.bgz

3. Utilisation du fichier .yml pour créer un environnement conda contenant les bibliothèques necessaires.

conda env create -f polymorph2_env.yml
conda activate polymorph2_env

4. Commande de lancement

python polymorph_v2_hg19.py 

5. Fichier d'entrée

Autant de fichiers que désiré doivent être placés dans le répertoire input.
Les fichiers doivent être au format texte (.txt) et formaté comme suit (séparateur = tabulation):

#GENE1
Primer1_F_ID	Primer1_F_Sequence	Primer1_R_ID 	Primer1_R_Sequence
Primer2_F_ID	Primer2_F_Sequence	Primer2_R_ID 	Primer2_R_Sequence
#GENE2
Primer3_F_ID	Primer3_F_Sequence	Primer3_R_ID	Primer3_R_Sequence
Primer3_F_ID	Primer3_F_Sequence	Primer4_R_ID 	Primer4_R_Sequence
...

Exemple :
#AR			
AR_seq_F	TCCAGAATCTGTTCCAGAGCGTGC	AR_R	GCTGTGAAGGTTGCTGTTCCTCAT
#TCF4			
PU_TCF4-2F	CCCCAGACTACCAGTACCATCT	PU_TCF4-2R	TTGGGACTACAGGCACATGC
...



Remarques :
    • Un couple de primers F/R est nécessaire à chaque ligne. Il est possible d’avoir une même amorce présente dans deux couples différents. 
    • Les identifiants des primers n’ont pas à être modifiés, la structure « Primer_F_ID » / « Primer_R_ID » n’est utilisée que pour illustration.
    • Eviter les espaces dans les noms de fichiers et dans les fichiers eux-mêmes (noms de gènes, identifiants...).
    • Attention aux noms de gènes indiqués dans le fichier input (#GENE) : cela doit se limiter -strictement- à un nom de gène existant, sinon le logiciel s'arrêtera.
    • Attention également au comportement d’Excel, qui transforme certaines chaînes de caractères en date (par exemple : 1-3 devient 01 mars).

6. Bases de données

Polymorph v2 utilise les dernières versions téléchargeables des bases de données dbSNP et GnomAD, téléchargées en local.

7. Version du génome

Polymorph v2 gère la versions hg19 génome.

8. Fonctionnement

- alignement des paires de primers sur le génome de référence humain par le logiciel in silico PCR.
- croisement des coordonnées génomiques des primers avec le fichier dbSNP afin de trouver les polymorphismes existants.
- calcul de la distance du polymorphisme par rapport à l'extrémité 3' du primer.
- annotation des fréquences GnomAD.

9. Résultats

Un fichier au format <nom du fichier de depart>_results_polymoprh.xlsx est généré dans le dossier output.
Il contient les champs suivants:

- Gene
- Contig
- Primer_name
- Primer_start
- Primer_end
- Score
- Strand
- dbSNP_rs
- dbSNP_position
- REF
- ALT
- Freq_gnomAD : Fréquence gnomAD + lien vers la page de gnomAD du variant.
- dist_to_3'
- Evaluation : Si FreqMax > 5% : valeur « A revoir ». Si FreqMax > 1% et dist_to_3'_end ≤ 10 : valeur « A revoir ». Sinon: valeur « OK ».
- Conclusion/Derogation : champ vide
- Biologiste : champ vide
- Date : champ vide

10. dist_to_3'_end

La dernière base en 3' du primer est considérée comme ayant une distance de 1 par rapport à l'extrémité 3'.
Dans le cas des indels, nous indiquons la distance la plus proximale de l'extrémité 3' ainsi que la distance la plus distale de l'extrémité 3'. En cas d’indel partiellement comprise dans la séquence du primer, seules les bases communes sont mentionnées.
