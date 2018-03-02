
The idea is to mine the literature to find associations between keywords such as:

- Quadruplex terms:
  - Direct terms: quadruplex, G-quadruplex, G-rich sequence, G4, guanine quartet, tetraplex, tetramer
  - Related terms: telomer, guanine, DNA binding, RNA binding, helicase, replication origin

- Lists of genes that Darcie has shared


## List of genes

In `uk-cri-lcst01`,

```bash
cd /lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data
mkdir -p 20170727/tables
cd 20170727/tables
rsync -arvuP martin03@sblab-srv001:~/Downloads/List_for_literature_search.xls .
```


## Ideas

I posted a question in Biostars to try to draw from people's expertise:

https://www.biostars.org/p/266920/


bio.tools:
- Text mining https://bio.tools/?page=1&q=Text%20mining&sort=score&ord=desc
- Literature search https://bio.tools/?page=1&q=Literature%20search&sort=score&ord=desc

Methods:
- PolySearch2 http://polysearch.cs.ualberta.ca/index

- Tool developed for STRING [database](https://string-db.org/)

- BEST http://best.korea.ac.kr/ , how to download results?

- Europe PMC
  - SciLite https://europepmc.org/Annotations

- NCBI https://www.ncbi.nlm.nih.gov/research/bionlp/tools/
  - PubTator https://www.ncbi.nlm.nih.gov/CBBresearch/Lu/Demo/PubTator/

- DeepLife https://gate.d5.mpi-inf.mpg.de/deeplife/en-health/

- Chilibot http://www.chilibot.net/

- GO G-quadruplex:
  - GO:0051880 G-quadruplex DNA binding https://www.ebi.ac.uk/QuickGO/term/GO:0051880 (function)
  - GO:0002151 G-quadruplex RNA binding https://www.ebi.ac.uk/QuickGO/term/GO:0002151 (function)
  - GO:0061849 telomeric G-quadruplex DNA binding https://www.ebi.ac.uk/QuickGO/term/GO:0061849 (function)
  - GO:0071919 G-quadruplex DNA formation https://www.ebi.ac.uk/QuickGO/term/GO:0071919 (biological process)
  - GO:0044806 G-quadruplex DNA unwinding https://www.ebi.ac.uk/QuickGO/term/GO:0044806 (biological process)
  - GO:1905493 regulation of G-quadruplex DNA binding https://www.ebi.ac.uk/QuickGO/term/GO:1905493 (biological process)
  - GO:1905494 negative regulation of G-quadruplex DNA binding https://www.ebi.ac.uk/QuickGO/term/GO:1905494 (biological process)
  - GO:1905495 positive regulation of G-quadruplex DNA binding https://www.ebi.ac.uk/QuickGO/term/GO:1905495 (biological process)
  - GO:1905465 regulation of G-quadruplex DNA unwinding https://www.ebi.ac.uk/QuickGO/term/GO:1905465 (biological process)
  - GO:1905466 negative regulation of G-quadruplex DNA unwinding https://www.ebi.ac.uk/QuickGO/term/GO:1905466 (biological process)
  - GO:1905467 positive regulation of G-quadruplex DNA unwinding https://www.ebi.ac.uk/QuickGO/term/GO:1905467 (biological process)

- GO G-rich:
  - GO:1990955 G-rich single-stranded DNA binding https://www.ebi.ac.uk/QuickGO/term/GO:1990955 (function)
  - GO:0098505 G-rich strand telomeric DNA binding https://www.ebi.ac.uk/QuickGO/term/GO:0098505 (function)


- Uniprot: scan human genes for quadruplex associated terms and rank


Help:
- Biostars https://www.biostars.org/local/search/page/?q=text+mining

People and organisations:
- Senay Kafkas
- Chen Li
- Dietrich Schumann
- Jo McEntire
- Martin Krallinger / Alfonso Valencia
- OpenTargets

Journals:
- Bioinformatics
- Briefings in bioinformatics

Previous work:
- G4IPDB http://bsbe.iiti.ac.in/bsbe/ipdb/index.php



## GO and Uniprot

1. *High confidence 1*.
List of Homo sapiens proteins from Uniprot annotated with at least one of the GO G-quadruplex terms above.
GO Evidence: Any assertion method.
Saved as `20170804_uniprot_go_G-quadruplex_homosapiens.tab`.
[GO G-quadruplex url](http://www.uniprot.org/uniprot/?query=goa:(go:%22G-quadruplex%20DNA%20binding%20[0051880]%22)%20OR%20goa:(go:%22G-quadruplex%20RNA%20binding%20[0002151]%22)%20OR%20goa:(go:%22telomeric%20G-quadruplex%20DNA%20binding%20[0061849]%22)%20OR%20goa:(go:%22G-quadruplex%20DNA%20formation%20[0071919]%22)%20OR%20goa:(go:%22G-quadruplex%20DNA%20unwinding%20[0044806]%22)%20OR%20goa:(go:%22regulation%20of%20G-quadruplex%20DNA%20binding%20[1905493]%22)%20OR%20goa:(go:%22negative%20regulation%20of%20G-quadruplex%20DNA%20binding%20[1905494]%22)%20OR%20goa:(go:%22positive%20regulation%20of%20G-quadruplex%20DNA%20binding%20[1905495]%22)%20OR%20goa:(go:%22regulation%20of%20G-quadruplex%20DNA%20unwinding%20[1905465]%22)%20OR%20goa:(go:%22negative%20regulation%20of%20G-quadruplex%20DNA%20unwinding%20[1905466]%22)%20OR%20goa:(go:%22positive%20regulation%20of%20G-quadruplex%20DNA%20unwinding%20[1905467]%22)&fil=organism%3A%22Homo+sapiens+%28Human%29+%5B9606%5D%22&sort=score).

2. *High confidence 2*.
List of Homo sapiens proteins from Uniprot annotated with at least one of the GO G-quadruplex terms above or where the term `quadruplex` appears somewhere in the uniprot entry.
GO Evidence: Any assertion method.
Saved as `20170804_uniprot_go_G-quadruplex_or_quadruplex_term_homosapiens.tab`.
[GO G-quadruplex or quadruplex term url](http://www.uniprot.org/uniprot/?query=goa:(go:%22G-quadruplex%20DNA%20binding%20[0051880]%22)%20OR%20goa:(go:%22G-quadruplex%20RNA%20binding%20[0002151]%22)%20OR%20goa:(go:%22telomeric%20G-quadruplex%20DNA%20binding%20[0061849]%22)%20OR%20goa:(go:%22G-quadruplex%20DNA%20formation%20[0071919]%22)%20OR%20goa:(go:%22G-quadruplex%20DNA%20unwinding%20[0044806]%22)%20OR%20goa:(go:%22regulation%20of%20G-quadruplex%20DNA%20binding%20[1905493]%22)%20OR%20goa:(go:%22negative%20regulation%20of%20G-quadruplex%20DNA%20binding%20[1905494]%22)%20OR%20goa:(go:%22positive%20regulation%20of%20G-quadruplex%20DNA%20binding%20[1905495]%22)%20OR%20goa:(go:%22regulation%20of%20G-quadruplex%20DNA%20unwinding%20[1905465]%22)%20OR%20goa:(go:%22negative%20regulation%20of%20G-quadruplex%20DNA%20unwinding%20[1905466]%22)%20OR%20goa:(go:%22positive%20regulation%20of%20G-quadruplex%20DNA%20unwinding%20[1905467]%22)%20OR%20quadruplex&fil=reviewed%3Ayes+AND+organism%3A%22Homo+sapiens+%28Human%29+%5B9606%5D%22&sort=score).

3. *Medium confidence 1*.
List of Homo sapiens proteins from Uniprot annotated with at least one of the GO G-quadruplex or GO G-rich terms above.
GO Evidence: Any assertion method.
Saved as `20170804_uniprot_go_G-quadruplex_or_go_G-rich_homosapiens.tab`
[GO G-quadruplex or GO G-rich url](http://www.uniprot.org/uniprot/?query=goa:(go:%22G-quadruplex%20DNA%20binding%20[0051880]%22)%20OR%20goa:(go:%22G-quadruplex%20RNA%20binding%20[0002151]%22)%20OR%20goa:(go:%22telomeric%20G-quadruplex%20DNA%20binding%20[0061849]%22)%20OR%20goa:(go:%22G-quadruplex%20DNA%20formation%20[0071919]%22)%20OR%20goa:(go:%22G-quadruplex%20DNA%20unwinding%20[0044806]%22)%20OR%20goa:(go:%22regulation%20of%20G-quadruplex%20DNA%20binding%20[1905493]%22)%20OR%20goa:(go:%22negative%20regulation%20of%20G-quadruplex%20DNA%20binding%20[1905494]%22)%20OR%20goa:(go:%22positive%20regulation%20of%20G-quadruplex%20DNA%20binding%20[1905495]%22)%20OR%20goa:(go:%22regulation%20of%20G-quadruplex%20DNA%20unwinding%20[1905465]%22)%20OR%20goa:(go:%22negative%20regulation%20of%20G-quadruplex%20DNA%20unwinding%20[1905466]%22)%20OR%20goa:(go:%22positive%20regulation%20of%20G-quadruplex%20DNA%20unwinding%20[1905467]%22)%20OR%20goa:(go:%22G-rich%20single-stranded%20DNA%20binding%20[1990955]%22)%20OR%20goa:(go:%22G-rich%20strand%20telomeric%20DNA%20binding%20[0098505]%22)&fil=organism%3A%22Homo+sapiens+%28Human%29+%5B9606%5D%22&sort=score).

4. *Medium confidence 2*.
List of Homo sapiens proteins from Uniprot annotated with at least one of the GO G-quadruplex or GO G-rich terms above or where the terms `quadruplex` or `G-rich` appear somewhere in the uniprot entry.
GO Evidence: Any assertion method.
Saved as `20170804_uniprot_go_G-quadruplex_or_go_G-rich_or_quadruplex_term_or_G-rich_term_homosapiens.tab`
[GO G-quadruplex or GO G-rich or quadruplex term or G-rich term url](http://www.uniprot.org/uniprot/?query=goa:(go:%22G-quadruplex%20DNA%20binding%20[0051880]%22)%20OR%20goa:(go:%22G-quadruplex%20RNA%20binding%20[0002151]%22)%20OR%20goa:(go:%22telomeric%20G-quadruplex%20DNA%20binding%20[0061849]%22)%20OR%20goa:(go:%22G-quadruplex%20DNA%20formation%20[0071919]%22)%20OR%20goa:(go:%22G-quadruplex%20DNA%20unwinding%20[0044806]%22)%20OR%20goa:(go:%22regulation%20of%20G-quadruplex%20DNA%20binding%20[1905493]%22)%20OR%20goa:(go:%22negative%20regulation%20of%20G-quadruplex%20DNA%20binding%20[1905494]%22)%20OR%20goa:(go:%22positive%20regulation%20of%20G-quadruplex%20DNA%20binding%20[1905495]%22)%20OR%20goa:(go:%22regulation%20of%20G-quadruplex%20DNA%20unwinding%20[1905465]%22)%20OR%20goa:(go:%22negative%20regulation%20of%20G-quadruplex%20DNA%20unwinding%20[1905466]%22)%20OR%20goa:(go:%22positive%20regulation%20of%20G-quadruplex%20DNA%20unwinding%20[1905467]%22)%20OR%20goa:(go:%22G-rich%20single-stranded%20DNA%20binding%20[1990955]%22)%20OR%20goa:(go:%22G-rich%20strand%20telomeric%20DNA%20binding%20[0098505]%22)%20OR%20goa:(quadruplex)%20OR%20%22g%20rich%22&fil=reviewed%3Ayes+AND+organism%3A%22Homo+sapiens+%28Human%29+%5B9606%5D%22&sort=score).


### Copy tables

```bash
cd /lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20170727
mkdir uniprot_go
cd uniprot_go
rsync -arvuP martin03@nm149s012925:/Users/martin03/Desktop/*.tab .
```


## PolySearch2

This will generate a *low confidence* list of associations between quadruplex and gene names.

PolySearch2 is cited in Krallinger2017.

MeSH term id for G-quadruplex: [D054856](https://www.ncbi.nlm.nih.gov/mesh/68054856).

Species: Coleophora quadruplex ???

- Search ID: **1501862729** (35)
  - Given: Text
  - Query Keyword: Quadruplex
  - Find ALL associated: Genes/Proteins
  - Automated synonym list: quadruplex;G-Quadruplex;G-quadruplex;Quadruplexes;quadruplexes;G-Quadruplexes;G-quadruplexes;Tetrad;tetrad;Tetrads;tetrads;Quartet;quartet;Quartets;quartets;Tetraplex;tetraplex;Tetraplexes;tetraplexes
  - Custom Negation Words, separate words with ";": (blank)
  - Select one or more Corpus to search: PubMed, PubMed Central
  - Document Limit: No Limit
  - Show borderline associations: (blank)

- Search ID: **1501863075** (633)
  - Given: Text
  - Query Keyword: Quadruplex
  - Find ALL associated: Genes/Proteins
  - Automated synonym list: quadruplex;G-Quadruplex;G-quadruplex;Quadruplexes;quadruplexes;G-Quadruplexes;G-quadruplexes;Tetrad;tetrad;Tetrads;tetrads;Quartet;quartet;Quartets;quartets;Tetraplex;tetraplex;Tetraplexes;tetraplexes
  - Custom Negation Words, separate words with ";": (blank)
  - Select one or more Corpus to search: PubMed, PubMed Central
  - Document Limit: No Limit
  - Show borderline associations: (ticked)

- Search ID: **1501864605** (35)
  - Given: Text
  - Query Keyword: Quadruplex
  - Find ALL associated: Genes/Proteins
  - Automated synonym list: quadruplex;G-Quadruplex;G-quadruplex;Quadruplexes;quadruplexes;G-Quadruplexes;G-quadruplexes;Tetrad;tetrad;Tetrads;tetrads;Quartet;quartet;Quartets;quartets;Tetraplex;tetraplex;Tetraplexes;tetraplexes
  - Custom Negation Words, separate words with ";": (blank)
  - Select one or more Corpus to search: PubMed, PubMed Central, Wikipedia, USPTO Patent, NCBI Books, MedlinePlus
  - Document Limit: No Limit
  - Show borderline associations: (blank)

- Search ID: **1501864803** (???)
  - Given: Text
  - Query Keyword: Quadruplex
  - Find ALL associated: Genes/Proteins
  - Automated synonym list: quadruplex;G-Quadruplex;G-quadruplex;Quadruplexes;quadruplexes;G-Quadruplexes;G-quadruplexes;Tetrad;tetrad;Tetrads;tetrads;Quartet;quartet;Quartets;quartets;Tetraplex;tetraplex;Tetraplexes;tetraplexes
  - Custom Negation Words, separate words with ";": (blank)
  - Select one or more Corpus to search: PubMed, PubMed Central, Wikipedia, USPTO Patent, NCBI Books, MedlinePlus
  - Document Limit: No Limit
  - Show borderline associations: (ticked)

- Search ID: **1502210485** (559)
  - Given: Text
  - Query Keyword: "G-quadruplex"
  - Find ALL associated: Genes/Proteins
  - Automated synonym list: (blank)
  - Custom Negation Words, separate words with ";": (blank)
  - Select one or more Corpus to search: PubMed, PubMed Central
  - Document Limit: No Limit
  - Show borderline associations: (ticked)

- Search ID: **1502210609** (???)
  - Given: Text
  - Query Keyword: "G-quadruplex"
  - Find ALL associated: Genes/Proteins
  - Automated synonym list: quadruplex;G-Quadruplex;G-quadruplex;Quadruplexes;quadruplexes;G-Quadruplexes;G-quadruplexes;Tetrad;tetrad;Tetrads;tetrads;Quartet;quartet;Quartets;quartets;Tetraplex;tetraplex;Tetraplexes;tetraplexes
  - Custom Negation Words, separate words with ";": (blank)
  - Select one or more Corpus to search: PubMed, PubMed Central
  - Document Limit: No Limit
  - Show borderline associations: (ticked)




### Copy tables

All the files are in `/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20170727/polysearch2`. I still need to copy the file `1501864803` and `1502210609` there. In `uk-cri-lcst01`,

```bash
cd /lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20170727/polysearch2
ls -lh
#-rw-r--r-- 1 martin03 sblab 4.4M Aug  4 10:07 POLYSEARCH-1501862729.json
#-rw-r--r-- 1 martin03 sblab 8.6M Aug  4 10:14 POLYSEARCH-1501863075.json
#-rw-r--r-- 1 martin03 sblab 4.5M Aug  4 10:38 POLYSEARCH-1501864605.json
#-rw-r--r-- 1 martin03 sblab 7.7M Aug  8 10:44 POLYSEARCH-1502210485.json
#-rw-r--r-- 1 martin03 sblab  25K Aug  7  2014 POLYSEARCH2_THESAURUS_gene_family.txt
#-rw-r--r-- 1 martin03 sblab 6.6M Aug  7  2014 POLYSEARCH2_THESAURUS_genes.txt
```

The files `POLYSEARCH2_THESAURUS_genes.txt` and `POLYSEARCH2_THESAURUS_gene_family.txt` links ids of genes (PS...) and gene families (PSF...) given in the result `.json` files


### Feedback from Wishart research group

My query:

`Dear Wishart research group,
I am using PolySearch2 to try to discover associations between the term: G-quadruplex and gene names in Pubmed abstracts. I have read the Documentation in detail and I have already performed some searches for testing (e.g. Search ID: 1501858201). I'd like to know how I can use the MeSH id for G-quadruplexes: D054856 in the "Given" for a "Find ALL associated": Genes/proteins. E.g.I have already tried Given:Text and Query Keyword:D054856 however I don't get any result and I am not sure how to search associations between a MeSH term and all gene names. It would be great if you could send some tips how to perform this search effectively.
Thanks in advance,
Sergio`

Response:

`Hi Sergio,
Try using "G-quadruplex", with double quotes, as your query keyword and "Genes/Protein" as your query type and target type.
For example, quicksearch results by doing the above: http://polysearch.cs.ualberta.ca/assoresult/1502206193
Regards,
Elvis`




## Analysis

### GO and Uniprot

I am going to work with the *High confidence 2* dataset.

```python

# Loading screen lists of genes
ifile1 = open("/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20170727/tables/List_for_literature_search.csv", "r")
ilines1 = ifile1.read().split("\r")
ifile1.close()

list_withoutFC = []
list_withFC = []

for line in ilines1[1:]:
  fields = line.split(",")
  if fields[0] != "":
    list_withoutFC.append(fields[0])
  if fields[2] != "":
    list_withFC.append(fields[2])


len(list_withoutFC) # 843 without FC threshold genes
len(set(list_withoutFC)) # 843 unique without FC threshold genes
len(list_withFC) # 763 with FC threshold genes
len(set(list_withFC)) # 763 unique with FC threshold genes


# Loading GO and Uniprot High confidence 2 dataset
ifile2 = open("/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20170727/uniprot_go/20170804_uniprot_go_G-quadruplex_or_quadruplex_term_homosapiens.tab", "r")
ilines2 = ifile2.readlines()
ifile2.close()

list_go_uniprot_hc2 = []

for line in ilines2[1:]:
  fields = line.split("\t")
  gene_names = fields[2].split()
  list_go_uniprot_hc2 += gene_names


len(list_go_uniprot_hc2) # 46 gene names GO and Uniprot High confidence 2
len(set(list_go_uniprot_hc2)) # 46 unique gene names GO and Uniprot High confidence 2
len(ilines2[1:]) # 18 proteins GO and Uniprot High confidence 2


# Intersection withoutFC
len(set(list_withoutFC).intersection(set(list_go_uniprot_hc2))) # 5 out of 18 are found in common
set(list_withoutFC).intersection(set(list_go_uniprot_hc2)) # set(['ATRX', 'DNA2', 'XRN1', 'DHX36', 'MCRS1'])

# Differences withoutFC
len(set(list_withoutFC).difference(set(list_go_uniprot_hc2))) # 838 out of 843 are without FC threshold genes only
len(set(list_go_uniprot_hc2).difference(set(list_withoutFC))) # 41 (13) out of 46 (18) are without FC threshold genes only

# Intersection withFC
len(set(list_withFC).intersection(set(list_go_uniprot_hc2))) # 5 out of 18 are found in common
set(list_withFC).intersection(set(list_go_uniprot_hc2)) # set(['ATRX', 'DNA2', 'XRN1', 'DHX36', 'MCRS1'])

# Differences withFC
len(set(list_withFC).difference(set(list_go_uniprot_hc2))) # 758 out of 763 are without FC threshold genes only
len(set(list_go_uniprot_hc2).difference(set(list_withFC))) # 41 (13) out of 46 (18) are without FC threshold genes only


# Output file
ofile = open("/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20170727/uniprot_go/20170808_uniprot_go_G-quadruplex_or_quadruplex_term_homosapiens_screen.txt", "w")

ofile.write("Entry\tEntry name\tGene names\twithout FC threshold?\twith FC threshold?\tFunction [CC]\tGene ontology IDs\tGene ontology (molecular function)\tGene ontology (biological process)\tGene ontology (cellular component)\tInvolvement in disease\n")

for line in ilines2[1:]:
  fields = line.split("\t")
  ofile.write("\t".join(fields[:3]) + "\t")
  gene_names = fields[2].split()
  if len(set(gene_names).intersection(set(list_withoutFC))) > 0:
    ofile.write("yes\t")
  else:
    ofile.write("no\t")
  if len(set(gene_names).intersection(set(list_withFC))) > 0:
    ofile.write("yes\t")
  else:
    ofile.write("no\t")
  ofile.write("\t".join(fields[3:]))


ofile.close()

```



### PolySearch2

I am going to use Search ID: **1501863075**. Exploring POLYSEARCH-1501863075.json:

- doc_count:
  - medline: 4233
  - pmc: 1244
- found_hits: 633
- num_hits: 5477

```python
# Loading screen lists of genes
ifile1 = open("/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20170727/tables/List_for_literature_search.csv", "r")
ilines1 = ifile1.read().split("\r")
ifile1.close()

list_withoutFC = []
list_withFC = []

for line in ilines1[1:]:
  fields = line.split(",")
  if fields[0] != "":
    list_withoutFC.append(fields[0])
  if fields[2] != "":
    list_withFC.append(fields[2])


len(list_withoutFC) # 843 without FC threshold genes
len(set(list_withoutFC)) # 843 unique without FC threshold genes
len(list_withFC) # 763 with FC threshold genes
len(set(list_withFC)) # 763 unique with FC threshold genes


# Loading thesauri into dictionaries
## POLYSEARCH2_THESAURUS_genes.txt
t_genes_ifile = open("/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20170727/polysearch2/POLYSEARCH2_THESAURUS_genes.txt", "r")
t_genes_ilines = t_genes_ifile.readlines()
t_genes_ifile.close()

t_genes = {}

for line in t_genes_ilines:
  fields = line.split("\t")
  entity = fields[0]
  t_genes[entity] = []
  for name in fields[1:]:
    name = name.replace("\n", "")
    t_genes[entity].append(name)


len(t_genes) # 27994
t_genes["PS12860"] # ['DNA topoisomerase 1', 'DNA topoisomerase I', 'TOP-1', 'TOP1', 'TOP1 protein', 'topoisomerase (DNA) I', 'DNA topoisomerase Is', 'TOP1 proteins', 'topoisomerase (DNA) Is', 'Type I DNA topoisomerase', 'EC 5.99.1.2']

## POLYSEARCH2_THESAURUS_gene_family.txt
t_genefamily_ifile = open("/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20170727/polysearch2/POLYSEARCH2_THESAURUS_gene_family.txt", "r")
t_genefamily_ilines = t_genefamily_ifile.readlines()
t_genefamily_ifile.close()

t_genefamily = {}

for line in t_genefamily_ilines:
  fields = line.split("\t")
  entity = fields[0]
  t_genefamily[entity] = []
  for name in fields[1:]:
    name = name.replace("\n", "")
    t_genefamily[entity].append(name)


len(t_genefamily) # 404
t_genefamily["PSF00377"] # ['Telomerase', 'Telomerases']


# Loading POLYSEARCH-1501863075.json and output file
import json
ifile2 = json.load(open("/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20170727/polysearch2/POLYSEARCH-1501863075.json"))

ofile = open("/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20170727/polysearch2/20170808_polysearch2_1501863075_quadruplex_pubmed_pmc_screen.txt", "w")

ofile.write("Gene names\tType\tScore\twithout FC threshold?\twith FC threshold?\n")

for entity in ifile2["hits"]["gene"]:
  i = entity["eid"]
  rscore = entity["rscore"]
  synonyms = entity["esynonyms"]
  if i.startswith("PSF"):
    synonyms2 = t_genefamily[i]
    ty = "Gene family"
  else:
    synonyms2 = t_genes[i]
    ty = "Gene"
  if len(set(synonyms2).intersection(set(list_withoutFC))) > 0:
    withoutFC = "yes"
  else:
    withoutFC = "no"
  if len(set(synonyms2).intersection(set(list_withFC))) > 0:
    withFC = "yes"
  else:
    withFC = "no"
  ofile.write("%s\n" % "\t".join([";".join(synonyms2), ty, str(rscore), withoutFC, withFC]))
#  print i, synonyms2, ty, rscore, withoutFC, withFC


ofile.close()

```

Further analysis:

```bash
cd /lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20170727/polysearch2
wc -l 20170808_polysearch2_1501863075_quadruplex_pubmed_pmc_screen.txt # 634
```

```r
library(data.table)

data <- fread("/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20170727/polysearch2/20170808_polysearch2_1501863075_quadruplex_pubmed_pmc_screen.txt")

table(data$Type)
#       Gene Gene family
#        582          51

table(data.frame(data)[,4]) # without FC threshold?
#no yes
#575  58

table(data.frame(data)[,5]) # with FC threshold?
#no yes
#581  52

```





## Analysis (including columns requested by Darcie)

### GO and Uniprot

I am going to work with the *High confidence 2* dataset here too.

In `uk-cri-lcst01`, copying tables across:

```bash
cd /lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20170727/tables
rsync -arvuP martin03@nm149s012925:/media/itshare/research/sblab/group/public_folders/martin03/20160418_shRNAscreen_katie_darcie/20170130_integration_all/tables/genes_all/*.txt .
```

Analysis,

```python

# Loading screen lists of genes
ifile1 = open("/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20170727/tables/List_for_literature_search.csv", "r")
ilines1 = ifile1.read().split("\r")
ifile1.close()

list_withoutFC = []
list_withFC = []

for line in ilines1[1:]:
  fields = line.split(",")
  if fields[0] != "":
    list_withoutFC.append(fields[0])
  if fields[2] != "":
    list_withFC.append(fields[2])


len(list_withoutFC) # 843 without FC threshold genes
len(set(list_withoutFC)) # 843 unique without FC threshold genes
len(list_withFC) # 763 with FC threshold genes
len(set(list_withFC)) # 763 unique with FC threshold genes


# Loading GO and Uniprot High confidence 2 dataset
ifile2 = open("/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20170727/uniprot_go/20170804_uniprot_go_G-quadruplex_or_quadruplex_term_homosapiens.tab", "r")
ilines2 = ifile2.readlines()
ifile2.close()

list_go_uniprot_hc2 = []

for line in ilines2[1:]:
  fields = line.split("\t")
  gene_names = fields[2].split()
  list_go_uniprot_hc2 += gene_names


len(list_go_uniprot_hc2) # 46 gene names GO and Uniprot High confidence 2
len(set(list_go_uniprot_hc2)) # 46 unique gene names GO and Uniprot High confidence 2
len(ilines2[1:]) # 18 proteins GO and Uniprot High confidence 2


# Intersection withoutFC
len(set(list_withoutFC).intersection(set(list_go_uniprot_hc2))) # 5 out of 18 are found in common
set(list_withoutFC).intersection(set(list_go_uniprot_hc2)) # set(['ATRX', 'DNA2', 'XRN1', 'DHX36', 'MCRS1'])

# Differences withoutFC
len(set(list_withoutFC).difference(set(list_go_uniprot_hc2))) # 838 out of 843 are without FC threshold genes only
len(set(list_go_uniprot_hc2).difference(set(list_withoutFC))) # 41 (13) out of 46 (18) are without FC threshold genes only

# Intersection withFC
len(set(list_withFC).intersection(set(list_go_uniprot_hc2))) # 5 out of 18 are found in common
set(list_withFC).intersection(set(list_go_uniprot_hc2)) # set(['ATRX', 'DNA2', 'XRN1', 'DHX36', 'MCRS1'])

# Differences withFC
len(set(list_withFC).difference(set(list_go_uniprot_hc2))) # 758 out of 763 are without FC threshold genes only
len(set(list_go_uniprot_hc2).difference(set(list_withFC))) # 41 (13) out of 46 (18) are without FC threshold genes only


# Loading PDS and PhenDC3 tables
ifile_pds = open("/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20170727/tables/20170130_PDSvst0_genes.txt", "r")
ilines_pds = ifile_pds.readlines()
ifile_pds.close()

pds_dict = {}

for line in ilines_pds[1:]:
  fields = line.split()
  gene = fields[0]
  n_detected_shRNAs = fields[2]
  n_significant_sensitiser_shRNAs = fields[11]
  logFC_significant_sensitiser_shRNAs_median = fields[14]
  pds_dict[gene] = [n_detected_shRNAs, n_significant_sensitiser_shRNAs, logFC_significant_sensitiser_shRNAs_median]


pds_dict["TOP1"] # ['8', '3', '-3.12857659580522']
pds_dict["MYC"] # ['8', '6', '-2.99796523521497']
pds_dict["ATRX"] # ['16', '1', '-1.55023350217559']

ifile_phendc3 = open("/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20170727/tables/20170130_PhenDC3vst0_genes.txt", "r")
ilines_phendc3 = ifile_phendc3.readlines()
ifile_phendc3.close()

phendc3_dict = {}

for line in ilines_phendc3[1:]:
  fields = line.split()
  gene = fields[0]
  n_detected_shRNAs = fields[2]
  n_significant_sensitiser_shRNAs = fields[11]
  logFC_significant_sensitiser_shRNAs_median = fields[14]
  phendc3_dict[gene] = [n_detected_shRNAs, n_significant_sensitiser_shRNAs, logFC_significant_sensitiser_shRNAs_median]


phendc3_dict["TOP1"] # ['8', '4', '-3.92126848474496']
phendc3_dict["MYC"] # ['8', '0', 'NA']
phendc3_dict["ATRX"] # ['16', '5', '-1.39558644062454']


# Output file
ofile = open("/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20170727/uniprot_go/20170810_uniprot_go_G-quadruplex_or_quadruplex_term_homosapiens_screen_summary.txt", "w")

ofile.write("Gene name\twithout FC threshold?\twith FC threshold?\tn_detected_shRNAs\tn_significant_sensitiser_shRNAs_PDS\tlogFC_significant_sensitiser_shRNAs_median_PDS\tn_significant_sensitiser_shRNAs_PhenDC3\tlogFC_significant_sensitiser_shRNAs_median_PhenDC3\tQuadruplex GO term / UniprotKB association\n")

for line in ilines2[1:]:
  fields = line.split("\t")
  gene_names = fields[2].split()
  gene_name = gene_names[0]
  ofile.write("%s\t" % gene_name)
  if len(set(gene_names).intersection(set(list_withoutFC))) > 0:
    ofile.write("yes\t")
  else:
    ofile.write("no\t")
  if len(set(gene_names).intersection(set(list_withFC))) > 0:
    ofile.write("yes\t")
  else:
    ofile.write("no\t")
  ofile.write("\t".join(pds_dict[gene_name]) + "\t")
  ofile.write("\t".join(phendc3_dict[gene_name][1:]) + "\t")
  go = [g for g in fields[5].split("; ") if "quadruplex" in g] + [g for g in fields[6].split("; ") if "quadruplex" in g]
  if go != []:
    ofile.write("; ".join(go) + '\n')
  else:
    ofile.write("quadruplex term present within UniprotKB:%s entry" % fields[0] + '\n')


ofile.close()

```



### PolySearch2

I am going to use Search ID: **1501863075**. Exploring POLYSEARCH-1501863075.json:

- doc_count:
  - medline: 4233
  - pmc: 1244
- found_hits: 633
- num_hits: 5477

```python
# Loading screen lists of genes
ifile1 = open("/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20170727/tables/List_for_literature_search.csv", "r")
ilines1 = ifile1.read().split("\r")
ifile1.close()

list_withoutFC = []
list_withFC = []

for line in ilines1[1:]:
  fields = line.split(",")
  if fields[0] != "":
    list_withoutFC.append(fields[0])
  if fields[2] != "":
    list_withFC.append(fields[2])


len(list_withoutFC) # 843 without FC threshold genes
len(set(list_withoutFC)) # 843 unique without FC threshold genes
len(list_withFC) # 763 with FC threshold genes
len(set(list_withFC)) # 763 unique with FC threshold genes


# Loading thesauri into dictionaries
## POLYSEARCH2_THESAURUS_genes.txt
t_genes_ifile = open("/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20170727/polysearch2/POLYSEARCH2_THESAURUS_genes.txt", "r")
t_genes_ilines = t_genes_ifile.readlines()
t_genes_ifile.close()

t_genes = {}

for line in t_genes_ilines:
  fields = line.split("\t")
  entity = fields[0]
  t_genes[entity] = []
  for name in fields[1:]:
    name = name.replace("\n", "")
    t_genes[entity].append(name)


len(t_genes) # 27994
t_genes["PS12860"] # ['DNA topoisomerase 1', 'DNA topoisomerase I', 'TOP-1', 'TOP1', 'TOP1 protein', 'topoisomerase (DNA) I', 'DNA topoisomerase Is', 'TOP1 proteins', 'topoisomerase (DNA) Is', 'Type I DNA topoisomerase', 'EC 5.99.1.2']

## POLYSEARCH2_THESAURUS_gene_family.txt
t_genefamily_ifile = open("/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20170727/polysearch2/POLYSEARCH2_THESAURUS_gene_family.txt", "r")
t_genefamily_ilines = t_genefamily_ifile.readlines()
t_genefamily_ifile.close()

t_genefamily = {}

for line in t_genefamily_ilines:
  fields = line.split("\t")
  entity = fields[0]
  t_genefamily[entity] = []
  for name in fields[1:]:
    name = name.replace("\n", "")
    t_genefamily[entity].append(name)


len(t_genefamily) # 404
t_genefamily["PSF00377"] # ['Telomerase', 'Telomerases']


# Loading PDS and PhenDC3 tables
ifile_pds = open("/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20170727/tables/20170130_PDSvst0_genes.txt", "r")
ilines_pds = ifile_pds.readlines()
ifile_pds.close()

pds_dict = {}

for line in ilines_pds[1:]:
  fields = line.split()
  gene = fields[0]
  n_detected_shRNAs = fields[2]
  n_significant_sensitiser_shRNAs = fields[11]
  logFC_significant_sensitiser_shRNAs_median = fields[14]
  pds_dict[gene] = [n_detected_shRNAs, n_significant_sensitiser_shRNAs, logFC_significant_sensitiser_shRNAs_median]


pds_dict["TOP1"] # ['8', '3', '-3.12857659580522']
pds_dict["MYC"] # ['8', '6', '-2.99796523521497']
pds_dict["ATRX"] # ['16', '1', '-1.55023350217559']

ifile_phendc3 = open("/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20170727/tables/20170130_PhenDC3vst0_genes.txt", "r")
ilines_phendc3 = ifile_phendc3.readlines()
ifile_phendc3.close()

phendc3_dict = {}

for line in ilines_phendc3[1:]:
  fields = line.split()
  gene = fields[0]
  n_detected_shRNAs = fields[2]
  n_significant_sensitiser_shRNAs = fields[11]
  logFC_significant_sensitiser_shRNAs_median = fields[14]
  phendc3_dict[gene] = [n_detected_shRNAs, n_significant_sensitiser_shRNAs, logFC_significant_sensitiser_shRNAs_median]


phendc3_dict["TOP1"] # ['8', '4', '-3.92126848474496']
phendc3_dict["MYC"] # ['8', '0', 'NA']
phendc3_dict["ATRX"] # ['16', '5', '-1.39558644062454']


# Loading POLYSEARCH-1501863075.json and writing output
import json
ifile2 = json.load(open("/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20170727/polysearch2/POLYSEARCH-1501863075.json"))

ofile = open("/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20170727/polysearch2/20170810_polysearch2_1501863075_quadruplex_pubmed_pmc_screen_summary.txt", "w")

ofile.write("Gene name\twithout FC threshold?\twith FC threshold?\tn_detected_shRNAs\tn_significant_sensitiser_shRNAs_PDS\tlogFC_significant_sensitiser_shRNAs_median_PDS\tn_significant_sensitiser_shRNAs_PhenDC3\tlogFC_significant_sensitiser_shRNAs_median_PhenDC3\tScore\n")

all_genes = [line.split()[0] for line in ilines_pds[1:]]

for entity in ifile2["hits"]["gene"]:
  i = entity["eid"]
  rscore = entity["rscore"]
  synonyms = entity["esynonyms"]
  if i.startswith("PSF"):
    synonyms2 = t_genefamily[i]
    ty = "Gene family"
  else:
    synonyms2 = t_genes[i]
    ty = "Gene"
  if len(set(synonyms2).intersection(set(list_withoutFC))) > 0:
    withoutFC = "yes"
  else:
    withoutFC = "no"
  if len(set(synonyms2).intersection(set(list_withFC))) > 0:
    withFC = "yes"
  else:
    withFC = "no"
  if ty == "Gene":
    inter = set(synonyms2).intersection(set(all_genes))
    if len(inter) > 0:
      if len(inter) == 1:
        gene_name = list(inter)[0]
        n_detected_shRNAs = pds_dict[gene_name][0]
        n_significant_sensitiser_shRNAs_PDS = pds_dict[gene_name][1]
        logFC_significant_sensitiser_shRNAs_median_PDS = pds_dict[gene_name][2]
        n_significant_sensitiser_shRNAs_PhenDC3 = phendc3_dict[gene_name][1]
        logFC_significant_sensitiser_shRNAs_median_PhenDC3 = phendc3_dict[gene_name][2]
      else:
        gene_name = ";".join(list(inter))
        n_detected_shRNAs = pds_dict[list(inter)[1]][0]
        n_significant_sensitiser_shRNAs_PDS = pds_dict[list(inter)[1]][1]
        logFC_significant_sensitiser_shRNAs_median_PDS = pds_dict[list(inter)[1]][2]
        n_significant_sensitiser_shRNAs_PhenDC3 = phendc3_dict[list(inter)[1]][1]
        logFC_significant_sensitiser_shRNAs_median_PhenDC3 = phendc3_dict[list(inter)[1]][2]
#      print gene_name, withoutFC, withFC, n_detected_shRNAs, n_significant_sensitiser_shRNAs_PDS, logFC_significant_sensitiser_shRNAs_median_PDS, n_significant_sensitiser_shRNAs_PhenDC3, logFC_significant_sensitiser_shRNAs_median_PhenDC3, rscore, ";".join(list(set(synonyms2 + synonyms)))
      ofile.write("\t".join([gene_name, withoutFC, withFC, n_detected_shRNAs, n_significant_sensitiser_shRNAs_PDS, logFC_significant_sensitiser_shRNAs_median_PDS, n_significant_sensitiser_shRNAs_PhenDC3, logFC_significant_sensitiser_shRNAs_median_PhenDC3, str(rscore)]) + "\n")


ofile.close()

```

Further analysis:

```bash
cd /lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20170727/polysearch2
wc -l 20170810_polysearch2_1501863075_quadruplex_pubmed_pmc_screen_summary.txt # 527
```

- 633 hits:
  - 51 gene families
  - 582 genes
    - 526 mapped to screen gene names
    - 56 that do not map (most likely from other species) e.g.


```bash
cd /lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20170727/polysearch2
cut 20170810_polysearch2_1501863075_quadruplex_pubmed_pmc_screen_summary.txt -f2 | sort | uniq -c
#    468 no
#      1 without FC threshold?
#     58 yes
```




## Write methods for literature / uniprot go and send to Darcie/Katie tomorrow

This is a draft to add to the computational methods section in the paper.

All code and scripts will be shared in a group's github [page](https://github.com/sblab-bioinformatics) dedicated to the paper, which we will need to create soon. The raw data also needs to be submitted to public repositories like [ArrayExpress](http://www.ebi.ac.uk/arrayexpress/) so it would be good to talk about this too.

**Exploring genes associated to G-quadruplexes in databases and biomedical literature**
We developed two different approaches to discover genes linked to G-quadruplexes:

*Scanning UniprotKB and Gene Ontology (GO) databases*
A list of 18 Homo sapiens genes associated to G-quadruplexes was obtained from UniprotKB (The UniProt Consortium, 2017) using searches with the following criteria:

- Genes annotated with at least one of the following eleven GO terms with any evidence assertion method (Ashburner et al., 2000):

GO:0051880 G-quadruplex DNA binding https://www.ebi.ac.uk/QuickGO/term/GO:0051880 (function)
GO:0002151 G-quadruplex RNA binding https://www.ebi.ac.uk/QuickGO/term/GO:0002151 (function)
GO:0061849 telomeric G-quadruplex DNA binding https://www.ebi.ac.uk/QuickGO/term/GO:0061849 (function)
GO:0071919 G-quadruplex DNA formation https://www.ebi.ac.uk/QuickGO/term/GO:0071919 (biological process)
GO:0044806 G-quadruplex DNA unwinding https://www.ebi.ac.uk/QuickGO/term/GO:0044806 (biological process)
GO:1905493 regulation of G-quadruplex DNA binding https://www.ebi.ac.uk/QuickGO/term/GO:1905493 (biological process)
GO:1905494 negative regulation of G-quadruplex DNA binding https://www.ebi.ac.uk/QuickGO/term/GO:1905494 (biological process)
GO:1905495 positive regulation of G-quadruplex DNA binding https://www.ebi.ac.uk/QuickGO/term/GO:1905495 (biological process)
GO:1905465 regulation of G-quadruplex DNA unwinding https://www.ebi.ac.uk/QuickGO/term/GO:1905465 (biological process)
GO:1905466 negative regulation of G-quadruplex DNA unwinding https://www.ebi.ac.uk/QuickGO/term/GO:1905466 (biological process)
GO:1905467 positive regulation of G-quadruplex DNA unwinding https://www.ebi.ac.uk/QuickGO/term/GO:1905467 (biological process)

- Genes where the corresponding UniprotKB entry is annotated with the term 'quadruplex'

From the total of eighteen genes, five of them (ATRX, DNA2, XRN1, DHX36, MCRS1) were found to be sensitisers in the genome-wide screen after PDS or PhenDC3 treatments (20170810_uniprot_go_G-quadruplex_or_quadruplex_term_homosapiens_screen_summary.txt).


*Mining PubMed and PubMed Central*
In order to expand the list of G-quadruplex-associated genes found in UniprotKB and GO, we explored associations between G-quadruplex terms and gene names in the biomedical literature using the text-mining tool PolySearch2 (Liu, 2015). This algorithm assumes that the greater the co-occurence frequency of terms within sentences or database records, the stronger the association is. The word span between co-occuring terms in the text also influences the association score.

The set of G-quadruplex entities and synonyms was defined using the corresponding MeSH term id [D054856](https://www.ncbi.nlm.nih.gov/mesh/68054856) and the thesaurus of gene names obtained from the PolySearch2 [website](http://polysearch.cs.ualberta.ca/). A total of 5477 pieces of text were identified in PubMed and PubMed Central where any of the G-quadruplex terms co-occur with more than 500 homo sapiens gene names, 58 of them were also found to be sensitisers in the genome-wide screen after PDS or PhenDC3 (20170810_polysearch2_1501863075_quadruplex_pubmed_pmc_screen_summary.txt).

Given that database records contain more carefully curated knowledge in comparison to free-text articles, we considered the UniprotKB and GO associations to be of higher quality than the PubMed searches.





## Why BRCA1 (PS07570) and BRCA2 (PS07592) are not present in the output table?

Especially because there are two papers:

- https://www.ncbi.nlm.nih.gov/pubmed/?term=26748828 > linked to "ADP ribosyl transferase" PARP (PS06071), TP53BP1 (PS12880), MAD2 beta (PS01994) rather than BRCA1

http://polysearch.cs.ualberta.ca/moredetails/1501863075?eindex=43&ename=ADP+ribosyl+transferase&query_term=Quadruplex&link=http%3A%2F%2Fwww.ncbi.nlm.nih.gov%2Fpubmed%2F26748828&eid=PS06071&etype=gene&pmid=26748828&ref=+%282016%29+Targeting+BRCA1+and+BRCA2+Deficiencies+with+G-Quadruplex-Interacting+Compounds.+Molecular+cell%3BMol.+Cell%3B2016+Feb%3B61%283%29%3A449-60

**Targeting BRCA1 and BRCA2 Deficiencies with G-Quadruplex-Interacting Compounds. Molecular cell;Mol. Cell;2016 Feb;61(3):449-60**
G-Quadruplex (G4) -forming genomic sequences, including telomeres, represent natural replication fork barriers. Stalled replication forks can be stabilized and restarted by homologous recombination (HR), which also repairs DNA double-strand breaks (DSBs) arising at collapsed forks. We have previously shown that HR facilitates telomere replication. Here, we demonstrate that the replication efficiency of guanine-rich (G-rich) telomeric repeats is decreased significantly in cells lacking HR. Treatment with theA G4-stabilizing compound pyridostatin (PDS) increases telomere fragility in BRCA2-deficient cells, suggesting that G4 formation drives telomere instability. Remarkably, PDS reduces proliferation of HR-defective cells by inducing DSB accumulation, checkpoint activation, and deregulated G2/M progression and by enhancing the replication defect intrinsic to HR deficiency. PDS toxicity extends to HR-defective cells that have acquired olaparib resistance through loss of 53BP1 or REV7. Altogether, these results highlight the therapeutic potential of G4-stabilizing drugs to selectively eliminate HR-compromised cells and tumors, including those resistant to PARP inhibition.

- https://www.ncbi.nlm.nih.gov/pubmed/?term=28211448 -> linked to PARP (PS06071) rather than BRCA1

http://polysearch.cs.ualberta.ca/moredetails/1501863075?eindex=43&ename=ADP+ribosyl+transferase&query_term=Quadruplex&link=http%3A%2F%2Fwww.ncbi.nlm.nih.gov%2Fpubmed%2F28211448&eid=PS06071&etype=gene&pmid=28211448&ref=+%282017%29+CX-5461+is+a+DNA+G-quadruplex+stabilizer+with+selective+lethality+in+BRCA1%2F2+deficient+tumours.+Nature+communications%3BNat+Commun%3B2017+Feb%3B8%3A14432

**CX-5461 is a DNA G-quadruplex stabilizer with selective lethality in BRCA1/2 deficient tumours. Nature communications;Nat Commun;2017 Feb;8:14432**
G-Quadruplex DNAs form four-stranded helical structures and are proposed to play key roles in different cellular processes. Targeting G-Quadruplex DNAs for cancer treatment is a very promising prospect. Here, we show that CX-5461 is a G-Quadruplex stabilizer, with specific toxicity against BRCA deficiencies in cancer cells and polyclonal patient-derived xenograft models, including tumours resistant to PARP inhibition. Exposure to CX-5461, and its related drug CX-3543, blocks replication forks and induces ssDNA gaps or breaks. The BRCA and NHEJ pathways are required for the repair of CX-5461 and CX-3543-induced DNA damage and failure to do so leads to lethality. These data strengthen the concept of G4 targeting as a therapeutic approach, specifically for targeting HR and NHEJ deficient cancers and other tumours deficient for DNA damage repair. CX-5461 is now in advanced phase I clinical trial for patients with BRCA1/2 deficient tumours (Canadian trial, NCT02719977, opened May 2016).






## There is something wrong with the genes that have two names, e.g. RTEL1, the logFC ... are not correct






## Update PolySearch2 above but including pubmed ids of papers associated to each link found

- ifile2["hits"]["gene"][0] -> first entity id
- ifile2["hits"]["gene"][0]["evidence"][0] -> first entity id, medline
- ifile2["hits"]["gene"][0]["evidence"][1] -> first entity id, pmc
- ifile2["hits"]["gene"][0]["evidence"][0]["entries"][0] -> first entity id, medline, first paper
- ifile2["hits"]["gene"][0]["evidence"][0]["entries"][0]["doc_id"] -> first entity id, medline, first paper and its pubmed id
- ifile2["hits"]["gene"][0]["evidence"][0]["entries"][0]["rscore"] -> first entity id, medline, first paper and its rscore
- ifile2["hits"]["gene"][0]["evidence"][0]["entries"][0]["snippets"] -> first entity id, medline, first paper and its abstract annotated
- ifile2["hits"]["gene"][0]["evidence"][1]["entries"][0]["doc_id"] -> first entity id, pmc, first paper and its pubmed id

I am going to use Search ID: **1501863075**. Exploring POLYSEARCH-1501863075.json:

- doc_count:
  - medline: 4233
  - pmc: 1244
- found_hits: 633
- num_hits: 5477

```python
# Loading screen lists of genes
ifile1 = open("/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20170727/tables/List_for_literature_search.csv", "r")
ilines1 = ifile1.read().split("\r")
ifile1.close()

list_withoutFC = []
list_withFC = []

for line in ilines1[1:]:
  fields = line.split(",")
  if fields[0] != "":
    list_withoutFC.append(fields[0])
  if fields[2] != "":
    list_withFC.append(fields[2])


len(list_withoutFC) # 843 without FC threshold genes
len(set(list_withoutFC)) # 843 unique without FC threshold genes
len(list_withFC) # 763 with FC threshold genes
len(set(list_withFC)) # 763 unique with FC threshold genes


# Loading thesauri into dictionaries
## POLYSEARCH2_THESAURUS_genes.txt
t_genes_ifile = open("/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20170727/polysearch2/POLYSEARCH2_THESAURUS_genes.txt", "r")
t_genes_ilines = t_genes_ifile.readlines()
t_genes_ifile.close()

t_genes = {}

for line in t_genes_ilines:
  fields = line.split("\t")
  entity = fields[0]
  t_genes[entity] = []
  for name in fields[1:]:
    name = name.replace("\n", "")
    t_genes[entity].append(name)


len(t_genes) # 27994
t_genes["PS12860"] # ['DNA topoisomerase 1', 'DNA topoisomerase I', 'TOP-1', 'TOP1', 'TOP1 protein', 'topoisomerase (DNA) I', 'DNA topoisomerase Is', 'TOP1 proteins', 'topoisomerase (DNA) Is', 'Type I DNA topoisomerase', 'EC 5.99.1.2']

## POLYSEARCH2_THESAURUS_gene_family.txt
t_genefamily_ifile = open("/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20170727/polysearch2/POLYSEARCH2_THESAURUS_gene_family.txt", "r")
t_genefamily_ilines = t_genefamily_ifile.readlines()
t_genefamily_ifile.close()

t_genefamily = {}

for line in t_genefamily_ilines:
  fields = line.split("\t")
  entity = fields[0]
  t_genefamily[entity] = []
  for name in fields[1:]:
    name = name.replace("\n", "")
    t_genefamily[entity].append(name)


len(t_genefamily) # 404
t_genefamily["PSF00377"] # ['Telomerase', 'Telomerases']


# Loading PDS and PhenDC3 tables
ifile_pds = open("/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20170727/tables/20170130_PDSvst0_genes.txt", "r")
ilines_pds = ifile_pds.readlines()
ifile_pds.close()

pds_dict = {}

for line in ilines_pds[1:]:
  fields = line.split()
  gene = fields[0]
  n_detected_shRNAs = fields[2]
  n_significant_sensitiser_shRNAs = fields[11]
  logFC_significant_sensitiser_shRNAs_median = fields[14]
  pds_dict[gene] = [n_detected_shRNAs, n_significant_sensitiser_shRNAs, logFC_significant_sensitiser_shRNAs_median]


pds_dict["TOP1"] # ['8', '3', '-3.12857659580522']
pds_dict["MYC"] # ['8', '6', '-2.99796523521497']
pds_dict["ATRX"] # ['16', '1', '-1.55023350217559']
pds_dict["RTEL1"] # ['4', '1', '-3.02109555091145']   RTEL1|TNFRSF6B
pds_dict["TNFRSF6B"] # ['2', '0', 'NA']   RTEL1|TNFRSF6B
pds_dict["BRIP1"] # ['7', '0', 'NA']   BRIP1;BACH1
pds_dict["BACH1"] # ['5', '0', 'NA']   BRIP1;BACH1


ifile_phendc3 = open("/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20170727/tables/20170130_PhenDC3vst0_genes.txt", "r")
ilines_phendc3 = ifile_phendc3.readlines()
ifile_phendc3.close()

phendc3_dict = {}

for line in ilines_phendc3[1:]:
  fields = line.split()
  gene = fields[0]
  n_detected_shRNAs = fields[2]
  n_significant_sensitiser_shRNAs = fields[11]
  logFC_significant_sensitiser_shRNAs_median = fields[14]
  phendc3_dict[gene] = [n_detected_shRNAs, n_significant_sensitiser_shRNAs, logFC_significant_sensitiser_shRNAs_median]


phendc3_dict["TOP1"] # ['8', '4', '-3.92126848474496']
phendc3_dict["MYC"] # ['8', '0', 'NA']
phendc3_dict["ATRX"] # ['16', '5', '-1.39558644062454']
phendc3_dict["RTEL1"] # ['4', '2', '-4.77796076453467']   RTEL1|TNFRSF6B
phendc3_dict["TNFRSF6B"] # ['2', '0', 'NA']   RTEL1|TNFRSF6B
phendc3_dict["BRIP1"] # ['7', '0', 'NA']   BRIP1;BACH1
phendc3_dict["BACH1"] # ['5', '0', 'NA']   BRIP1;BACH1


# Loading POLYSEARCH-1501863075.json and writing output
import json
from operator import itemgetter

ifile2 = json.load(open("/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20170727/polysearch2/POLYSEARCH-1501863075.json"))

all_genes = [line.split()[0] for line in ilines_pds[1:]]

polysearch2 = []

for entity in ifile2["hits"]["gene"]:
  i = entity["eid"]
  n = str(entity["ename"])
  rscore = str(entity["rscore"])
  synonyms = [str(syn) for syn in entity["esynonyms"]]
  evidence = entity['evidence']
  if not i.startswith("PSF"):
    synonyms2 = t_genes[i]
    inter = set(synonyms + synonyms2).intersection(set(all_genes))
    if len(inter) > 0:
      if len(set(synonyms + synonyms2).intersection(set(list_withoutFC))) > 0:
        withoutFC = "yes"
      else:
        withoutFC = "no"
      if len(set(synonyms + synonyms2).intersection(set(list_withFC))) > 0:
        withFC = "yes"
      else:
        withFC = "no"
      pubs = []
      for e in evidence:
        for entry in e['entries']:
          pubs = pubs + [(str(entry['doc_id']), str(entry['doc_type']), int(entry['rscore']))]
      pubs.sort(key=itemgetter(2), reverse=True)
      pubs2 = [",".join(list((p[0], p[1], str(p[2])))) for p in pubs]
      polysearch2.append(["|".join(list(inter)), n, withoutFC, withFC, rscore, "|".join(pubs2), "|".join(list(set(synonyms + synonyms2))).replace("\r", " ")])


len(polysearch2) # 526

ofile = open("/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20170727/polysearch2/20171013_polysearch2_1501863075_quadruplex_medline_pmc_screen_summary.txt", "w")

ofile.write("Gene name screen\tGene name polysearch2\twithout FC threshold?\twith FC threshold?\tRelevancy Score\tPubmed ids\tSynonyms polysearch2\n")

for gene in polysearch2:
  ofile.write("\t".join(gene) + "\n")

ofile.close()

```

Further analysis:

```bash
cd /lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20170727/polysearch2
wc -l 20171013_polysearch2_1501863075_quadruplex_medline_pmc_screen_summary.txt # 527
cut 20171013_polysearch2_1501863075_quadruplex_medline_pmc_screen_summary.txt -f3-4 | sort | uniq -c
#    468 no	no
#      4 yes	no
#     54 yes	yes

```





## 20170220 PolySearch2

Request from Katie, what happens with the following papers? Pubmed ids:
  - 29129743: TLN1, RPS29, hRNPU, GRSF1 - Pubmed abstract but not in PMC, association with these genes are only present in the main text, associations with other genes are present in the abstract
  - 29434328: FUS - Both Pubmed abstract and PMC, association is present in the abstract
  - 28760773: TIMELESS - Pubmed abstract but not in PMC, association is present in the abstract but involving G-Quadruplexes and TIMELESS/TIPIN


- Search ID: **1519121447** (2547)
  - Given: Text
  - Query Keyword: Quadruplex
  - Find ALL associated: Genes/Proteins
  - Automated synonym list: quadruplex;G-Quadruplex;G-quadruplex;Quadruplexes;quadruplexes;G-Quadruplexes;G-quadruplexes;Tetrad;tetrad;Tetrads;tetrads;Quartet;quartet;Quartets;quartets;Tetraplex;tetraplex;Tetraplexes;tetraplexes
  - Custom Negation Words, separate words with ";": (blank)
  - Select one or more Corpus to search: PubMed, PubMed Central
  - Document Limit: No Limit
  - Show borderline associations: (ticked)

- Search ID: **1519125198** (???) - same parameters as 1519121447, just to double-check
  - Given: Text
  - Query Keyword: Quadruplex
  - Find ALL associated: Genes/Proteins
  - Automated synonym list: quadruplex;G-Quadruplex;G-quadruplex;Quadruplexes;quadruplexes;G-Quadruplexes;G-quadruplexes;Tetrad;tetrad;Tetrads;tetrads;Quartet;quartet;Quartets;quartets;Tetraplex;tetraplex;Tetraplexes;tetraplexes
  - Custom Negation Words, separate words with ";": (blank)
  - Select one or more Corpus to search: PubMed, PubMed Central
  - Document Limit: No Limit
  - Show borderline associations: (ticked)

- Search ID: **1519923987** (491)
  - Given: Text
  - Query Keyword: Quadruplex
  - Find ALL associated: Genes/Proteins
  - Automated synonym list: quadruplex;G-Quadruplex;G-quadruplex;Quadruplexes;quadruplexes;G-Quadruplexes;G-quadruplexes;Tetrad;tetrad;Tetrads;tetrads;Quartet;quartet;Quartets;quartets;Tetraplex;tetraplex;Tetraplexes;tetraplexes
  - Custom Negation Words, separate words with ";": (blank)
  - Select one or more Corpus to search: **PubMed**
  - Document Limit: No Limit
  - Show borderline associations: (ticked)

- Search ID: **1519129866** (411)
  - Given: Text
  - Query Keyword: **G-quadruplex**
  - Find ALL associated: Genes/Proteins
  - Automated synonym list: **quadruplex;Quadruplex;G-Quadruplex;Quadruplexes;quadruplexes;G-Quadruplexes;G-quadruplexes;Tetrad;tetrad;Tetrads;tetrads;Quartet;quartet;Quartets;quartets;Tetraplex;tetraplex;Tetraplexes;tetraplexes**
  - Custom Negation Words, separate words with ";": (blank)
  - Select one or more Corpus to search: **PubMed**
  - Document Limit: No Limit
  - Show borderline associations: (ticked)

- Search ID: **1519135535** (233)
  - Given: Text
  - Query Keyword: **G-quadruplexes**
  - Find ALL associated: Genes/Proteins
  - Automated synonym list: **G Quadruplexes, DNA; Guanine Quartets; Guanine-Tetrads; Guanine Quadruplexes; Guanine-Quartets; RNA G-Quadruplexes; DNA, Quadruplex; G-Quadruplexes RNA; G Quadruplexes, RNA; DNA G Quadruplexes; Guanine-Quadruplexes; Guanine Tetrads; DNA, Tetraplex; G-Quadruplexes, RNA; RNA, G-Quadruplexes; Tetraplex DNA; G-Quadruplexes, DNA; DNA G-Quadruplexes; Quadruplex DNA; Guanine-Tetrad; RNAs, G-Quadruplexes; RNA, G Quadruplexes; Guanine-Quartet; G-Quadruplexes; G-Quadruplexes RNAs; G Quadruplexes**
  - Custom Negation Words, separate words with ";": (blank)
  - Select one or more Corpus to search: **PubMed**
  - Document Limit: No Limit
  - Show borderline associations: (ticked)

The Automated synonym list appears automatically when introducing G-quadruplexes as query word. Now 28760773 is picked up too but in association with MIDAS, Cyclin E and TERC - not TIMELESS. I think it is because it appears as TIMELESS/TIPIN.

- Search ID: **1519922917** (???)
  - Given: Text
  - Query Keyword: **G-quadruplexes**
  - Find ALL associated: Genes/Proteins
  - Automated synonym list: **G Quadruplexes, DNA; Guanine Quartets; Guanine-Tetrads; Guanine Quadruplexes; Guanine-Quartets; RNA G-Quadruplexes; DNA, Quadruplex; G-Quadruplexes RNA; G Quadruplexes, RNA; DNA G Quadruplexes; Guanine-Quadruplexes; Guanine Tetrads; DNA, Tetraplex; G-Quadruplexes, RNA; RNA, G-Quadruplexes; Tetraplex DNA; G-Quadruplexes, DNA; DNA G-Quadruplexes; Quadruplex DNA; Guanine-Tetrad; RNAs, G-Quadruplexes; RNA, G Quadruplexes; Guanine-Quartet; G-Quadruplexes; G-Quadruplexes RNAs; G Quadruplexes**
  - Custom Negation Words, separate words with ";": (blank)
  - Select one or more Corpus to search: **PubMed** **PubMed Central**
  - Document Limit: No Limit
  - Show borderline associations: (ticked)

- Search ID: **1519128827** (411)
  - Given: Text
  - Query Keyword: **G-quadruplex**
  - Find ALL associated: Genes/Proteins
  - Automated synonym list: **(blank)**
  - Custom Negation Words, separate words with ";": (blank)
  - Select one or more Corpus to search: **PubMed**
  - Document Limit: No Limit
  - Show borderline associations: (ticked)


Just to follow the consistency with previous searches, the idea would be to:

- 1st: 1519121447 (pubmed + pubmed central) (broad search) - the number of hits has increased dramatically compared to 1501863075 above
- 2nd: 1519922917 (pubmed + pubmed central) (stringent search using the MeSH) - by using the MeSH term, we seem to be picking the 28760773 reference requested by Katie, not linked to TIMELESS though but MIDAS, Cyclin E and TERC
- 3rd: 1519135535 (pubmed only) (very stringent search using the MeSH) - by using the MeSH term, we seem to be picking the 28760773 reference requested by Katie, not linked to TIMELESS though but MIDAS, Cyclin E and TERC




### Copy tables

In `clust1-headnode`,

```bash
cd /scratchb/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data
mkdir -p 20180301/polysearch2
cd 20180301
mkdir tables
cd tables
rsync -arvuP martin03@143.65.169.61:/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20170727/tables/List_for_literature_search.csv .
rsync -arvuP martin03@143.65.169.61:/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20170727/tables/20170130_PDSvst0_genes.txt .
rsync -arvuP martin03@143.65.169.61:/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20170727/tables/20170130_PhenDC3vst0_genes.txt .
```

In `C02Q70MUFVH8`,

```bash
cd /Users/martin03/Desktop
rsync -arvuP POLYSEARCH* martin03@10.20.236.34:/scratchb/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20180301/polysearch2/
```

In `clust1-headnode`,

```bash
cd /scratchb/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20180301/polysearch2
ls -lh
#-rw-r--r-- 1 martin03 sblab  41M Feb 20 03:25 POLYSEARCH-1519121447.json
#-rw-r--r-- 1 martin03 sblab  25K Aug  7  2014 POLYSEARCH2_THESAURUS_gene_family.txt
#-rw-r--r-- 1 martin03 sblab 6.6M Aug  7  2014 POLYSEARCH2_THESAURUS_genes.txt
```

The files `POLYSEARCH2_THESAURUS_genes.txt` and `POLYSEARCH2_THESAURUS_gene_family.txt` links ids of genes (PS...) and gene families (PSF...) given in the result `.json` files.



### 1st: 1519121447

- doc_count:
  - medline: 4387
  - pmc: 1244
- found_hits: 2547

```python
# Loading screen lists of genes
ifile1 = open("/scratchb/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20180301/tables/List_for_literature_search.csv", "r")
ilines1 = ifile1.read().split("\r")
ifile1.close()

list_withoutFC = []
list_withFC = []

for line in ilines1[1:]:
  fields = line.split(",")
  if fields[0] != "":
    list_withoutFC.append(fields[0])
  if fields[2] != "":
    list_withFC.append(fields[2])


len(list_withoutFC) # 843 without FC threshold genes
len(set(list_withoutFC)) # 843 unique without FC threshold genes
len(list_withFC) # 763 with FC threshold genes
len(set(list_withFC)) # 763 unique with FC threshold genes


# Loading thesauri into dictionaries
## POLYSEARCH2_THESAURUS_genes.txt
t_genes_ifile = open("/scratchb/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20180301/polysearch2/POLYSEARCH2_THESAURUS_genes.txt", "r")
t_genes_ilines = t_genes_ifile.readlines()
t_genes_ifile.close()

t_genes = {}

for line in t_genes_ilines:
  fields = line.split("\t")
  entity = fields[0]
  t_genes[entity] = []
  for name in fields[1:]:
    name = name.replace("\n", "")
    t_genes[entity].append(name)


len(t_genes) # 27994
t_genes["PS12860"] # ['DNA topoisomerase 1', 'DNA topoisomerase I', 'TOP-1', 'TOP1', 'TOP1 protein', 'topoisomerase (DNA) I', 'DNA topoisomerase Is', 'TOP1 proteins', 'topoisomerase (DNA) Is', 'Type I DNA topoisomerase', 'EC 5.99.1.2']

## POLYSEARCH2_THESAURUS_gene_family.txt
t_genefamily_ifile = open("/scratchb/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20180301/polysearch2/POLYSEARCH2_THESAURUS_gene_family.txt", "r")
t_genefamily_ilines = t_genefamily_ifile.readlines()
t_genefamily_ifile.close()

t_genefamily = {}

for line in t_genefamily_ilines:
  fields = line.split("\t")
  entity = fields[0]
  t_genefamily[entity] = []
  for name in fields[1:]:
    name = name.replace("\n", "")
    t_genefamily[entity].append(name)


len(t_genefamily) # 404
t_genefamily["PSF00377"] # ['Telomerase', 'Telomerases']


# Loading PDS and PhenDC3 tables
ifile_pds = open("/scratchb/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20180301/tables/20170130_PDSvst0_genes.txt", "r")
ilines_pds = ifile_pds.readlines()
ifile_pds.close()

pds_dict = {}

for line in ilines_pds[1:]:
  fields = line.split()
  gene = fields[0]
  n_detected_shRNAs = fields[2]
  n_significant_sensitiser_shRNAs = fields[11]
  logFC_significant_sensitiser_shRNAs_median = fields[14]
  pds_dict[gene] = [n_detected_shRNAs, n_significant_sensitiser_shRNAs, logFC_significant_sensitiser_shRNAs_median]


pds_dict["TOP1"] # ['8', '3', '-3.12857659580522']
pds_dict["MYC"] # ['8', '6', '-2.99796523521497']
pds_dict["ATRX"] # ['16', '1', '-1.55023350217559']
pds_dict["RTEL1"] # ['4', '1', '-3.02109555091145']   RTEL1|TNFRSF6B
pds_dict["TNFRSF6B"] # ['2', '0', 'NA']   RTEL1|TNFRSF6B
pds_dict["BRIP1"] # ['7', '0', 'NA']   BRIP1;BACH1
pds_dict["BACH1"] # ['5', '0', 'NA']   BRIP1;BACH1


ifile_phendc3 = open("/scratchb/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20180301/tables/20170130_PhenDC3vst0_genes.txt", "r")
ilines_phendc3 = ifile_phendc3.readlines()
ifile_phendc3.close()

phendc3_dict = {}

for line in ilines_phendc3[1:]:
  fields = line.split()
  gene = fields[0]
  n_detected_shRNAs = fields[2]
  n_significant_sensitiser_shRNAs = fields[11]
  logFC_significant_sensitiser_shRNAs_median = fields[14]
  phendc3_dict[gene] = [n_detected_shRNAs, n_significant_sensitiser_shRNAs, logFC_significant_sensitiser_shRNAs_median]


phendc3_dict["TOP1"] # ['8', '4', '-3.92126848474496']
phendc3_dict["MYC"] # ['8', '0', 'NA']
phendc3_dict["ATRX"] # ['16', '5', '-1.39558644062454']
phendc3_dict["RTEL1"] # ['4', '2', '-4.77796076453467']   RTEL1|TNFRSF6B
phendc3_dict["TNFRSF6B"] # ['2', '0', 'NA']   RTEL1|TNFRSF6B
phendc3_dict["BRIP1"] # ['7', '0', 'NA']   BRIP1;BACH1
phendc3_dict["BACH1"] # ['5', '0', 'NA']   BRIP1;BACH1


# Loading POLYSEARCH-1519121447.json and writing output
import json
from operator import itemgetter

ifile2 = json.load(open("/scratchb/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20180301/polysearch2/POLYSEARCH-1519121447.json"))

all_genes = [line.split()[0] for line in ilines_pds[1:]]

polysearch2 = []

for entity in ifile2["hits"]["gene"]:
  i = entity["eid"]
  n = str(entity["ename"])
  rscore = str(entity["rscore"])
  synonyms = [str(syn) for syn in entity["esynonyms"]]
  evidence = entity['evidence']
  if not i.startswith("PSF"):
    synonyms2 = t_genes[i]
    inter = set(synonyms + synonyms2).intersection(set(all_genes))
    if len(inter) > 0:
      if len(set(synonyms + synonyms2).intersection(set(list_withoutFC))) > 0:
        withoutFC = "yes"
      else:
        withoutFC = "no"
      if len(set(synonyms + synonyms2).intersection(set(list_withFC))) > 0:
        withFC = "yes"
      else:
        withFC = "no"
      pubs = []
      for e in evidence:
        for entry in e['entries']:
          pubs = pubs + [(str(entry['doc_id']), str(entry['doc_type']), int(entry['rscore']))]
      pubs.sort(key=itemgetter(2), reverse=True)
      pubs2 = [",".join(list((p[0], p[1], str(p[2])))) for p in pubs]
      polysearch2.append(["|".join(list(inter)), n, withoutFC, withFC, rscore, "|".join(pubs2), "|".join(list(set(synonyms + synonyms2))).replace("\r", " ")])


len(polysearch2) # 2365

ofile = open("/scratchb/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20180301/polysearch2/20180301_polysearch2_1519121447_quadruplex_medline_pmc_screen_summary.txt", "w")

ofile.write("Gene name screen\tGene name polysearch2\twithout FC threshold?\twith FC threshold?\tRelevancy Score\tPubmed ids\tSynonyms polysearch2\n")

for gene in polysearch2:
  ofile.write("\t".join(gene) + "\n")

ofile.close()

```

Further analysis:

```bash
cd /scratchb/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20180301/polysearch2
wc -l 20180301_polysearch2_1519121447_quadruplex_medline_pmc_screen_summary.txt # 2366
cut 20180301_polysearch2_1519121447_quadruplex_medline_pmc_screen_summary.txt -f3-4 | sort | uniq -c
#   2147 no	no
#     15 yes	no
#    203 yes	yes
```



### 2nd: 1519922917

Waiting for the polysearch process to finish ...



### 3rd: 1519135535

- doc_count:
  - medline: 1627
- found_hits: 233

```python
# Loading screen lists of genes
ifile1 = open("/scratchb/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20180301/tables/List_for_literature_search.csv", "r")
ilines1 = ifile1.read().split("\r")
ifile1.close()

list_withoutFC = []
list_withFC = []

for line in ilines1[1:]:
  fields = line.split(",")
  if fields[0] != "":
    list_withoutFC.append(fields[0])
  if fields[2] != "":
    list_withFC.append(fields[2])


len(list_withoutFC) # 843 without FC threshold genes
len(set(list_withoutFC)) # 843 unique without FC threshold genes
len(list_withFC) # 763 with FC threshold genes
len(set(list_withFC)) # 763 unique with FC threshold genes


# Loading thesauri into dictionaries
## POLYSEARCH2_THESAURUS_genes.txt
t_genes_ifile = open("/scratchb/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20180301/polysearch2/POLYSEARCH2_THESAURUS_genes.txt", "r")
t_genes_ilines = t_genes_ifile.readlines()
t_genes_ifile.close()

t_genes = {}

for line in t_genes_ilines:
  fields = line.split("\t")
  entity = fields[0]
  t_genes[entity] = []
  for name in fields[1:]:
    name = name.replace("\n", "")
    t_genes[entity].append(name)


len(t_genes) # 27994
t_genes["PS12860"] # ['DNA topoisomerase 1', 'DNA topoisomerase I', 'TOP-1', 'TOP1', 'TOP1 protein', 'topoisomerase (DNA) I', 'DNA topoisomerase Is', 'TOP1 proteins', 'topoisomerase (DNA) Is', 'Type I DNA topoisomerase', 'EC 5.99.1.2']

## POLYSEARCH2_THESAURUS_gene_family.txt
t_genefamily_ifile = open("/scratchb/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20180301/polysearch2/POLYSEARCH2_THESAURUS_gene_family.txt", "r")
t_genefamily_ilines = t_genefamily_ifile.readlines()
t_genefamily_ifile.close()

t_genefamily = {}

for line in t_genefamily_ilines:
  fields = line.split("\t")
  entity = fields[0]
  t_genefamily[entity] = []
  for name in fields[1:]:
    name = name.replace("\n", "")
    t_genefamily[entity].append(name)


len(t_genefamily) # 404
t_genefamily["PSF00377"] # ['Telomerase', 'Telomerases']


# Loading PDS and PhenDC3 tables
ifile_pds = open("/scratchb/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20180301/tables/20170130_PDSvst0_genes.txt", "r")
ilines_pds = ifile_pds.readlines()
ifile_pds.close()

pds_dict = {}

for line in ilines_pds[1:]:
  fields = line.split()
  gene = fields[0]
  n_detected_shRNAs = fields[2]
  n_significant_sensitiser_shRNAs = fields[11]
  logFC_significant_sensitiser_shRNAs_median = fields[14]
  pds_dict[gene] = [n_detected_shRNAs, n_significant_sensitiser_shRNAs, logFC_significant_sensitiser_shRNAs_median]


pds_dict["TOP1"] # ['8', '3', '-3.12857659580522']
pds_dict["MYC"] # ['8', '6', '-2.99796523521497']
pds_dict["ATRX"] # ['16', '1', '-1.55023350217559']
pds_dict["RTEL1"] # ['4', '1', '-3.02109555091145']   RTEL1|TNFRSF6B
pds_dict["TNFRSF6B"] # ['2', '0', 'NA']   RTEL1|TNFRSF6B
pds_dict["BRIP1"] # ['7', '0', 'NA']   BRIP1;BACH1
pds_dict["BACH1"] # ['5', '0', 'NA']   BRIP1;BACH1


ifile_phendc3 = open("/scratchb/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20180301/tables/20170130_PhenDC3vst0_genes.txt", "r")
ilines_phendc3 = ifile_phendc3.readlines()
ifile_phendc3.close()

phendc3_dict = {}

for line in ilines_phendc3[1:]:
  fields = line.split()
  gene = fields[0]
  n_detected_shRNAs = fields[2]
  n_significant_sensitiser_shRNAs = fields[11]
  logFC_significant_sensitiser_shRNAs_median = fields[14]
  phendc3_dict[gene] = [n_detected_shRNAs, n_significant_sensitiser_shRNAs, logFC_significant_sensitiser_shRNAs_median]


phendc3_dict["TOP1"] # ['8', '4', '-3.92126848474496']
phendc3_dict["MYC"] # ['8', '0', 'NA']
phendc3_dict["ATRX"] # ['16', '5', '-1.39558644062454']
phendc3_dict["RTEL1"] # ['4', '2', '-4.77796076453467']   RTEL1|TNFRSF6B
phendc3_dict["TNFRSF6B"] # ['2', '0', 'NA']   RTEL1|TNFRSF6B
phendc3_dict["BRIP1"] # ['7', '0', 'NA']   BRIP1;BACH1
phendc3_dict["BACH1"] # ['5', '0', 'NA']   BRIP1;BACH1


# Loading POLYSEARCH-1519135535.json and writing output
import json
from operator import itemgetter

ifile2 = json.load(open("/scratchb/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20180301/polysearch2/POLYSEARCH-1519135535.json"))

all_genes = [line.split()[0] for line in ilines_pds[1:]]

polysearch2 = []

for entity in ifile2["hits"]["gene"]:
  i = entity["eid"]
  n = str(entity["ename"])
  rscore = str(entity["rscore"])
  synonyms = [str(syn) for syn in entity["esynonyms"]]
  evidence = entity['evidence']
  if not i.startswith("PSF"):
    synonyms2 = t_genes[i]
    inter = set(synonyms + synonyms2).intersection(set(all_genes))
    if len(inter) > 0:
      if len(set(synonyms + synonyms2).intersection(set(list_withoutFC))) > 0:
        withoutFC = "yes"
      else:
        withoutFC = "no"
      if len(set(synonyms + synonyms2).intersection(set(list_withFC))) > 0:
        withFC = "yes"
      else:
        withFC = "no"
      pubs = []
      for e in evidence:
        for entry in e['entries']:
          pubs = pubs + [(str(entry['doc_id']), str(entry['doc_type']), int(entry['rscore']))]
      pubs.sort(key=itemgetter(2), reverse=True)
      pubs2 = [",".join(list((p[0], p[1], str(p[2])))) for p in pubs]
      polysearch2.append(["|".join(list(inter)), n, withoutFC, withFC, rscore, "|".join(pubs2), "|".join(list(set(synonyms + synonyms2))).replace("\r", " ")])


len(polysearch2) # 211

ofile = open("/scratchb/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20180301/polysearch2/20180301_polysearch2_1519135535_quadruplex_medline_pmc_screen_summary.txt", "w")

ofile.write("Gene name screen\tGene name polysearch2\twithout FC threshold?\twith FC threshold?\tRelevancy Score\tPubmed ids\tSynonyms polysearch2\n")

for gene in polysearch2:
  ofile.write("\t".join(gene) + "\n")

ofile.close()
```

Further analysis:

```bash
cd /scratchb/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20180301/polysearch2
wc -l 20180301_polysearch2_1519135535_quadruplex_medline_pmc_screen_summary.txt # 212
cut 20180301_polysearch2_1519135535_quadruplex_medline_pmc_screen_summary.txt -f3-4 | sort | uniq -c
#    184 no	no
#      4 yes	no
#     23 yes	yes
```





## TODO

- Develop Europe PMC API engine
