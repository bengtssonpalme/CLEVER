# CLEVER
## Version 2.0

This repository contains the code to build and maintain the CLEVER database.

Building instructions can be found in build-scripts/README.md

## Quick introduction to the CLEVER gene header format
This is a CLEVER FASTA header:

`C2|mexF|mex-6_1|ECV|CARD|AAG05882.1`

The fields are delimited by 'bars' (**|**), and have the following meanings:

CLEVER version (C2) | Gene name (mexF) | CLEVER identifier (mex-6_1) | evidence codes (ECV) | source database (CARD) | Identifier in source database (AAG05882.1)

Note that the gene name often is the same as the CLEVER identifier (mex-6_1), but with a wave '**~**' sign added, which means that CLEVER has inferred this name, i.e. there is (as far as we know) no established gene name for this ARG.

So for this gene the gene class is “mex”, the family number is 6 and the variant within the family is 1

Evidence codes are always three letters:

E or L (Established or Latent)
C or M (Chromosomal or Mobile)
V or P or S (Validated or Predicted or Structure-predicted)

The definitions of these concepts are as follows:
- **Established ARG (E)**: An ARG that is experimentally verified to confer antibiotic resistance and is present in human pathogens
- **Latent ARG (L)**: An ARG that confers a resistance function (or is predicted to do so), but does not exist in pathogens
- **Mobile ARG (M)**: An ARG which is present on a mobile genetic element, which could be plasmids, integrons, transposons or integrative conjugative elements.
- **Chromosomal ARG (C)**: Any ARG that does not meet the criteria for a mobile ARG above.
- **Validated ARG (V)**: An ARG for which the resistance function has been verified in laboratory experiments, by showing that the presence of the gene increases the MIC of the host compared to an otherwise isogenic strain that does not carry the gene, or alternatively that over-expression of the gene induces a higher MIC compared to an isogenic reference strain.
- **Predicted ARG (P)**: An ARG for which its function has not been verified experimentally (see above), but has been predicted to be an ARG by fARGene (Berglund et al. 2019)
- **Structure-predicted ARG (S)**: An ARG for which its function has not been verified experimentally (see above), but has been predicted to be an ARG based on 3D-structure similarity by MUSTARD (Ruppé et al. 2018)

## ARG classes covered by CLEVER
| CLEVER Class | Class name                                                       | Resistance profile                     |
|--------------|------------------------------------------------------------------|----------------------------------------|
| aac          | aminoglycoside acetyltransferase (no subclass)                   | aminoglycosides                        |
| aac2p        | aminoglycoside N(2')-acetyltransferase                           | aminoglycosides                        |
| aac3         | aminoglycoside N(3)-acetyltransferase                            | aminoglycosides                        |
| aac6p        | aminoglycoside N(6')-acetyltransferase                           | aminoglycosides                        |
| aad          | aminoglycoside 6-adenylyltransferase                             | aminoglycosides                        |
| adt          | aminoglycoside resistance gene (no class)                        | aminoglycosides                        |
| ant          | aminoglycoside O-nucleotidyltransferase                          | aminoglycosides                        |
| aph          | aminoglycoside O-phosphotransferase (no subclass)                | aminoglycosides                        |
| aph2b        | aminoglycoside 2''-O-phosphotransferase                          | aminoglycosides                        |
| aph3b        | aminoglycoside 3''-O-phosphotransferase                          | aminoglycosides                        |
| aph3p        | aminoglycoside 3'-O-phosphotransferase                           | aminoglycosides                        |
| aph4         | aminoglycoside 4-O-phosphotransferase                            | aminoglycosides                        |
| aph6         | aminoglycoside 6-O-phosphotransferase                            | aminoglycosides                        |
| aph7b        | aminoglycoside 7''-O-phosphotransferase                          | aminoglycosides                        |
| aph9         | aminoglycoside 9-O-phosphotransferase                            | aminoglycosides                        |
| arr          | rifampin ADP-ribosylating transferase                            | rifampin                               |
| blaA         | class A beta-lactamase (no subclass)                             | beta-lactams                           |
| blaA1        | class A1 beta-lactamase                                          | beta-lactams                           |
| blaB1        | class B1 beta-lactamase                                          | beta-lactams                           |
| blaB3        | class B3 beta-lactamase                                          | beta-lactams                           |
| blaC         | class C beta-lactamase (no subclass)                             | beta-lactams                           |
| blaC1        | class C1 beta-lactamase                                          | beta-lactams                           |
| blaC2        | class C2 beta-lactamase                                          | beta-lactams                           |
| blaC3        | class C3 beta-lactamase                                          | beta-lactams                           |
| blaD         | class D beta-lactamase (no subclass)                             | beta-lactams                           |
| blaD1        | class D1 beta-lactamase                                          | beta-lactams                           |
| blaD2        | class D2 beta-lactamase                                          | beta-lactams                           |
| blaX         | beta-lactamase (no class)                                        | beta-lactams                           |
| cat          | chloramphenicol acetyltransferase                                | chloramphenicol                        |
| cfr          | 23S ribosomal RNA methyltransferase                              | macrolides,streptogramines,lincosamide |
| cml          | chloramphenicol exporter                                         | chloramphenicol                        |
| dfr          | trimethoprim-resistant dihydrofolate reductase                   | trimethoprims                          |
| dhr          | trimethoprim-resistant dihydrofolate reductase                   | trimethoprims                          |
| dpr          | sulfamethoxazole resistance gene                                 | sulfonamides                           |
| ere          | erythromycin esterase                                            | macrolides                             |
| erm          | rRNA methyltransferase                                           | macrolides,streptogramines,lincosamide |
| fos          | fosfomycin thiol transferase                                     | fosfomycin                             |
| lnu          | lincosamide nucleotidyltransferase                               | macrolides,streptogramines,lincosamide |
| lsa          | LSA ABC-F subfamily protein                                      | macrolides,streptogramines,lincosamide |
| mcr          | MCR phosphoethanolamine transferase                              | polymyxins                             |
| mec          | methicillin resistant PBP                                        | beta-lactams                           |
| mex          | resistance-nodulation-cell division (RND) antibiotic efflux pump | multidrug-exporter                     |
| mph          | macrolide 2'-phosphotransferase                                  | macrolides,streptogramines,lincosamide |
| msr          | ABC-F binding cassette ribosomal protection protein              | macrolides,streptogramines,lincosamide |
| nim          | nitroimidazole reductase                                         | nitroimidazoles                        |
| omp          | outer membrane porin                                             | beta-lactams                           |
| opr          | outer membrane factor protein                                    | multidrug-exporter                     |
| pbp          | penicillin-binding protein                                       | beta-lactams                           |
| qnr          | quinolone resistance protein (qnr)                               | quinolones                             |
| rmt          | 16S rRNA methyltransferase                                       | aminoglycosides                        |
| sul          | sulfonamide resistant dihydropteroate synthase                   | sulfonamides                           |
| tet          | tetracycline resistance gene (no class)                          | tetracyclines                          |
| tetEFF       | tetracycline efflux pump                                         | tetracyclines                          |
| tetRPG       | tetracycline ribosomal protection gene                           | tetracyclines                          |
| tetX         | tetracycline inactivation enzyme                                 | tetracyclines                          |
| tran         | transporter (no class)                                           | multidrug-exporter                     |
| van          | vancomycin resistance gene                                       | vancomycin                             |
| vat          | streptogramin vat acetyltransferase                              | macrolides,streptogramines,lincosamide |
| vga          | VGA ABC-F subfamily protein                                      | macrolides,streptogramines,lincosamide |
| vgb          | streptogramin vgb lyase                                          | macrolides,streptogramines,lincosamide |
| misc         | miscellaneous ARGs                                               | variable                               |
| unc          | unclassified ARG                                                 | unclassified                           |

## Version 2.0 updates
### Version 2.0 adds new additional sources of fARGene predicted ARGs:
- Inda-Díaz JS, Lund D, Parras-Moltó M, Johnning A, Bengtsson-Palme J, Kristiansson E. Latent antibiotic resistance genes are abundant, diverse, and mobile in human, animal, and environmental microbiomes. Microbiome. 2023;11(1):44. Published 2023 Mar 8. doi:10.1186/s40168-023-01479-0 (already present in CLEVER 1.0)
- Victor MP, Radisic V, Grevskott DH, Marathe NP. Hospital effluent in a low-resistance setting is responsible for dissemination of novel antibiotic resistance genes into the marine environment. Ecotoxicol Environ Saf. 2025;301:118390. doi:10.1016/j.ecoenv.2025.118390
- Victor MP, Øvreås L, Marathe NP. Characterization of known and novel clinically important antibiotic resistance genes and novel microbes from wastewater-impacted high Arctic fjord sediments. Sci Total Environ. 2025;985:179699. doi:10.1016/j.scitotenv.2025.179699
- Li B, Jiang L, Johnson T, et al. Global health risks lurking in livestock resistome. Sci Adv. 2025;11(26):eadt8073. doi:10.1126/sciadv.adt8073 (https://zenodo.org/records/15025586)
- Sommerville V, Meola M, Nunes-Richards A, et al. Microbial community dynamics in a traditional Swiss mountain cheese over 142 years of cheesemaking. bioRxiv, 2026.02.26.708305 (2026). doi:10.64898/2026.02.26.708305

### Version 2.0 also adds additional sources of MUSTARD/PCM predicted ARGs:
- Ruppé E, Ghozlane A, Tap J, et al. Prediction of the intestinal resistome by a three-dimensional structure-based method. Nat Microbiol. 2019;4(1):112-123. doi:10.1038/s41564-018-0292-6
- Wang K, Xu J, Li X, et al. Evolutionary selection of trimethoprim-resistant dfrA genes in lytic phages affects phage and host fitness during infection. Sci Adv. 2025;11(39):eadt4817. doi:10.1126/sciadv.adt4817
