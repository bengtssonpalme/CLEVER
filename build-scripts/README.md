## CLEVER BUILD INSTRUCTIONS ##

Note that this is mostly a guide for internal use; in most cases there are little motivation to deploy a copy of CLEVER in your own space (as long as it is being actively maintained by the Bengtsson-Palme lab and hosted at https://microbiology.se).

Note also that this building process will only include updates to ResFinder, CARD and ResFinderFG. To add new other sources to the database, new parsing code has to be added to the `annotate_args.pl` and `create_annotation_file.pl` Perl scripts.

### To create a new CLEVER build from an existing version of CLEVER, follow the steps below. ###

1. Update the build.sh script within this folder. Enter the following variables:

  - Find the right *raw* data directory for RESFINDER at bitbucket. The final URL should look like below:
    `RESFINDER='https://bitbucket.org/genomicepidemiology/resfinder_db/raw/eecf0aa207594fe6d51badf808473de62b28cb06/'`
  - Identify the download link for the bz2 archive of the lastest version of CARD. The URL should look like this:
    `CARD='https://card.mcmaster.ca/download/0/broadstreet-v4.0.1.tar.bz2'`
  - Find the link to the ResFinderFG amino acid sequences. The final URL should look like below:
    `RESFINDERFG='https://raw.githubusercontent.com/RemiGSC/ResFinder_FG_Construction/refs/heads/main/output/RFG_db/ResFinder_FG_AA.faa'`
  - Find the download link to the latest version of PLSDB. The URL should look like this:
    `PLSDB='https://ccb-microbe.cs.uni-saarland.de/plsdb2025/download_fasta'`
  - Change the version to the *new* version of CLEVER (increment by one, or a subversion). For example:
    `VERSION=1.1`
  - Make sure the correct annotation file for the manually curated annotations for the *previous version* of CLEVER is named:
    `CLEVER_ANNOTATION='CLEVER_annotation_curated_2026-03-18'`

2. Run the updated build script using bash. Make sure that you are starting the script from within the CLEVER directory, so that the `build-scripts` directory is directly under it. It should then be possible to run the build script like this:
   
   `bash build-scripts/build.sh`

If there is a problem and you need to rerun the script, but you do not need to re-download all of the databases, it is possible to add the option 'ONLYBUILD' to skip all the downloads, e.g.:

  `bash build-scripts/build.sh ONLYBUILD`

3. The build script might request some user curation for cases where it cannot figure out the correct annotation itself. Provide this as necessary (following the on-screen instructions).

4. The script will (if successful) generate a file called `CLEVER_annotation.txt`. A copy of this file can be opened in e.g. Excel to manually curate CLEVER entries. Sequences marked with an asterisk `*` in the last column are cases which require some user input to resolve them. Save the copy after finishing manual curation (in this example it would be saved under the name "CLEVER_annotation_curated.txt").

5. Run the `update_database_post_curation.pl` script. Make sure that you specify the right CLEVER-build directory, as some of the files will be overwritten based on the curated data! Run the script like this:

  `perl update_database_post_curation.pl CLEVER_annotation.txt CLEVER_annotation_curated.txt CLEVER-build`

6. Once this step is completed (which may take a few minutes), the new CLEVER build is ready. The final directory will contain the files:
 
   - CLEVER.families.faa : All CLEVER families (90% identity clusters) in protein FASTA format – great for annotation and e.g. Diamond searches
   - CLEVER.families.tsv : Annotation for CLEVER families in tsv format
   - CLEVER.lineages.faa : All CLEVER lineages (70% identity clusters) in protein FASTA format
   - CLEVER.lineages.tsv : Annotation for CLEVER lineages in tsv format
   - CLEVER.variants.faa : All CLEVER variants (non-redundant proteins) in protein FASTA format - this is the full set of ARG variants in CLEVER
   - CLEVER.variants.tsv : Annotation for CLEVER variants in tsv format
   - CLEVER.variant-family-mapping.txt : A file mapping all variants to a specific gene name in CLEVER

  In addition, the curated annotation file ("CLEVER_annotation_curated.txt" or similar) is essential for building *the next* CLEVER update, so that also needs to be saved.
   
