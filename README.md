# Capstone Project

# Data

+ `data/api` contains data retreived through the LOVD api
+ `data/scraped` contains data retreived through webscraping the LOVD
+ `data/variant_summary.tsv` contains comprehensive ClinVar data, summarized, retreived from ncbi's ftp server
+ `data/GRCh37.gff` is the latest annotation of the GRCh37 asssembly.

# Scripts

- The `get.py` script retrieves genetic variant data from the LOVD database using three main functions:

  + `getDataBS(gene)`: Scrapes XML data with BeautifulSoup, extracts fields, and saves as a TSV file in data/BSData.

  + `getData(gene)`: Parses XML data with xml.etree.ElementTree, extracts fields, and saves as a TSV file in data.

  + `scrape(gene)`: Scrapes paginated HTML tables using pandas.read_html, combines results, logs incomplete pages, and saves as a TSV file in data/scraped.
      + **FINAL ONE USED!!!!**

- The `parse.py` script processes and merges two TSV files `-C` from the CinVar data and `-L` from the LOVD data and merges them. It tries to extract HGVS standard notation for each variation.

  + `format_contig(chrom)`: gets you from a chromosome number into a standardized contig format (16 $\to$ NC_000016.9) based on GCh37.

  + `process_c_file(path, contig, ignored_writer)`:
    + Extracts transcript and transcript HGVS data, assigns unique `CV#####` IDs, and logs invalid rows to an ignored file.
    + Returns mutation relative to transcript in genomic context. 
    + **THIS IS NO LONGER USEFUL!!!**
  + `process_l_file(path, contig, ignored_writer)`:
    + Extracts genomic coordinate of LOVD mutation data, assigns unique `LOVD#####` IDs, and logs invalid rows to an ignored file.
    + Returns mutations relative to chromosome in the hg17/Ch37 assembly. Since the transcripts the mutations are reported on are not always the RefSeq latest version.  

- The `find.sh` script takes a gtf file, coordinate, and chromosome number and returns all features that contain that coordinate on that chromosome.

- The `validateMutalyzer.py` script attempts to noramlize hgvs notation using mutalizer checking if any fail. **CAN NOT USE THIS DUE TO API LIMITS**

**SCRAPING DATA FILES WERE NOT ENTIRELY CONSISTENT**:
+ this affected downstream parsing
+ `head -n 1 data/scraped/* | grep "Effect" | sort | uniq -c` returned 3 different kinds of headers, 5 files were effected
+ fixed with `fixScrape.py`, now all headers are consistent


overlaps and why overlaps

lets focus on the top reported mutations 