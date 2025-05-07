


import requests
import xml.etree.ElementTree as ET
import pandas as pd
from io import StringIO
import requests
import xml.etree.ElementTree as ET
import pandas as pd
import os
import requests
import pandas as pd
from bs4 import BeautifulSoup
import os

def getDataBS(gene):
    url = f"https://databases.lovd.nl/shared/api/rest.php/variants/{gene}?search_VariantOnTranscript/Haplotype=FA"
    response = requests.get(url)

    soup = BeautifulSoup(response.content, "xml")  # or "html.parser" for max leniency

    data = []
    for entry in soup.find_all("entry"):
        entry_data = {
            "title": entry.find("title").text if entry.find("title") else None,
            "id": entry.find("id").text if entry.find("id") else None,
            "published": entry.find("published").text if entry.find("published") else None,
            "updated": entry.find("updated").text if entry.find("updated") else None
        }

        content = entry.find("content")
        if content and content.text:
            for line in content.text.strip().split("\n"):
                if ":" in line:
                    key, value = line.split(":", 1)
                    entry_data[key.strip()] = value.strip()

        data.append(entry_data)

    df = pd.DataFrame(data)
    os.makedirs("data", exist_ok=True)
    df.to_csv(f"data/BSData/{gene}_variants.tsv", sep="\t", index=False)

    return df


def getData(gene):
    url = f"https://databases.lovd.nl/shared/api/rest.php/variants/{gene}?search_VariantOnTranscript/Haplotype=FA"
    response = requests.get(url)
    root = ET.fromstring(response.content)

    ns = {"atom": "http://www.w3.org/2005/Atom"}

    data = []
    for entry in root.findall("atom:entry", ns):
        entry_data = {
            "title": entry.find("atom:title", ns).text,
            "id": entry.find("atom:id", ns).text,
            "published": entry.find("atom:published", ns).text,
            "updated": entry.find("atom:updated", ns).text
        }

        # Parse key-value lines in the <content>
        content = entry.find("atom:content", ns).text
        if content:
            for line in content.strip().split("\n"):
                if ":" in line:
                    key, value = line.split(":", 1)
                    entry_data[key.strip()] = value.strip()

        data.append(entry_data)

    df = pd.DataFrame(data)

    os.makedirs("data", exist_ok=True)
    df.to_csv(f"data/{gene}_variants.tsv", sep="\t", index=False)

    return df

import pandas as pd

def scrape(gene):
    base_url = f"https://databases.lovd.nl/shared/variants/{gene}/unique"
    all_variants = []
    page = 1

    os.makedirs("data/scraped", exist_ok=True)
    log_file = "data/scraped/incomplete_pages.log"

    while True:
        url = f"{base_url}?page={page}"
        try:
            tables = pd.read_html(url)
            variant_table = tables[8]
        except Exception as e:
            print(f"Failed to load page {page} for {gene}: {e}")
            break

        n_rows = variant_table.shape[0]

        if n_rows == 0:
            break

        print(f"Fetched page {page} with {n_rows} variants")

        #  Log if a non-terminal page has < 100 entries
        if n_rows < 100:
            # Check if this is the last page
            next_url = f"{base_url}?page={page+1}"
            try:
                next_tables = pd.read_html(next_url)
                next_variant_table = next_tables[8]
                if next_variant_table.shape[0] > 0:
                    with open(log_file, "a") as f:
                        f.write(f"{gene}, page {page} only had {n_rows} entries\n")
            except:
                pass  # Assume it's the terminal page if loading fails

        all_variants.append(variant_table)
        page += 1

    # Combine all pages
    df_all = pd.concat(all_variants, ignore_index=True)
    df_all.to_csv(f"data/scraped/{gene}_variants_full.tsv", sep="\t", index=False)


genes = ["FANCA", "FANCB", "FANCC", "BRCA2", "FANCD2", "FANCE",
         "FANCF", "FANCG", "FANCI", "BRIP1", "FANCL", "FANCM",
         "PALB2", "RAD51C", "SLX4", "ERCC4", "RAD51", "BRCA1",
         "UBE2T", "XRCC2", "MAD2L2", "RFWD3"]

for gene in genes:
    try:
        df = scrape(gene)
        print(f"Data for {gene} retrieved successfully.")
    except Exception as e:
        print(f"Error retrieving data for {gene}: {e}")