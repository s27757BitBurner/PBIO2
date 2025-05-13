#!/usr/bin/env python3
"""
NCBI GenBank Data Retriever with Filtering, CSV Report, and Data Visualization

This script connects to NCBI GenBank using Biopython’s Entrez module to retrieve genetic sequence records
for a given taxonomic ID. The user is prompted for minimum and maximum sequence length thresholds for filtering.
After retrieval, a CSV report is generated containing the GenBank accession number, sequence length, and description.
Additionally, a line chart is produced that charts sequence lengths (sorted from longest to shortest) and saves it as a PNG.

Required Libraries: Biopython, pandas, matplotlib
"""

from Bio import Entrez, SeqIO
import time
import os
import pandas as pd
import matplotlib.pyplot as plt


class NCBIRetriever:
    def __init__(self, email, api_key):
        """Initialize with NCBI credentials and set up Entrez."""
        self.email = email
        self.api_key = api_key
        Entrez.email = email
        Entrez.api_key = api_key
        Entrez.tool = 'BioScriptEx10'

    def search_taxid(self, taxid):
        """
        Search for records with the given taxonomic ID.
        This function first retrieves taxonomy info (such as the organism name) then does an esearch.
        """
        print(f"Searching for records with TaxID: {taxid}")
        try:
            # Get taxonomy info
            handle = Entrez.efetch(db="taxonomy", id=taxid, retmode="xml")
            records = Entrez.read(handle)
            organism_name = records[0]["ScientificName"]
            print(f"Organism: {organism_name} (TaxID: {taxid})")

            # Search for nucleotide records using the taxonomy id
            search_term = f"txid{taxid}[Organism]"
            handle = Entrez.esearch(db="nucleotide", term=search_term, usehistory="y")
            search_results = Entrez.read(handle)
            count = int(search_results["Count"])

            if count == 0:
                print(f"No records found for {organism_name}")
                return None

            print(f"Found {count} records")

            # Save parameters for fetching records in batches
            self.webenv = search_results["WebEnv"]
            self.query_key = search_results["QueryKey"]
            self.count = count

            return count
        except Exception as e:
            print(f"Error during taxonomy search for TaxID {taxid}: {e}")
            return None

    def fetch_records(self, start=0, max_records=10):
        """
        Fetch a batch of GenBank records starting at the given index.
        Returns a list of SeqRecord objects.
        """
        if not hasattr(self, 'webenv') or not hasattr(self, 'query_key'):
            print("No search results available. Please run search_taxid() first.")
            return []
        try:
            # Limit each batch for safety
            batch_size = min(max_records, 500)
            handle = Entrez.efetch(
                db="nucleotide",
                rettype="gb",
                retmode="text",
                retstart=start,
                retmax=batch_size,
                webenv=self.webenv,
                query_key=self.query_key
            )
            # Parse fetched GenBank records
            records = list(SeqIO.parse(handle, "genbank"))
            return records
        except Exception as e:
            print(f"Error fetching records: {e}")
            return []


def retrieve_and_filter_records(retriever, max_records, min_len, max_len):
    """
    Retrieve records in batches (up to max_records) and filter them based on sequence length.
    Returns a list of dictionaries with keys: Accession, Length, Description.
    """
    records_data = []
    fetched_count = 0
    while fetched_count < max_records:
        # Adjust batch_size in case remaining records is less than typical batch size
        batch_size = min(100, max_records - fetched_count)
        records = retriever.fetch_records(start=fetched_count, max_records=batch_size)
        if not records:
            break
        for rec in records:
            seq_length = len(rec.seq)
            if min_len <= seq_length <= max_len:
                # Access the accession number from the record annotations if available;
                # if not, fall back to using rec.id.
                accession = rec.annotations.get("accessions", [rec.id])[0] if rec.annotations.get(
                    "accessions") else rec.id
                description = rec.description
                records_data.append({
                    "Accession": accession,
                    "Length": seq_length,
                    "Description": description
                })
        fetched_count += batch_size
        # Pause if needed to respect NCBI's rate limits (optional)
        # time.sleep(0.4)
    return records_data


def generate_csv_report(records_data, filename):
    """
    Generate a CSV report from records_data containing accession number, sequence length, and description.
    """
    if not records_data:
        print("No records data available to write to CSV.")
        return
    df = pd.DataFrame(records_data)
    df.to_csv(filename, index=False)
    print(f"CSV report generated: {filename}")


def generate_line_chart(records_data, filename):
    """
    Generate a line chart showing sequence lengths (from longest to shortest).
    X-axis: GenBank accession number, Y-axis: sequence length.
    The chart is saved as a PNG image.
    """
    if not records_data:
        print("No records data available for visualization.")
        return

    # Sort records by sequence length in descending order.
    sorted_data = sorted(records_data, key=lambda x: x["Length"], reverse=True)
    accessions = [d["Accession"] for d in sorted_data]
    lengths = [d["Length"] for d in sorted_data]

    plt.figure(figsize=(10, 6))
    plt.plot(accessions, lengths, marker="o", linestyle="-")
    plt.xlabel("GenBank Accession Number")
    plt.ylabel("Sequence Length")
    plt.title("Sequence Lengths of GenBank Records (Longest to Shortest)")
    plt.xticks(rotation=45, ha="right")
    plt.tight_layout()
    plt.savefig(filename)
    plt.close()
    print(f"Line chart saved as: {filename}")


def main():
    print("NCBI GenBank Data Retriever with Filtering, CSV Report, and Data Visualization")

    # Get user credentials and create retriever
    email = input("Enter your email address for NCBI: ")
    api_key = input("Enter your NCBI API key: ")
    retriever = NCBIRetriever(email, api_key)

    # Get taxonomic ID and perform initial search
    taxid = input("Enter taxonomic ID (taxid) of the organism: ")
    count = retriever.search_taxid(taxid)
    if not count:
        print("No records found. Exiting.")
        return

    # Get sequence length filtering parameters from the user
    try:
        min_seq_length = int(input("Enter minimum sequence length: "))
        max_seq_length = int(input("Enter maximum sequence length: "))
    except ValueError:
        print("Invalid input for sequence length thresholds. Exiting.")
        return

    # Ask for maximum number of records to retrieve. If 0 is entered, fetch all available records.
    try:
        max_records_input = int(input("Enter maximum number of records to retrieve (0 for all available): "))
        if max_records_input == 0:
            max_records = count
        else:
            max_records = max_records_input
    except ValueError:
        print("Invalid input for maximum record count. Defaulting to 10 records.")
        max_records = 10

    print("\nRetrieving and filtering records based on sequence length criteria...")
    records_data = retrieve_and_filter_records(retriever, max_records, min_seq_length, max_seq_length)
    print(f"After filtering, {len(records_data)} records matched the criteria.")

    if not records_data:
        print("No records matched the filter criteria. Exiting.")
        return

    # Generate CSV report
    csv_filename = f"taxid_{taxid}_report.csv"
    generate_csv_report(records_data, csv_filename)

    # Generate data visualization (line chart)
    chart_filename = f"taxid_{taxid}_chart.png"
    generate_line_chart(records_data, chart_filename)


if __name__ == "__main__":
    main()
