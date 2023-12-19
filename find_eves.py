"""
please find attached a file, seqs.fasta, containing a series of exogenous and endogenous viruses. 
We would like you to write a script/programme, ideally in Python but any language of your choice, 
that would do the following: 
(1) determine which sequences are known exogenous viruses and which are EVEs, 
(2) approximately classify the EVEs based on the exogenous virus that they are most similar to 
(in terms of taxonomy and identity of the gene).
"""

from Bio.Blast import NCBIWWW as ncbi
from Bio.Blast import NCBIXML as xml
from Bio import SeqIO, Entrez

import argparse
import csv
import time
import sys

def get_args(parser):
    parser.add_argument('file', 
        help = 'Path to fasta file with query sequences'),
    parser.add_argument('out', 
        help = 'Label for outputs'), 
    parser.add_argument("--email", type=str, default="lmuoz@uwo.ca",
        help="<option> Required for Entrez database transactions."), 
    parser.add_argument("--coverage", type=int, default=80,
        help="<option> Minimum coverage to be considered exogenous virus."),
    parser.add_argument("--identity", type=int, default=60,
        help="<option> Minimum identity to be considered exogenous.")

    return parser.parse_args()

def blast_n(query):
    """
    Nucleotide blast against viruses
    :param query: str, nucleotide sequence
    :return hit: alignments from first hit
    """
    response = ncbi.qblast("blastn", "nt", query.seq, 
                           hitlist_size=10, alignments = 10, 
                           entrez_query='txid10239[ORGN]')
    
    hits = xml.parse(response)
    hit = next(hits)
    time.sleep(10)  # wait before next attempt
    return hit

def parse_aln(aln):
    """
    Extract statistics from alignment
    Use only the first scoring pair
    :param hit: Bio.Blast.Record parsed from xml response
    :param record_length: int, to calculate hit coverage
    :return results: list, alignment stats
    """
    accn = aln.accession
    desc = aln.hit_def
    h0 = aln.hsps[0]

    # extract statistics
    alen = h0.align_length
    ascore = h0.score
    evalue = h0.expect
    per_ident = h0.identities / alen * 100
    gaps = h0.gaps / alen * 100
    strand = '/'.join(h0.strand)
    stats = [accn, alen, per_ident, ascore, evalue, gaps, strand, desc]
    
    return stats

def get_hit_metadata(accn, email):
    """
    Extract metadata from an accession number
    :param accn: str, accession number to efetch
    :param email: str, required for Entrez database transactions
    :return records: parsed records
    """
    Entrez.email = email
    print(f"Getting info of hit against: {accn}")
    handle = Entrez.efetch(db='nuccore', id=accn, rettype='gb', 
                           retmode='text')
    
    records = SeqIO.parse(handle, 'gb')
    return records

def get_record_info(record):
    """
    Extract taxonomic information associated with a genbank record
    :return record_info: list, record ID, host, and taxonomic information
    """
    taxon = ':'.join(record.annotations["taxonomy"])
    source = list(filter(lambda f: f.type=='source', record.features))
    quals = source[0].qualifiers
    record_info = [
        record.id,  # record.name is LOCUS, not always accession
        quals.get('host', [''])[0],
        taxon,
        # quals.get('isolate', [''])[0],
        # quals.get('country', [''])[0],
        # quals.get('collection_date', [''])[0]
        ]
    
    return record_info

if __name__ == "__main__":
    parser=argparse.ArgumentParser(
        description='Classify sequences as exogenous or endogenous viruses'
    )

    args = get_args(parser)
    email = args.email
    queries = SeqIO.parse(args.file, "fasta")

    virus = open(f'{args.out}_viruses.csv', 'w')
    virus_w = csv.writer(virus)
    virus_w.writerow(['q_name', 'coverage', 'hit_accn', 'align_len', 
                      'per_ident', 'ascore', 'evalue', 'gaps', 
                      'strand', 'desc'])
    
    eves = open(f'{args.out}_eves.csv', 'w')
    eves_w = csv.writer(eves)
    eves_w.writerow(['q_name', 'coverage', 'hit_accn', 'align_len', 
                      'per_ident', 'ascore', 'evalue', 'gaps', 
                      'strand', 'desc', 'hit_id', 'host', 'taxonomy'])
    # sys.exit()
    for query in queries:
        name = query.name
        print(f'\nBlasting: {name}\n')
        hit = blast_n(query)

        for aln in hit.alignments:
            res = parse_aln(aln)  # alignment stats
            coverage = res[1] / len(query.seq) * 100
            identity = res[2]
            accn = res[0]
        
            # Alignment covers most of the sequence, with high identity
            if coverage > args.coverage and identity > args.identity:
                # Likely a virus
                row = [name, coverage] + res
                virus_w.writerow(row)        
                break  # Only get the first hit

            # Low coverage or low identity, is it EVE?
            else:
                records = get_hit_metadata(accn, email)
                for record in records:
                    info = get_record_info(record)
                    row = [name, coverage] + res + info
                    eves_w.writerow(row)

    eves.close()
    virus.close()