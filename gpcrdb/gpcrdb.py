import json

import requests

# from bisect import bisect_left
# import os
# from Bio import SeqIO
# from Bio.PDB import PDBList, MMCIFParser, PDBParser


BASE_URL = "https://gpcrdb.org/services"


def fetch_generic_numbers(uniprotid):
    """Fetch JSON object with generic numbers for the specified uniprot ID
    (e.g. "opsd_bovin")"""
    url = f"{BASE_URL}/residues/{uniprotid.lower()}"
    return get(url)


def fetch_structures(uniprotid, representative=False):
    """Fetch JSON object with PDB structures for the specified uniprot ID
    (e.g. "opsd_bovin")

    representative: restrict to one structure per activation state
    """
    url = f"{BASE_URL}/structure/protein/{uniprotid.lower()}/"
    if representative:
        url += "representative/"
    return get(url)


def fetch_protein(uniprotid):
    """Get general information about a protein."""
    url = f"{BASE_URL}/protein/{uniprotid.lower()}/"
    return get(url)


def get(url):
    # Retrive all entries and parse into JSON object
    req = requests.get(url)
    req.raise_for_status()
    entries = json.loads(req.content)

    return entries
