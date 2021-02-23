import json
import requests

# from bisect import bisect_left
# import os
# from Bio import SeqIO
# from Bio.PDB import PDBList, MMCIFParser, PDBParser


class GPCRdb(object):
    BASE_URL = "https://gpcrdb.org/services"

    @classmethod
    def fetch_generic_numbers(cls, uniprotid):
        """Fetch JSON object with generic numbers for the specified uniprot ID
        (e.g. "opsd_bovin")"""
        url = f"{cls.BASE_URL}/residues/{uniprotid.lower()}"
        return cls.get(url)

    @classmethod
    def fetch_structures(cls, uniprotid, representative=False):
        """Fetch JSON object with PDB structures for the specified uniprot ID
        (e.g. "opsd_bovin")

        representative: restrict to one structure per activation state
        """
        url = f"{cls.BASE_URL}/structure/protein/{uniprotid.lower()}/"
        if representative:
            url += "representative/"
        return cls.get(url)

    @classmethod
    def get(cls, url):
        # Retrive all entries and parse into JSON object
        req = requests.get(url)
        req.raise_for_status()
        entries = json.loads(req.content)

        return entries
