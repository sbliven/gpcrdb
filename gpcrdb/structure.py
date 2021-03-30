"""Working with protein structures

This module requires biopython
"""

import warnings
from io import StringIO
from typing import TYPE_CHECKING, Dict, Tuple

import requests
from Bio import BiopythonWarning  # type: ignore
from Bio.PDB import Dice  # type: ignore
from Bio.PDB.PDBIO import PDBIO  # type: ignore
from Bio.PDB.PDBParser import PDBParser  # type: ignore

from .gpcrdb import BASE_URL

if TYPE_CHECKING:
    from Bio.PDB.Structure import Structure   # type: ignore


class FullChainSelector(Dice.ChainSelector):
    "Extracts a full chain from a structure"

    def __init__(self, chain_id, model_id=0):
        """Initialize the class."""
        super().__init__(chain_id, None, None, model_id)

    def accept_residue(self, residue):
        """Verify if a residue sequence is between the start and end sequence."""
        # residue - between start and end
        hetatm_flag, resseq, icode = residue.get_id()
        if hetatm_flag != " ":
            # skip HETATMS
            return 0
        if icode != " ":
            warnings.warn(
                "WARNING: Icode %s at position %s" % (icode, resseq), BiopythonWarning
            )
        return 1


def assign_generic_numbers(structure: "Structure", chain=None) -> Dict[Tuple, float]:
    """Assigns generic residue numbers to an arbitrary BioPython structure"""

    io_pdb = PDBIO()
    io_pdb.set_structure(structure)
    io_str = StringIO()
    if chain is not None:
        sel = FullChainSelector(chain)
        io_pdb.save(io_str, select=sel)
    else:
        io_pdb.save(io_str)

    url = f"{BASE_URL}/structure/assign_generic_numbers"
    content = post(url, io_str.getvalue())

    # numbers are stored in b-values
    parser = PDBParser()
    content_fp = StringIO(content.decode("utf8"))
    content_fp.seek(0)
    new_structure = parser.get_structure(structure.id, content_fp)
    assert len(new_structure) == 1
    chains = list(new_structure.get_chains())
    assert len(chains) == 1

    generics = {res.full_id: res["CA"].bfactor for res in chains[0]}
    return generics


def post(url, filecontents):
    req = requests.post(url, files={"pdb_file": filecontents})
    req.raise_for_status()
    return req.content
