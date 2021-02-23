from gpcrdb import GPCRdb
import json


def test_fetch_generic_numbers():
    entries = GPCRdb.fetch_generic_numbers("opsd_bovin")
    print(json.dumps(entries, indent=2))

    assert len(entries) > 0

    assert "sequence_number" in entries[0]
    assert "amino_acid" in entries[0]
    assert "protein_segment" in entries[0]
    assert "display_generic_number" in entries[0]

    res33 = [r for r in entries if r["sequence_number"] == 33]
    assert len(res33) == 1
    res33 = res33[0]

    assert res33["protein_segment"] == "TM1"
    assert res33["display_generic_number"] == "1.28x28"


def test_fetch_protein():
    entry = "opsd_bovin"
    prot = GPCRdb.fetch_protein(entry)

    assert prot["entry_name"] == entry
    assert prot["name"] == "Rhodopsin"
    assert prot["accession"] == "P02699"
    assert prot["family"] == "001_009_001_001"
    assert len(prot["sequence"]) > 0
    assert prot["genes"] == ["RHO"]
