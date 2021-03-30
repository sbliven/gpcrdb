from Bio.PDB.mmtf import MMTFParser  # type: ignore

from gpcrdb.structure import assign_generic_numbers


def test_assign_residue_numbers():
    structure = MMTFParser.get_structure_from_url("7DHI")
    nums = assign_generic_numbers(structure, "R")

    def full_id(resi):
        return (structure.id, 0, "R", (" ", resi, " "))

    # https://github.com/protwis/gpcrdb_data/blob/master/residue_data/reference_positions/adrb2_human.yaml
    # 1x50: 51
    # 12x50: 64
    # 2x50: 79
    # 23x50: 99
    # 3x50: 131
    # 34x50: 138
    # 4x50: 158
    # 45x50: 191
    # 5x50: 211
    # 6x50: 288
    # 7x50: 323
    # 8x50: 332
    assert nums[full_id(51)] == 1.50
    assert nums[full_id(64)] == 12.50
    assert nums[full_id(79)] == 2.50
    assert nums[full_id(99)] == 23.50
    assert nums[full_id(131)] == 3.50
    assert nums[full_id(138)] == 34.50
    assert nums[full_id(158)] == 4.50
    assert nums[full_id(191)] == 45.50
    assert nums[full_id(211)] == 5.50
    assert nums[full_id(288)] == 6.50
    assert nums[full_id(323)] == 7.50
    assert nums[full_id(332)] == 8.50
