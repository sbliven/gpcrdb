# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:percent
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.10.2
#   kernelspec:
#     display_name: Python [conda env:gpcrdb]
#     language: python
#     name: conda-env-gpcrdb-py
# ---

# %% [markdown]
# # gpcrdb basics

# %%
import gpcrdb

# %%
# BioPython & mmtf-python are optional
from Bio.PDB.mmtf import MMTFParser, MMTFIO
import gpcrdb.structure as gstruc

# %% [markdown]
# ## Single-protein information

# %%
gpcrdb.fetch_protein("opsd_bovin")

# %% [markdown]
# ## Generic numbers
#
# Fetching gpcrdb entries directly:

# %%
gpcrdb.fetch_generic_numbers("opsd_bovin")[200:210]

# %% [markdown]
# Assigning numbers to your own structures:

# %%
structure = MMTFParser.get_structure_from_url("7DHI")


# %%
gstruc.assign_generic_numbers(structure, "R")

# %% [markdown]
# ## Get all structures in gpcrdb
#
# Eventually this will have a nice python wrapper too. For now it has to be done with the get method

# %%
# All families
gpcrdb.get(f"{gpcrdb.BASE_URL}/proteinfamily")[:3]

# %%
# Navigate family tree
slug = "001_001"
gpcrdb.get(f"{gpcrdb.BASE_URL}/proteinfamily/children/{slug}")[:3]

# %%
gpcrdb.get(f"{gpcrdb.BASE_URL}/proteinfamily/descendants/{slug}")[:3]

# %%
[(s['pdb_code'],s['state']) for s in gpcrdb.fetch_structures("opsd_bovin", representative=False)]

# %%
gpcrdb.get(f"{gpcrdb.BASE_URL}/proteinfamily/proteins/{slug}")[:10]
