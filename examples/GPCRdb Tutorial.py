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
# # GPCRdb basics

# %%
from gpcrdb import GPCRdb

# %% [markdown]
# ## Single-protein information

# %%
GPCRdb.fetch_protein("opsd_bovin")

# %% [markdown]
# ## Generic numbers

# %%
GPCRdb.fetch_generic_numbers("opsd_bovin")[:10]

# %% [markdown]
# ## Get all structures in GPCRdb
#
# Eventually this will have a nice python wrapper too. For now it has to be done with the get method

# %%
# All families
GPCRdb.get(f"{GPCRdb.BASE_URL}/proteinfamily")[:10]

# %%
# Navigate family tree
slug = "001_001"
GPCRdb.get(f"{GPCRdb.BASE_URL}/proteinfamily/children/{slug}")[:10]

# %%
GPCRdb.get(f"{GPCRdb.BASE_URL}/proteinfamily/descendants/{slug}")[:10]

# %%
GPCRdb.fetch_structures("opsd_bovin", representative=False)[:10]

# %%
GPCRdb.get(f"{GPCRdb.BASE_URL}/proteinfamily/proteins/{slug}")[:10]
