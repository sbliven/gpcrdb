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

# %%
from gpcrdb import GPCRdb

# %% [markdown]
# ## GPCRdb basics
#
# ### Generic numbers

# %%
GPCRdb.fetch_generic_numbers("opsd_bovin")[:10]

# %% [markdown]
# ## Get all structures in GPCRdb

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

# %% [markdown]
# ## Access PDBe-KB
#
# This isn't really part of GPRCdb, but can be useful for fetching additional info.

# %%
from neo4j import GraphDatabase


# %%
class PDBeKB(object):
    def __init__(self, uri, user, password):
        self._driver = GraphDatabase.driver(uri, auth=(user, password))

    def close(self):
        self._driver.close()

    def execute(self, cypher, params={}):
        with self._driver.session() as session:
            result = session.write_transaction(self._run_transaction, cypher, params)
            return result

    @staticmethod
    def _run_transaction(tx, cypher, params={}):
        result = tx.run(cypher, **params)
        return result.data()


# %%
uri = "bolt://18.189.182.229:7687"

# %%
db = PDBeKB(uri, "", "")

# %%
db.execute(
    """
MATCH (n:Entry)-[:HAS_ENTITY]->(m)
RETURN n.ID, m.ID
LIMIT 10;
"""
)

# %%

# %%
"""
// Problem 10
WITH '6cmo' AS pdbId, 'A' AS chain
// Get the UniRef90 cluster for a given chain
MATCH (p:Entry {ID: pdbId})-[:HAS_ENTITY]->(e)-[:CONTAINS_CHAIN]->(c:Chain {AUTH_ASYM_ID: chain}),
(e)-[:HAS_UNIPROT]->()-[:IS_IN_UNIREF90_CLUSTER]->(uref:UniRef_Reference)
MATCH (uref)<-[:IS_IN_UNIREF90_CLUSTER]-(u:UniProt)
WITH u
// Get residues from this cluster that interact with proteins
MATCH 
(u)<-[seg:HAS_UNIREF90_SEGMENT]-(:Entity)
-[:HAS_PDB_RESIDUE]->(r1:PDB_Residue)
-[:HAS_PISA_BOND]-(r2:PDB_Residue)
<-[:HAS_PDB_RESIDUE]-(partnerEnt:Entity)
-[:HAS_UNIPROT]->(partner:UniProt),
(partnerEnt)<-[:HAS_ENTITY]-(partnerEntry:Entry)
WHERE toInteger(r1.ID) IN RANGE(toInteger(seg.PDB_START),toInteger(seg.PDB_END))
WITH u.ACCESSION AS source, (toInteger(r1.ID) - toInteger(seg.PDB_START)) AS alignPos, partner.ACCESSION as partner, partnerEntry.ID as partnerEntry
RETURN DISTINCT alignPos, source, partner, partnerEntry
ORDER BY alignPos
"""

# %%
d = {"pdbId": "6cmo", "chain": "A"}
db.execute(
    """
MATCH (p:Entry {ID: $pdbId})-[:HAS_ENTITY]->(e)-[:CONTAINS_CHAIN]->(c:Chain {AUTH_ASYM_ID: $chain}),
(e)-[:HAS_UNIPROT]->()-[:IS_IN_UNIREF90_CLUSTER]->(uref:UniRef_Reference)
MATCH (uref)<-[:IS_IN_UNIREF90_CLUSTER]-(u:UniProt)
WITH u
// Get residues from this cluster that interact with proteins
MATCH 
(u)<-[seg:HAS_UNIREF90_SEGMENT]-(:Entity)
-[:HAS_PDB_RESIDUE]->(r1:PDB_Residue)
-[:HAS_PISA_BOND]-(r2:PDB_Residue)
<-[:HAS_PDB_RESIDUE]-(partnerEnt:Entity)
-[:HAS_UNIPROT]->(partner:UniProt),
(partnerEnt)<-[:HAS_ENTITY]-(partnerEntry:Entry)
WHERE toInteger(r1.ID) IN RANGE(toInteger(seg.PDB_START),toInteger(seg.PDB_END))
WITH u.ACCESSION AS source, (toInteger(r1.ID) - toInteger(seg.PDB_START)) AS alignPos, partner.ACCESSION as partner, partnerEntry.ID as partnerEntry
RETURN DISTINCT alignPos, source, partner, partnerEntry
ORDER BY alignPos
""",
    d,
)

# %%
db.close()

# %% [markdown]
# ### With py2neo

# %%
import py2neo
import py2neo.ogm


# %%
class PDBeKB2neo(object):
    def __init__(self, uri, **args):
        self._graph = py2neo.Graph(uri, **args)

    def close(self):
        self._graph.close()

    def get_UniRef90_siblings(self, pdbid, auth_chain):
        self.execute(
            """
MATCH
(entry:Entry {ID: $pdbId})-[:HAS_ENTITY]->
    (seed:Entity)-[:CONTAINS_CHAIN]->
    (seedChain:Chain {AUTH_ASYM_ID: $chain}),
(seed)-[:HAS_UNIPROT]->
    (seedU:UniProt)-[:IS_IN_UNIREF90_CLUSTER]->
    (uref:UniRef_Reference)
MATCH (uref)<-[:IS_IN_UNIREF90_CLUSTER]-(member:UniProt)
RETURN member.ACCESSION"""
        )

    def execute(self, cypher, params={}):
        results = self._graph.run(cypher, params)
        return results


# %%
db2 = PDBeKB2neo(uri)

# %%
[r[0] for r in db2.execute("MATCH (u:UniProt) RETURN u.ACCESSION LIMIT 5")]


# %%
class UniProt(py2neo.ogm.GraphObject):
    ACCESSION = py2neo.ogm.Property()


# %%
[a.ACCESSION for a in UniProt.match(db2._graph).limit(3)]

# %%
u = next(iter(UniProt.match(db2._graph).limit(1)))

# %%

# %%

# %%
