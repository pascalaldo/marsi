# Copyright 2016 Chr. Hansen A/S and The Novo Nordisk Foundation Center for Biosustainability, DTU.

# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at

# http://www.apache.org/licenses/LICENSE-2.0

# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

from functools import partial

import pubchempy
from bioservices.chebi import ChEBI
from bioservices.kegg import KEGG
from bioservices.uniprot import UniProt

from cachetools import cached, LRUCache
from cachetools.keys import hashkey

from marsi.chemistry.openbabel import mol_str_to_inchi

lru_cache = LRUCache(maxsize=1024)


try:
    __all__ = ['find_chebi_id', 'find_pubchem_id', "map_uniprot_from_pdb_ids", "find_best_chebi_structure"]

    chebi_client = ChEBI()
    kegg_client = KEGG()
    uniprot_client = UniProt()

    @cached(lru_cache, key=partial(hashkey, 'uniprot'))
    def map_uniprot_from_pdb_ids(pdb_ids):
        return uniprot_client.mapping(fr="PDB_ID", to="ACC", query=pdb_ids)

    @cached(lru_cache, key=partial(hashkey, 'uniprot'))
    def inchi_from_chebi(chebi_id):
        try:
            return chebi_client.getCompleteEntity(chebi_id).inchi.strip()
        except AttributeError:
            return None

    @cached(lru_cache, key=partial(hashkey, 'uniprot'))
    def inchi_from_kegg(kegg_id):
        try:
            return mol_str_to_inchi(kegg_client.get(kegg_id, 'mol'))
        except Exception:
            return None

    def find_chebi_id(metabolite):
        """
        Queries ChEBI using InChI Key.

        Parameters
        ----------
        metabolite : Metabolite

        Returns
        -------
        list
            With tuples of (id, name)
        """
        try:
            query = chebi_client.getLiteEntity(metabolite.inchi_key, searchCategory="INCHI KEY")
            return [(e.chebiId, e.chebiAsciiName) for e in query][0]
        except IndexError:
            raise KeyError("%s not found in ChEBI" % metabolite.inchi_key)

    def find_pubchem_id(metabolite):
        """
        Queries PubChem Compound using InChI Key.

        Parameters
        ----------
        metabolite : Metabolite

        Returns
        -------
        list
            With tuples of (id, name)
        """
        try:
            return [(str(e.cid), e.iupac_name) for e in pubchempy.get_compounds(metabolite.inchi_key,
                                                                                namespace="inchikey")][0]
        except:
            raise KeyError("%s not found in PubChem Compound" % metabolite.inchi_key)

except Exception as e:
    from warnings import warn
    __all__ = ['no_services_available', 'find_best_chebi_structure']
    no_services_available = "Please check your internet connection"

    def map_uniprot_from_pdb_ids(pdb_ids):
        raise RuntimeError(no_services_available)

    def inchi_from_chebi(chebi_id):
        raise RuntimeError(no_services_available)

    def inchi_from_kegg(kegg_id):
        raise RuntimeError(no_services_available)

    def find_chebi_id(metabolite):
        raise RuntimeError(no_services_available)

    def find_pubchem_id(metabolite):
        raise RuntimeError(no_services_available)

    warn(no_services_available + " because of " + str(e))


def find_best_chebi_structure(entity):
    structure = entity.ChemicalStructures[0]
    for i, struct in enumerate(entity.ChemicalStructures):
        if struct.dimension == "3D":
            if structure.dimension == "3D" and structure.defaultStructure:
                return structure.structure
            else:
                structure = struct
        else:
            if struct.defaultStructure and not structure.dimension == "3D":
                structure = struct

    return structure.structure
