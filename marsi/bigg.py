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
"""
BiGG database API v2
"""
import os
import requests

BASE_URL = "http://bigg.ucsd.edu/api/v2/"


class DBVersion(object):
    def __init__(self, version, api_version, last_update):
        self.version = version
        self.api_version = api_version
        self.last_update = last_update

    def __repr__(self):
        return "BiGG Database version %s (API %s). Last update: %s" % (self.version, self.api_version, self.last_update)


def database_version():
    """
    Retrieves the current version of BiGG database
    """
    response = requests.get(BASE_URL + "database_version")
    if response.ok:
        data = response.json()
        return DBVersion(data['bigg_models_version'], data['api_version'], data['last_updated'])


def download_model(model_id, file_format="json", save=True, path="."):
    """
    Download models from BiGG. You can chose to save the file or to return the JSON data.

    Arguments
    ---------

    model_id: str
        A valid id for a model in BiGG.
    file_format: str
        If you want to save the file, you can import the model in the following formats:
            1. json (JSON format)
            2. xml (SBML)
            3. xml.gz (SBML compressed)
            4. mat (MATLAB)
    save: bool
        If True, writes the model to a file with the model name (the path can be specified).
    path: str
        Specifies in which folder the model should be written if *save* is True.
    """

    if save:
        response = requests.get("http://bigg.ucsd.edu/static/models/%s.%s" % (model_id, file_format), stream=True)
        if response.ok:
            with open(os.path.join(path, "%s.%s" % (model_id, file_format)), "w") as model_file:
                for block in response.iter_content(1024):
                    model_file.write(block)

    else:
        response = requests.get(BASE_URL + "models/%s/%s" % (model_id, file_format))
        if response.ok:
            return response.json()


def model_details(model_id):
    """
    Summarize the model.

    Arguments
    ---------
    model_id: str
        A valid id for a model in BiGG.
    """
    response = requests.get(BASE_URL + "models/%s" % model_id)
    if response.ok:
        return response.json()


def list_models():
    """
    Listsall models available in BiGG.
    """
    response = requests.get(BASE_URL + "models/")
    if response.ok:
        return response.json()


def list_reactions():
    """
    List all reactions available in BiGG.
    """
    response = requests.get(BASE_URL + "universal/reactions")
    if response.ok:
        return response.json()


def list_model_reactions(model_id):
    """
    List all reactions in a model.

    Arguments
    ---------
    model_id: str
        A valid id for a model in BiGG.
    """
    response = requests.get(BASE_URL + "models/%s/reactions" % model_id)
    if response.ok:
        return response.json()


def get_reaction(reaction_id):
    """
    Retrieve a reaction from BiGG.

    Arguments
    ---------
    reaction_id: str
        A valid id for a reaction in BiGG.
    """
    response = requests.get(BASE_URL + "universal/reactions" % reaction_id)
    if response.ok:
        return response.json()


def get_model_reaction(model_id, reaction_id):
    """
    Retrieve a reaction in the context of a model from BiGG.

    Arguments
    ---------
    model_id: str
        A valid id for a model in BiGG.
    reaction_id: str
        A valid id for a reaction in BiGG.
    """
    response = requests.get(BASE_URL + "models/%s/reactions/%s" % (model_id, reaction_id))
    if response.ok:
        return response.json()


def list_metabolites():
    """
    List all metabolites in BiGG.
    """
    response = requests.get(BASE_URL + "universal/metabolites")
    if response.ok:
        return response.json()


def list_model_metabolites(model_id):
    """
    List all metabolites in a model.

    Arguments
    ---------
    model_id: str
        A valid id for a model in BiGG.
    """
    response = requests.get(BASE_URL + "models/%s/metabolites" % model_id)
    if response.ok:
        return response.json()


def get_metabolite(metabolite_id):
    """
    Retrieve a metabolite from BiGG.

    Arguments
    ---------
    metabolite_id: str
        A valid id for a reaction in BiGG.
    """
    response = requests.get(BASE_URL + "universal/metabolites/%s" % metabolite_id)
    if response.ok:
        return response.json()


def list_model_genes(model_id):
    """
    List all genes in a model.

    Arguments
    ---------
    model_id: str
        A valid id for a model in BiGG.
    """
    response = requests.get(BASE_URL + "models/%s/genes" % model_id)
    if response.ok:
        return response.json()


def get_model_gene(model_id, gene_id):
    """
    Retrieve a gene in the context of a model from BiGG.

    Arguments
    ---------
    model_id: str
        A valid id for a model in BiGG.
    gene_id: str
        A valid id for a gene in BiGG.
    """
    response = requests.get(BASE_URL + "models/%s/metabolites/%s" % (model_id, gene_id))
    if response.ok:
        return response.json()


def search(query, search_type):
    """
    Fuzzy search the BiGG database.

    Arguments
    ---------
    query: str
        Whatever you are searching for (ids, names, etc...).
    search_type: str
        Search domain. One of "models", "genes", "reactions", "metabolites".

    """
    response = requests.get(BASE_URL + "search", params=dict(query=query, search_type=search_type))
    if response.ok:
        return response.json()
