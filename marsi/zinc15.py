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

import requests


BASE_URL = "http://zinc15.docking.org/"


def _format(request_format, format_options):
    if request_format is None:
        return ""
    else:
        f = ".format"
        if format_options:
            f += ":" + "+".join(format_options)

    return f


def _get(resource, identifier=None, relation=None, subsets=None, request_format="json", format_options=None):
    request_format = _format(request_format, format_options)

    url = BASE_URL + resource + request_format

    if identifier:
        url += "/%s" % identifier
    if relation:
        url += "/%s" % relation

    return requests.get(url)


def get_substance(identifier):
    response = _get("substances", identifier=identifier)
    print(response.text)


def get_substance_analogs(identifier):
    response = _get("substances", relation="ecfp4_fp-tanimoto-50=%s" % identifier)
    print(response.text)
