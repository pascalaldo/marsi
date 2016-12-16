# Copyright 2016 The Novo Nordisk Foundation Center for Biosustainability, DTU.

# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at

# http://www.apache.org/licenses/LICENSE-2.0

# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

from pandas import DataFrame

from marsi.design.operations import ReplaceKnockoutFinder
from IProgress import ProgressBar, Bar, Percentage


def replace_knockouts(model, knockouts, objective_function, simulation_method, simulation_kwargs, currency_metabolites):
        knockout_replace_finder = ReplaceKnockoutFinder(model,
                                                        simulation_method=simulation_method,
                                                        simulation_kwargs=simulation_kwargs,
                                                        currency_metabolites=currency_metabolites)

        result = DataFrame(columns=['reactions', 'replaced_reaction', 'metabolites', 'original_fitness', 'anti_metabolite_fitness'])
        pbar = ProgressBar(maxval=len(knockouts), widgets=["Replacing KOs: ", Bar(), Percentage()])
        replacements = [knockout_replace_finder(row.reactions, row.fitness, objective_function, None)
                        for row in pbar(knockouts.itertuples(index=False))]
        i = 0
        for replace, original in zip(replacements, knockouts.itertuples(index=False)):
            if len(replace) == 0:
                continue

            for row in replace.itertuples(index=False):
                reactions = list(original.reactions)
                reactions.remove(row.reaction)
                for metabolite in row.metabolites:
                    result.loc[i] = [reactions, row.reaction, metabolite, original.fitness, row.fitness]
                    i += 1

        return result
