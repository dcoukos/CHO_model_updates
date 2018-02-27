import os
import sys
import json
import copy
import statistics
import cobra
from enum import Enum, auto  # Enum breaks spyder autocomplete
'''
    This file contains functions pertaining to the treatment of the data that
    is returned from BRENDA, including functions that find the KEGG Ids of that
    data, as well as functions that organize that data and simplify it in order
    to keep only the necessary information for improving the model.
'''
# pylint: disable=too-few-public-methods
# pylint: disable=too-many-branches, too-many-return-statements
# pylint: disable=W0603,R1704,R0902,C0103


def importWeights(model_update):
    '''Adds molecular weights to model enzymes'''
    mol_info = openJson('JSONs/molecular_weights.json')
    for bigg_id in mol_info:
        if 'molecular_weights' in mol_info[bigg_id]:
            if mol_info[bigg_id]['molecular_weights'] != []:
                max_weight = max(mol_info[bigg_id]['molecular_weights'])
                model_update[bigg_id].molecular_weight = max_weight
                # TODO: should this be the min weight?


def exportData(model_update):
    simple_data = {}
    for bigg_id in model_update:
        simple_data[bigg_id] = {}
        simple_data[bigg_id]['forward'] = model_update[bigg_id].f_spec_act
        simple_data[bigg_id]['backward'] = model_update[bigg_id].b_spec_act
    write('JSONs/model_updates.json', simple_data)


def fillData(model_update):
    fillEmptyValues(model_update)
    for bigg_id in model_update:
        if model_update[bigg_id].f_spec_act is None:
            model_update[bigg_id].f_spec_act = \
                float(model_update[bigg_id].forward_turnover) / \
                float(model_update[bigg_id].molecular_weight)
        if model_update[bigg_id].b_spec_act is None:
            model_update[bigg_id].b_spec_act = \
                float(model_update[bigg_id].backward_turnover) / \
                float(model_update[bigg_id].molecular_weight)

# CHANGED: I am now calculating the median molecular weights and turnovers
    # separately, to have access to more data.


def medianTurnovers(model_update):
    turnovers = []
    for bigg_id in model_update:
        if model_update[bigg_id].forward_turnover is not None:
            turnovers.append(model_update[bigg_id].forward_turnover)
        if model_update[bigg_id].backward_turnover is not None:
            turnovers.append(model_update[bigg_id].backward_turnover)
    return statistics.median(turnovers)


def medianWeights(model_update):
    importWeights(model_update)
    weights = []
    for bigg_id in model_update:
        if model_update[bigg_id].molecular_weight is not None:
            weights.append(model_update[bigg_id].molecular_weight)
    return statistics.median(weights)


def fillEmptyValues(model_update):
    for bigg_id in model_update:
        if model_update[bigg_id].forward_turnover == {}:
            model_update[bigg_id].forward_turnover is None
        if model_update[bigg_id].backward_turnover == {}:
            model_update[bigg_id].backward_turnover is None
    median_turnover = medianTurnovers(model_update)
    median_weight = medianWeights(model_update)
    for bigg_id in model_update:
        if model_update[bigg_id].forward_turnover is None:
            model_update[bigg_id].forward_turnover = median_turnover
        if model_update[bigg_id].backward_turnover is None:
            model_update[bigg_id].backward_turnover = median_turnover
        if model_update[bigg_id].molecular_weight is None:
            model_update[bigg_id].molecular_weight = median_weight


def initializeModelUpdate(model):
    '''Creates deep copy of BiGG Model representation, to which updated
        information will be added. '''

    return copy.deepcopy(model)


def treatTurnoverData(path_to_brenda_output, path_to_keggs,
                      path_to_iCHO_keggs):
    '''Filters and eliminates unnecessary data, choosing best turnover data.

    This function is meant specifically for data of the DataType.turnover Enum.
    All of the selection processes (filtering by organism, direction,
    wild-type, and magnitude) occur in functions called by
    treatTurnoverData().
    '''
    # TODO: should the model be opened as a dict or as a cobra model?
    model = openModelAsDict('iCHOv1.xml')

    potential_updates = {}

    treated_output = openJson(path_to_brenda_output)
    # brenda_keggs = correctJson(path_to_keggs)
    brenda_keggs = openJson('JSONs/brenda_keggs.json')
    loadKeggsIntoModel(model, path_to_iCHO_keggs)
    matchById(model, potential_updates, brenda_keggs, treated_output,
              DataType.turnover)
    updates = selectBestData(potential_updates)
    model_update = applyBestData(model, updates, DataType.turnover)
    fillData(model_update)
    exportData(model_update)


def loadKeggsIntoModel(model, path_to_keggs):
    '''Which model is this supposed to use again?'''
    kegg_dict = openJson(path_to_keggs)
    corrected_dict = {}
    for enzyme in kegg_dict:
        corrected_dict[enzyme] = {}
        corrected_dict[enzyme]['reactants'] = {}
        corrected_dict[enzyme]['products'] = {}
        for name in kegg_dict[enzyme]['reactants']:
            if kegg_dict[enzyme]['reactants'][name] is not None:
                corrected_dict[enzyme]['reactants'].update(
                    {code: name for code in kegg_dict
                     [enzyme]['reactants'][name]})
        for name in kegg_dict[enzyme]['products']:
            if kegg_dict[enzyme]['products'][name] is not None:
                corrected_dict[enzyme]['products'].update(
                    {code: name for code in kegg_dict
                     [enzyme]['products'][name]})
    for enzyme in model:
        if model[enzyme].bigg in corrected_dict:
            model[enzyme].with_kegg = corrected_dict[enzyme]


def matchById(model, potential_updates, brenda_keggs, treated_brenda_output,
              data_type):
    '''Tries to match metabolites by KEGG Id.

    Args:
        potential_updates: empty dict which will contain Enzyme entries. This
            dict is where the incoming data will be processed and selected. Is
            updated by this function.
        brenda_keggs: dict of KEGG codes (keys) and corresponding names
            (values)
        treated_brenda_output: This is a dict containing data from BRENDA.
            Keys are reaction Ids. First order of keys are BiGG reaction Ids,
            and the second order keys are the names of metabolites in BRENDA.
            Thereafter there is a list of outputs, obtained programatically
            from BRENDA.
        data_type: DataType which tells matchById how to add new data to the
            potential_updates dict.
        model: variable containing the representation of the model.
            Necessary as an argument for testing purposes.

    Return:
        unmatched: dict of BRENDA metabolites that could not be matched by name
    '''
    unmatched = {}

    # TODO: Things to check:

    #   3. How does the program handle empty data so as to not add to handle
    #       empty data correctly and not throw an error. (Remember not throwing
    #       an error does not mean that the situation is handled correctly.)
    for bigg_ID in model:
        unmatched[bigg_ID] = []
        potential_updates[bigg_ID] = Enzyme(bigg_ID)
        if bigg_ID in treated_brenda_output and treated_brenda_output[bigg_ID]\
                != {} and bigg_ID in brenda_keggs:
            for kegg in brenda_keggs[bigg_ID]:
                if kegg in model[bigg_ID].with_kegg['reactants'] and\
                        treated_brenda_output[bigg_ID][brenda_keggs[
                        bigg_ID][kegg]] != []:
                    name = brenda_keggs[bigg_ID][kegg]
                    bigg_name = model[bigg_ID].with_kegg['reactants'][kegg]
                    # TODO: what if the name in the xml file does not match the
                    #       name in the cobra model?
                    potential_updates[bigg_ID].forward[bigg_name] = []
                    data = getData(treated_brenda_output[bigg_ID],
                                   name, data_type)
                    for index, entry in enumerate(treated_brenda_output[
                            bigg_ID][brenda_keggs[bigg_ID][kegg]]):
                        potential_updates[bigg_ID].forward[
                            bigg_name].append(
                            MetaboliteCandidate(
                                    brenda_keggs[bigg_ID][kegg],
                                    kegg=kegg,
                                    **data[index]))
                elif kegg in model[bigg_ID].with_kegg['products'] and\
                        treated_brenda_output[bigg_ID][brenda_keggs[
                        bigg_ID][kegg]] != []:
                    name = brenda_keggs[bigg_ID][kegg]
                    bigg_name = model[bigg_ID].with_kegg['products'][kegg]
                    potential_updates[bigg_ID].backward[bigg_name] = []
                    data = getData(treated_brenda_output[bigg_ID],
                                   name, data_type)
                    for index, entry in enumerate(treated_brenda_output[
                            bigg_ID][brenda_keggs[bigg_ID][kegg]]):
                        potential_updates[bigg_ID].backward[
                            bigg_name].append(
                            MetaboliteCandidate(
                                    brenda_keggs[bigg_ID][kegg],
                                    kegg=kegg,
                                    **data[index]))
        # CHANGED: deleted returning unmatched. No longer interested because
            # we don't match by name.


def selectBestData(model_updates):
    '''Selects best data to add to model based on selection criteria

    This function selects the best data based on Organism, wild-type, and the
    magnitude of the data.

    Args:
        model_updates: contains new data for the model, is filtered by this
            function. The final candidates are returned through this arg.
        data_type: instance of DataType Enum, describing where new data is to
            be placed.

    '''
    # FIXME: eliminate entries with turnover = -999
    filtered_by_organism = selectBestOrganismEntries(model_updates)

    filtered_by_wild_type = selectWildTypeEntries(filtered_by_organism)

    filtered_by_magnitude = selectByTurnover(filtered_by_wild_type)

    return filtered_by_magnitude
    # FIXME: Somewhere below, the elminated data is being recopied.
    # -----Choose data according to magnitude preference.-----------------
    # if data_type is DataType.turnover:


def selectByTurnover(filtered_by_wild_type):
    filtered_by_magnitude = {}
    for reaction in filtered_by_wild_type:
        filtered_by_magnitude[reaction] = Enzyme(reaction)
        if filtered_by_wild_type[reaction].forward != {}:
            data = chooseHighestTurnover(
                filtered_by_wild_type[reaction].forward)
            filtered_by_magnitude[reaction].forward = data
            # If it prints, then it means that there are some Returns
            # with multiple output values. (i.e. dict with two keys. )
            try:
                print(len(data))
            except TypeError:
                pass
        if filtered_by_magnitude[reaction].forward is None:
            filtered_by_magnitude[reaction].forward = {}
        if filtered_by_magnitude[reaction].backward != {}:
            data = chooseHighestTurnover(
                filtered_by_wild_type[reaction].backward)
            try:
                print(len(data))
            except TypeError:
                pass
            filtered_by_magnitude[reaction].backward = data
        if filtered_by_magnitude[reaction].backward is None:
            filtered_by_magnitude[reaction].backward = {}
    return filtered_by_magnitude


def selectBestOrganismEntries(model_updates):
    filtered_by_organism = {}
    for reaction in model_updates:
        filtered_by_organism[reaction] = model_updates[
                                            reaction].copyOnlyEnzyme()
        filtered_by_organism[reaction].forward = {}
        filtered_by_organism[reaction].backward = {}
        for metabolite_name in model_updates[reaction].forward:
            if metabolite_name not in filtered_by_organism[reaction].forward:
                filtered_by_organism[reaction].forward[metabolite_name] = []
            (closest_organism, indices) = findClosestOrganism(model_updates[
                reaction].forward[metabolite_name])
            for organism in Organism:
                if organism == closest_organism:
                    for index in indices[organism]:
                        filtered_by_organism[reaction].forward[
                                metabolite_name].append(model_updates[
                                    reaction].forward[metabolite_name][index])

        for metabolite_name in model_updates[reaction].backward:
            if metabolite_name not in filtered_by_organism[reaction].backward:
                filtered_by_organism[reaction].backward[metabolite_name] = []
            (closest_organism, indices) = findClosestOrganism(model_updates[
                reaction].backward[metabolite_name])
            for organism in Organism:
                if organism == closest_organism:
                    for index in indices[organism]:
                        filtered_by_organism[reaction].backward[
                                metabolite_name].append(model_updates[
                                    reaction].backward[metabolite_name][index])
    return filtered_by_organism


def selectWildTypeEntries(filtered_by_organism):
    filtered_by_wild_type = {}
    for reaction in filtered_by_organism:
        filtered_by_wild_type[reaction] = filtered_by_organism[
                reaction].copyOnlyEnzyme()
        for metabolite_name in filtered_by_organism[reaction].forward:
            wild_type = []
            for entry in filtered_by_organism[reaction].forward[
                    metabolite_name]:
                if entry.wild_type is True:
                    wild_type.append(entry)
            if wild_type != []:
                filtered_by_wild_type[reaction].forward[
                    metabolite_name] = wild_type
        for metabolite_name in filtered_by_organism[reaction].backward:
            wild_type = []
            for entry in filtered_by_organism[reaction].backward[
                    metabolite_name]:
                if entry.wild_type is True:
                    wild_type.append(entry)
            if wild_type != []:
                filtered_by_wild_type[reaction].backward[
                    metabolite_name] = wild_type
    return filtered_by_wild_type


def chooseHighestTurnover(metabolite_entries):
    '''Selects entry with highest turnover

    Args:
        metabolite_entries: dict of selected entries for a reaction direction.

    Returns:
        entry with the highest turnover
    '''
    highest_index = 0
    highest_turnover = 0
    highest_metabolite = ''
    for metabolite in metabolite_entries:
        for index, entry in enumerate(metabolite_entries[metabolite]):
            if float(entry.turnover) > highest_turnover:
                highest_index = index
                highest_turnover = float(entry.turnover)
                highest_metabolite = metabolite

    if highest_metabolite == '':
        return None

    best_metabolite = metabolite_entries[highest_metabolite][highest_index]
    assert isinstance(best_metabolite, MetaboliteCandidate)
    return best_metabolite


def findClosestOrganism(metabolite_entries):
    '''finds which organism closest to hamster is in list of MetaboliteCandidate

    Args:
        metabolite_entries: list of MetaboliteCandidate

    Returns:
        Organism: enum of closest organism closest to hamster
        dict: dictionary of Organism and indices matching to the organism
    '''
    has_hamster = False
    has_mouse = False
    has_rat = False
    has_human = False
    has_fly = False
    has_yeast = False
    has_coli = False
    indices = {
        Organism.hamster: [],
        Organism.mouse: [],
        Organism.rat: [],
        Organism.human: [],
        Organism.fly: [],
        Organism.yeast: [],
        Organism.coli: []}
    for index, metabolite in enumerate(metabolite_entries):
        if metabolite.organism == 'Cricetulus griseus':
            indices[Organism.hamster].append(index)
            has_hamster = True
        elif metabolite.organism == 'Mus musculus':
            has_mouse = True
            indices[Organism.mouse].append(index)
        elif metabolite.organism == 'Rattus norvegicus':
            has_rat = True
            indices[Organism.rat].append(index)
        elif metabolite.organism == 'Homo sapiens':
            has_human = True
            indices[Organism.human].append(index)
        elif metabolite.organism == 'Drosophila melanogaster':
            has_fly = True
            indices[Organism.fly].append(index)
        elif metabolite.organism == 'Saccharomyces cerevisiae':
            has_yeast = True
            indices[Organism.yeast].append(index)
        elif metabolite.organism == 'Escherichia coli':
            has_coli = True
            indices[Organism.coli].append(index)

    if has_hamster:
        return (Organism.hamster, indices)
    elif has_mouse:
        return (Organism.mouse, indices)
    elif has_rat:
        return (Organism.rat, indices)
    elif has_human:
        return (Organism.human, indices)
    elif has_fly:
        return (Organism.fly, indices)
    elif has_yeast:
        return (Organism.yeast, indices)
    elif has_coli:
        return (Organism.coli, indices)

    return None, None


def correctJson(path):
    ''' Ensures that KEGG codes are the keys and metabolite names are the values.

    When running brendaToKegg() in Jupyter Notebook I noticed a funny error:
    some, but not all, of the key:value entries were inverted. This appeared to
    occur 'at random'. This function ensures that the KEGG code is the value.

    Args:
        path: filepath to KEGG ids and BRENDA metabolite names.

    Return:
        corrected_json: JSON object with KEGG codes as the dict keys.
    '''

    brenda_keggs = openJson(path)
    corrected_json = {}

    for reaction in brenda_keggs:
        corrected_json[reaction] = {}
        for supposed_code, name in brenda_keggs[reaction].items():
            if isNumber(supposed_code[2:]):
                corrected_json[reaction][supposed_code] = name
            elif isNumber(name[2:]):
                corrected_json[reaction][brenda_keggs[reaction]
                                         [supposed_code]] = supposed_code
            else:
                raise BadDataError

    return corrected_json


def openJson(path):
    '''Shortened call to open JSON files.'''
    if os.stat(path) == 0:
        raise NoDataError
    with open(path) as r:
        return json.load(r)


def write(path, data):
    '''Shortened call to JSON dumps with indent = 4'''
    with open(path, 'w') as wr:
        json.dump(data, wr, indent=4)


def writeEnzymes(path, data):
    '''Writes json representation of a dict of Enzymes'''
    data_repr = {}
    for entry in data:
        data_repr[entry] = data[entry].getDict()
    write(path, data_repr)


def isNumber(s):
    try:
        float(s)
        return True
    except ValueError:
        pass

    return False


def getData(reaction, metabolite, data_type):
    '''Chooses data to give to matchById.

    This function is intended to improve code reusability so that the parent
    functions can treat the data without knowing what type of data it is.

    Args:
        reaction: treated brenda output for one EC number. This is a dict
            where the keys are metabolites, and the values are lists of entries
            pertaining to that metabolite name.
        metabolite: metabolite name
        data_type: instance of DataType Enum. Tells the function which data is
            pertinent.

    Returns:
        metabolite_data: a list containing data which is used to construct
            MetaboliteCandidate instances.
    '''

    metabolite_data = []
    for entry in reaction[metabolite]:
        metabolite_entry = {}
        if data_type is DataType.turnover:
            metabolite_entry['turnover'] = entry['turnoverNumber']
        elif data_type is DataType.specific_activity:
            metabolite_entry['specific activity'] = entry['specific activity']
        else:
            metabolite_entry['molecular weight'] = entry['molecular weight']
        try:
            metabolite_entry['wild-type'] = entry['wild-type']
        except KeyError:
            pass
        metabolite_entry['organism'] = entry['organism']
        metabolite_data.append(metabolite_entry)
    return metabolite_data


def applyBestData(model, updates, data_type):
    '''Adds updated data to global MODEL_UPDATES variable
        The best data is already pared down at this point.

    Args:
        updates: updated model information. dict of Enzymes.
        data_type: type of data being added to model.
    '''
    # TODO: next line suspect.
    '''model_update = initializeModelUpdate(model)
    if model_update is None or model_update == {}:
        raise NoDataError'''
    model_update = {}
    for reaction in updates:
        model_update[reaction] = updates[reaction].copyOnlyEnzyme()
        if updates[reaction].forward != {}:
            print(reaction)
            model_update[reaction].forward = \
                updates[reaction].forward.returnAttributes()
        if updates[reaction].backward != {}:
            model_update[reaction].backward = \
                updates[reaction].backward.returnAttributes()
        if data_type is DataType.turnover:
            model_update[reaction].applyHighestTurnover()
    return model_update


def getEcNumber(reaction):
    '''Gets EC number to make new Enzyme object. "Flat is better than nested"
    '''
    for metabolite in reaction:
        for entry in reaction[metabolite]:
            try:
                return entry['ecNumber']
            except KeyError:
                continue


def cleanMetaboliteNames(brenda_metabolites):
    '''
        TODO: fill out function
    '''
    pass


def openModelAsDict(path):
    '''Creates and stores a representation of the BiGG Model.

        This creates and pickles a representation of the BiGG Model which
        should be a little bit easier to work with.
    '''
    iCHO_model = cobra.io.read_sbml_model(path)
    model = {}
    for reaction in iCHO_model.reactions:
        model[reaction.id] = Enzyme(reaction.id)

        for reactant in reaction.reactants:
            model[reaction.id].forward[reactant.name] = Metabolite(
                reactant.name, bigg=reactant.id)
        for reactant in reaction.products:
            model[reaction.id].backward[reactant.name] = Metabolite(
                reactant.name, bigg=reactant.id)
    return model


class DataType(Enum):
    turnover = auto()
    specific_activity = auto()
    molecular_weight = auto()


class Organism(Enum):
    hamster = auto()
    mouse = auto()
    rat = auto()
    human = auto()
    fly = auto()
    yeast = auto()
    coli = auto()


class Enzyme():
    '''Class containing reactions and metabolites pertaining to an enzyme.

    Attributes:
        ID: the BiGG idenitifier of the enzyme
        EC: Enzyme commision number
        metabolites: a dict of metabolites
        forward: a dict of metabolites participating in the forward reaction
        backward: a dict of metabolites participating in the backward reaction
        with_kegg: a dict of metabolites with a KEGG identifier, where the KEGG
            is the key, and the metabolite is the valye.
    '''
    def __init__(self, bigg, ec=None):
        self.bigg = bigg
        if ec is not None:
            numerals = ec.split('.')
            for number in numerals[:-1]:
                try:
                    int(number)
                except ValueError:
                    raise NotNumericError
            self.EC = ec
        else:
            self.EC = None
        self.forward = {}
        self.backward = {}
        self.with_kegg = {}
        self.forward_turnover = None
        self.backward_turnover = None
        self.has_f_wild_type = None
        self.has_b_wild_type = None
        self.molecular_weight = None
        self.f_spec_act = None
        self.b_spec_act = None

    def copyOnlyEnzyme(self):
        if self.EC:
            enz = Enzyme(self.bigg, self.EC)
        else:
            enz = Enzyme(self.bigg)
        return enz

    def applyHighestTurnover(self):
        '''Chooses best value from forward metabolites.

            This function sets the value of forward_turnover and
            backward_turnover.

        Note: The data must already be narrowed down to one entry per
            metabolite. Otherwise, this function will break
        '''
        if self.forward is not None and self.forward != {}:
            self.forward_turnover = self.forward['turnover']
        if self.backward is not None and self.backward != {}:
            self.backward_turnover = self.backward['turnover']

    def getDict(self):
        '''
            Returns a represenation of this enzyme as a dict. This allows
            testing for equality, as well as writing to a JSON format.
        '''
        return_dict = {}
        try:
            for metabolite in self.forward:
                return_dict[metabolite] = []
                for entry in self.forward[metabolite]:
                    met_representation = {}
                    met_representation['name'] = entry.name
                    met_representation['bigg'] = entry.bigg
                    met_representation['kegg'] = entry.kegg
                    met_representation['turnover'] = entry.turnover
                    met_representation['specific activity'] = \
                        entry.specific_activity
                    met_representation['organism'] = entry.organism
                    met_representation['wild-type'] = entry.wild_type
        except TypeError:
            pass
        try:
            for metabolite in self.backward:
                return_dict[metabolite] = []
                for entry in self.backward[metabolite]:
                    met_representation = {}
                    met_representation['name'] = entry.name
                    met_representation['bigg'] = entry.bigg
                    met_representation['kegg'] = entry.kegg
                    met_representation['turnover'] = entry.turnover
                    met_representation['specific activity'] = \
                        entry.specific_activity
                    met_representation['organism'] = entry.organism
                    met_representation['wild-type'] = entry.wild_type
        except TypeError:
            pass
        return return_dict

    def getSimpleDict(self):
        '''Function only meant for TESTING purposes, to be compared against
        JSON files.'''
        return_dict = {}
        try:
            return_dict[self.forward.name] = {
                'name': self.forward.name,
                'bigg': self.forward.bigg,
                'kegg': self.forward.kegg,
                'turnover': self.forward.turnover,
                'specific activity': self.forward.specific_activity,
                'organism': self.forward.organism
                }
        except KeyError:
            pass
        try:
            return_dict[self.backward.name] = {
                'name': self.backward.name,
                'bigg': self.backward.bigg,
                'kegg': self.backward.kegg,
                'turnover': self.backward.turnover,
                'specific activity': self.backward.specific_activity,
                'organism': self.backward.organism
                }
        except KeyError:
            pass
        return return_dict


class Metabolite():
    '''Class containing all information of interest pertaining to a metabolite.


    Attributes:
        name: the name of the metabolite
        bigg: BiGG identifier
        kegg: kegg identifier
        turnover: turnover number
        specific_activity: specific activity (turnover number divided by the
              molecular weight)
        molecular_weight: molecular weight
    '''

    def __init__(self, name, bigg=None, kegg=None, **kwargs):
        # CHANGED: molecular_weight is now an attribute of the enzyme, not of
        #           metabolite.
        '''TODO: Some info on this constructor'''
        self.turnover = None
        self.specific_activity = None
        self.name = name
        if bigg is None:
            self.bigg = None
        else:
            self.bigg = bigg
        if kegg is None:
            self.kegg = None
        else:
            self.kegg = kegg
        if kwargs is not None:
            for identifier, data in kwargs.items():
                if identifier == 'turnover':
                    self.turnover = data
                if identifier == 'specific activity':
                    self.specific_activity = data

    def initFromDict(self, **kwargs):

        if 'bigg' not in kwargs:
            kwargs['bigg'] = None
        if 'kegg' not in kwargs:
            kwargs['kegg'] = None
        self.__init__(kwargs['name'], kwargs['bigg'], kwargs['kegg'], kwargs)


class MetaboliteCandidate(Metabolite):
    '''Metabolite class intended to be used for data selection.

    On top of the actual values, selection of the potential metabolites data
        is based on the organism that data is provenant from, and whether the
        molecule tested was wild-type or not
    '''

    def __init__(self, name, bigg=None, kegg=None, **kwargs):
        self.organism = None
        self.wild_type = None
        super().__init__(name, bigg, kegg, **kwargs)
        if kwargs is not None:
            for identifier, value in kwargs.items():
                if identifier == 'organism':
                    self.organism = value
                if identifier == 'wild-type':
                    self.wild_type = value

    def returnAttributes(self):
        '''Only intended to be used by the best selection - Does not
            return organism or wild_type'''
        return {
            'name': self.name,
            'bigg': self.bigg,
            'kegg': self.kegg,
            'turnover': self.turnover,
            'specific activity': self.specific_activity,
            }


class DataMissingError(Exception):
    pass


class DataNotRefinedError(AttributeError):
    def __init__(self):
        print("Data not selected correctly. Expected a metabolite, found a"
              "list (Probably)")
        super().__init__(self)


class BadDataError(ValueError):
    '''The function expected a different data structure.'''
    def __init__(self):
        print("Expected a different data structure to the one found.")
        super().__init__(self)


class NoDataError(BadDataError):
    '''The function expected data, but found none.'''
    def __init__(self):
        print('Expected data, but found none.')
        super().__init__(self)


class NotNumericError(ValueError):
    '''Not numeric.'''


if __name__ == '__main__':
    if len(sys.argv) == 1:
        new_data = treatTurnoverData('JSONs/treated_brenda_output.json',
                                     'JSONs/brenda_keggs.json',
                                     'JSONs/iCHOv1_keggs.json')
        # write new data as dict --> add constraints to Cobra model.
