import os
import sys
import json
import copy
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


def initializeModelUpdate(model):
    '''Creates deep copy of BiGG Model representation, to which updated
        information will be added. '''

    return copy.deepcopy(model)


def treatTurnoverData(path_to_brenda_output, path_to_keggs):
    '''Filters and eliminates unnecessary data, choosing best turnover data.

    This function is meant specifically for data of the DataType.turnover Enum.
    All of the selection processes (filtering by organism, direction,
    wild-type, and magnitude) occur in functions called by
    treatTurnoverData().
    '''
    model = openModelAsDict('iCHOv1.xml')

    potential_updates = {}

    treated_output = openJson(path_to_brenda_output)
    brenda_keggs = correctJson(path_to_keggs)

    loadKeggsIntoModel(model, path_to_keggs)
    matchById(model, potential_updates, treated_output, brenda_keggs,
              DataType.turnover)
    print(potential_updates['10FTHF5GLUtl'].forward)
    # Is this an empty metabolite?  --> yes
    selectBestData(potential_updates, DataType.turnover)

    # TODO: potential_updates must have only one metabolite per enzyme
    # direction at this point.
    applyBestData(model, potential_updates, DataType.turnover)

    # TODO: apply appropriate return type.


def loadKeggsIntoModel(model, path_to_keggs):
    '''Which model is this supposed to use again?'''
    kegg_dict = openJson(path_to_keggs)
    for enzyme in model:
        if model[enzyme].bigg in kegg_dict:
            model[enzyme].with_kegg = kegg_dict[enzyme]


def matchById(model, potential_updates, brenda_keggs, treated_brenda_output,
              data_type):
    # TODO: Make sure that 'BIGG_MODEL' is the iCHOv1 model and not the bigg
        # global model.
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
    # Somewhere here the metabolite is not being added to the forward
    #    dictionary
    # TODO: Should have had an error. Is someone reading iCHOv1_keggs?
    # Nobody is reading that file.

    # TODO: Things to check:
    #   1. Where are the keggs from the iCHOv1 model being matched? Should be
    #        in addKeggsToModel or some such shazam.
    #   2. What happens to the entries that have null Kegg codes. How is the
    #       new data going to be added to the iCHOv1 model?
    #   3. How does the program handle empty data so as to not add to handle
    #       empty data correctly and not throw an error. (Remember not throwing
    #       an error does not mean that the situation is handled correctly.)
    for bigg_ID in model:
        unmatched[bigg_ID] = []
        potential_updates[bigg_ID] = Enzyme(bigg_ID)
        if bigg_ID in treated_brenda_output and treated_brenda_output[bigg_ID]\
                != {} and bigg_ID in brenda_keggs:
            for kegg in brenda_keggs[bigg_ID]:
                if kegg in model[bigg_ID].with_kegg:
                    if model[bigg_ID].with_kegg[kegg] in model[
                            bigg_ID].forward and treated_brenda_output[
                            bigg_ID][brenda_keggs[bigg_ID][kegg]] != []:
                        name = brenda_keggs[bigg_ID][kegg]
                        bigg_name = model[bigg_ID].with_kegg[kegg]
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
                    elif model[bigg_ID].with_kegg[kegg] in model[
                            bigg_ID].backward and treated_brenda_output[
                            bigg_ID][brenda_keggs[bigg_ID][kegg]] != []:
                        name = brenda_keggs[bigg_ID][kegg]
                        bigg_name = model[bigg_ID].with_kegg[kegg]
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
    # CHECK: could this be because the value is null in iCHOv1_keggs
    # for 10FTHF5GLUtl?
    assert '10FTHF5GLUtl' in model
    print(potential_updates['10FTHF5GLUtl'].forward)


def selectBestData(model_updates, data_type):
    '''Selects best data to add to model based on selection criteria

    This function selects the best data based on Organism, wild-type, and the
    magnitude of the data.

    Args:
        model_updates: contains new data for the model, is filtered by this
            function. The final candidates are returned through this arg.
        data_type: instance of DataType Enum, describing where new data is to
            be placed.

    '''
    # TODO: but model updates is empty!!!!!
    # first the data must be filtered by organism.
    filtered_by_organism = {}
    for reaction in model_updates:
        # checked
        filtered_by_organism[reaction] = model_updates[
                reaction].copyOnlyEnzyme()
        # try:
        for metabolite_name in model_updates[reaction].forward:
            print(metabolite_name)
            filtered_by_organism[reaction].forward[metabolite_name] = []
            (closest_organism, indices) = findClosestOrganism(model_updates[
                reaction].forward[metabolite_name])
            best_found = False
            for organism in Organism:
                if best_found is False and closest_organism is not None:
                    for index in indices[organism]:
                        filtered_by_organism[reaction].forward[
                                metabolite_name].append(model_updates[
                                    reaction].forward[metabolite_name][index])
                        print(model_updates[
                            reaction].forward[metabolite_name][index].turnover)
                        if organism is closest_organism:
                            best_found = True
        # except TypeError:
        #     pass
        # SAME as above but for reverse reaction
        # TODO: why are there exception blocks here?
        # try:
        for metabolite_name in model_updates[reaction].backward:
            print(metabolite_name)
            filtered_by_organism[reaction].backward[metabolite_name] = []
            closest_organism, indices = findClosestOrganism(model_updates[
                reaction].backward[metabolite_name])
            best_found = False
            for organism in Organism:
                if best_found is False and closest_organism is not None:
                    for index in indices[organism]:
                        filtered_by_organism[reaction].backward[
                            metabolite_name].append(model_updates[
                                reaction].backward[metabolite_name][index])
                        if organism is closest_organism:
                            best_found = True
        # except TypeError:
            # pass
    print(filtered_by_organism['10FTHF5GLUtl'].forward)  # TODO: Porque Pablo?
    filtered_by_wild_type = {}
    for reaction in filtered_by_organism:
        filtered_by_wild_type[reaction] = filtered_by_organism[
                reaction].copyOnlyEnzyme()
    # try:
        for metabolite_name in filtered_by_organism[reaction].forward:
            wild_type = []
            for entry in filtered_by_organism[reaction].forward[
                    metabolite_name]:
                if entry.wild_type is not None:
                    wild_type.append(entry)
            if wild_type != []:
                filtered_by_wild_type[reaction].forward[
                    metabolite_name] = wild_type
            else:
                filtered_by_wild_type[reaction].forward[
                    metabolite_name] = copy.copy(filtered_by_organism
                                                 [reaction].forward[
                                                     metabolite_name])
    # except TypeError:
        # pass
        # SAME as above but for reverse reactions
    # TODO: why are there exception blocks here?
        # try:
        for metabolite_name in filtered_by_organism[reaction].backward:
            wild_type = []
            for entry in filtered_by_organism[reaction].backward[
                    metabolite_name]:
                if entry.wild_type is not None:
                    wild_type.append(entry)
            if wild_type != []:
                filtered_by_wild_type[reaction].backward[
                    metabolite_name] = wild_type
            else:
                filtered_by_wild_type[reaction].backward[
                    metabolite_name] = copy.copy(filtered_by_organism
                                                 [reaction].backward[
                                                     metabolite_name])
        # except TypeError:
            # pass

    # -----Choose data according to magnitude preference.-----------------
    # work on filtered_by_wild_type)
# TODO: below block should return a single metabolite per direction
# TODO: Why is there a key for the metabolite that was being returned?
    if data_type is DataType.turnover:
        print('Selecting highest turnover')
        for reaction in filtered_by_wild_type:
            if filtered_by_wild_type[reaction].forward != {}:
                data = chooseHighestTurnover(
                    filtered_by_wild_type[reaction].forward)
                model_updates[reaction].forward = data
            if filtered_by_wild_type[reaction].backward != {}:
                data = chooseHighestTurnover(
                    filtered_by_wild_type[reaction].backward)
                model_updates[reaction].backward = data

    print(model_updates['10FTHF5GLUtl'].forward)


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

    # NOTE: Returns a dict.
    return metabolite_entries[highest_metabolite][highest_index]


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
    # TODO: ensure that incoming model udpate only has one metabolite per
    # enzyme direction.

    model_update = initializeModelUpdate(model)
    if model_update is None or model_update == {}:
        raise NoDataError
    for reaction in updates:
        # TODO: figure out why the forward and backward return values are {}
        print(reaction)
        print(str(updates[reaction].forward))
        print(str(updates[reaction].backward))
        model_update[reaction].forward = \
            model[reaction].forward.returnAttributes()
        model_update[reaction].backward = \
            model[reaction].backward.returnAttributes()
        if data_type is DataType.turnover:
            # WARNING: applyHighestTurnover requires that highest metabolites
            # be chosen and forward and backward be replaced with dicts from
            # those metabolites
            model_update[reaction].applyHighestTurnover()


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

    def copyOnlyEnzyme(self):
        if self.EC:
            enz = Enzyme(self.bigg, self.EC)
        else:
            enz = Enzyme(self.bigg)
        return enz

    def applyHighestTurnover(self):
        '''Chooses best value from forward metabolites.

            This function sets the value of forward_tunover and
            backward_turnover.

        Note: The data must already be narrowed down to one entry per
            metabolite. Otherwise, this function will break
        '''
        counter = 0
        try:
            counter += 1
            print('Turnovers applied: ' + str(counter))
            if self.forward is not None:
                print(len(self.forward))
                # BUG: Sometimes this returns with 2 metabolties
                print(self.forward.keys())
                self.forward_turnover = self.forward['turnover']
            if self.backward is not None:
                self.backward_turnover = self.backward['turnover']
        except KeyError:
            print(self.forward)

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
                    met_representation['molecular weight'] = \
                        entry.molecular_weight
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
                    met_representation['molecular weight'] = \
                        entry.molecular_weight
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
                'molecular weight': self.forward.molecular_weight,
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
                'molecular weight': self.backward.molecular_weight,
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
        '''TODO: Some info on this constructor'''
        self.turnover = None
        self.molecular_weight = None
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
                if identifier == 'molecular weight':
                    self.molecular_weight = data

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
            'molecular weight': self.molecular_weight
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
                                     'JSONs/brenda_keggs.json')
        # write new data as dict --> add constraints to Cobra model.
