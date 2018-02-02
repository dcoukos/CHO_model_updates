import cobra
import json


def main():
    model_updates = openJson('JSONs/model_updates.json')
    # TODO: Dont pickle the whole model! Create a simpler dictionary.
    model = cobra.io.read_sbml_model('iCHOv1.xml')

    flux_expression = 0
    for bigg_id in model_updates:
        assert bigg_id in model.reactions
        assert 0 <= model_updates[bigg_id]['forward'] < float('inf')
        assert 0 <= model_updates[bigg_id]['backward'] < float('inf')
        flux_expression += model_updates[bigg_id]['forward'] * \
            bigg_id.backward_variable
        flux_expression += model_updates[bigg_id]['backward'] * \
            bigg_id.forward_variable
    model.problem.Constraint(flux_expression, lb=0, ub=1)


def openJson(path):
    '''Shortened call to open JSON files.'''
    with open(path) as r:
        return json.load(r)


if __name__ == '__main__':
    main()
