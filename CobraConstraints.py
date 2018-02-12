import cobra
import json


def main():
    model_updates = openJson('JSONs/model_updates.json')
    model = cobra.io.read_sbml_model('iCHOv1.xml')

    reactions = []
    for reaction in model.reactions:
        reactions.append(reaction.id)

    coefficients_forward = {}
    coefficients_backward = {}
    for reaction in model_updates:
        coefficients_forward[reaction] = model_updates[reaction]['forward']
        coefficients_backward[reaction] = model_updates[reaction]['backward']

    flux_constraint(model, reactions, coefficients_forward,
                    coefficients_backward, bound=0.055)


def openJson(path):
    '''Shortened call to open JSON files.'''
    with open(path) as r:
        return json.load(r)


def flux_constraint(cobra_model, coefficients_forward, coefficients_reverse,
                    bound=1):
    """
    Adds a linear constrain of the form
        0 <= cf[1] vf[1] + cb[1] vb[1] + ... <= 1
    where vf[i], vb[i] are the forward and reverse fluxes of reaction i, and
    cf = coefficients_forward, cb = coefficients_backward are dictionaries

    Author: Jorge Fernandez-de-Cossio-Diaz
    """
    coefficients = dict()
    for (bigg_id, cf) in coefficients_forward.items():
        rxn = cobra_model.reactions.get_by_id(bigg_id)
        coefficients[rxn.forward_variable] = cf
    for (bigg_id, cr) in coefficients_reverse.items():
        rxn = cobra_model.reactions.get_by_id(bigg_id)
        coefficients[rxn.reverse_variable] = cr
    constraint = cobra_model.problem.Constraint(0, lb=0, ub=bound)
    cobra_model.add_cons_vars(constraint)
    cobra_model.solver.update()
    constraint.set_linear_coefficients(coefficients=coefficients)


if __name__ == '__main__':
    main()
