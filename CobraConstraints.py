import cobra
import json
from progress.bar import Bar


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
                    coefficients_backward)


def openJson(path):
    '''Shortened call to open JSON files.'''
    with open(path) as r:
        return json.load(r)


def enzymatic_expression(cobra_model, reactions, coefficients_forward,
                         coefficients_reverse):
    """
    A linear expression of the form,
        cf[1] vf[1] + cb[1] vb[1] + ...
    where vf[i], vb[i] are the forward and reverse fluxes of reaction i, and
    cf = coefficients_forward, cb = coefficients_reverse are dictionaries
    containing the coefficients of each reaction.

    Author: Jorge Fernandez-de-Cossio-Diaz
    """
    total = len(reactions)
    bar = Bar('Adding reaction coefficients to expression: ', max=total)
    expr = 0
    for bigg_id in reactions:
        rxn = cobra_model.reactions.get_by_id(bigg_id)
        assert 0 <= coefficients_forward[bigg_id] < float('inf')
        assert 0 <= coefficients_reverse[bigg_id] < float('inf')
        expr += coefficients_forward[bigg_id] * rxn.forward_variable
        expr += coefficients_reverse[bigg_id] * rxn.reverse_variable
        bar.next()
    return expr


def flux_constraint(cobra_model, reactions, coefficients_forward,
                    coefficients_reverse):
    """
    Adds a linear constrain of the form
        0 <= cf[1] vf[1] + cb[1] vb[1] + ... <= 1
    where vf[i], vb[i] are the forward and reverse fluxes of reaction i, and
    cf = coefficients_forward, cb = coefficients_backward are dictionaries.

    Author: Jorge Fernandez-de-Cossio-Diaz
    """
    expr = enzymatic_expression(cobra_model, reactions, coefficients_forward,
                                coefficients_reverse)
    cobra_model.problem.Constraint(expr, lb=0, ub=1)


if __name__ == '__main__':
    main()
