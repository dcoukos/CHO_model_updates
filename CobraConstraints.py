import json
import pandas
import cobra
from optlang.symbolics import Zero


def main():
    k1_model = cobra.io.read_sbml_model('iCHOv1_K1_final.xml')
    k1_updates = openJson('JSONs/k1_updates.json')
    coef_forward = {}
    coef_backward = {}
    release_bounds(k1_model)
    sol = k1_model.optimze()
    # TODO: solution used before it is defined? I have added a first FBA to
    #           estimate the enzyme mass.
    enzyme_mass = set_enzyme_mass(sol, coef_forward, coef_backward)
    flux_constraint(k1_model, coef_forward, coef_backward, enzyme_mass)
    get_coefficients(k1_updates, coef_forward, coef_backward)
    fba_and_min_enzyme(k1_model, coef_forward, coef_backward)
    k1_model.summary()


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


def get_coefficients(model_updates, coef_forward, coef_backward):
    for reaction in model_updates:
        coef_forward[reaction] = model_updates[reaction]['forward']
        coef_backward[reaction] = model_updates[reaction]['backward']


def fba_and_min_enzyme(cobra_model, coefficients_forward,
                       coefficients_reverse):
    """
    Performs FBA follows by minimization of enzyme content
    """

    with cobra_model as model:
        model.optimize()
        cobra.util.fix_objective_as_constraint(model)
        set_enzymatic_objective(model, coefficients_forward,
                                coefficients_reverse)
        sol = cobra_model.optimize()
        return sol


def release_bounds(model):
    for rxn in model.reactions:
        if rxn.upper_bound == rxn.lower_bound:
            if rxn.upper_bound < 0:
                rxn.upper_bound = 0
            if rxn.upper_bound > 0:
                rxn.lower_bound = 0


def set_enzymatic_objective(cobra_model, coefficients_forward,
                            coefficients_reverse):
    coefficients = {}
    for (bigg_id, cf) in coefficients_forward.items():
        rxn = cobra_model.reactions.get_by_id(bigg_id)
        coefficients[rxn.forward_variable] = cf
    for (bigg_id, cr) in coefficients_reverse.items():
        rxn = cobra_model.reactions.get_by_id(bigg_id)
        coefficients[rxn.reverse_variable] = cr

    cobra_model.objective = \
        cobra_model.problem.Objective(Zero, direction='min', sloppy=True,
                                      name="min_enzymatic")

    cobra_model.objective.set_linear_coefficients(coefficients=coefficients)


def constrain_uptakes(model, xi):
    # TODO: Write mathematical representation of constraint.
    '''Uptake is defined by culture conditions, but this should depend _only_
    on xi. We consider that the metabolite concentrations and such in the
    culture are at a steady state. '''
    v_glc = 0.9  # mM/min
    v_aa = 0.09  # mM/min
    # TODO: On what time interval are uptakes defined in the model.
    IMDM = pandas.read_table('IMDM.txt', comment='#', sep='\s+')
    conc = {}
    for index, row in IMDM.iterrows():
        conc[row['id']] = row['mM']
    for met in conc:
        if met == 'glc_D_e':
            if xi == 0:
                return v_glc
            model.reactions.get_by_id('EX_glc_e_').lower_bound = \
                -min(v_glc, conc[met]/xi)
        else:
            if xi == 0:
                return v_aa
            name = met.split('_')[0]
            model.reactions.get_by_id('EX_' + name + '_e_').lower_bound = \
                -min(v_aa, conc[met]/xi)


def set_enzyme_mass(solution, coef_forward, coef_backward):
    enzyme_mass = 0
    for rxn, cf in coef_forward.items():
        enzyme_mass += cf * max(solution.fluxes[rxn], 0)
    for rxn, cb in coef_backward.items():
        enzyme_mass -= cb * min(solution.fluxes[rxn], 0)
    return enzyme_mass


def exchanges_secretion(cobra_model):
    "Subset of exchanges that can secrete."
    return [rxn for rxn in cobra_model.exchanges if
            not rxn.products and rxn.upper_bound > 0 or
            not rxn.reactants and rxn.lower_bound < 0]


def exchanges_consumption(cobra_model):
    "Subset of exchanges that can consume."
    return [rxn for rxn in cobra_model.exchanges if
            not rxn.products and rxn.lower_bound < 0 or
            not rxn.reactants and rxn.upper_bound > 0]


def osmolarity(imports, exports):
    osmo = float(0)
    permeable = ['o2_e', 'co2_e', 'h2o_e']
    for mol in imports:
        if mol not in permeable:
            osmo -= imports[mol]
    for mol in exports:
        if mol not in permeable:
            osmo += exports[mol]
    return osmo


if __name__ == '__main__':
    main()
