import os
import sys
import json
import pandas
import numpy
import cobra
from multiprocessing import Pool
from matplotlib.pyplot import figure
from optlang.symbolics import Zero
from cobra.exceptions import Infeasible
from progress.bar import Bar

# FIXME: solver status is infeasible before any minimization
# FIXME: python memory consumption out of hand. change backend?


def main():
    k1_model = cobra.io.read_sbml_model('iCHOv1_K1_final.xml')
    k1_updates = openJson('JSONs/k1_updates.json')
    population_osmolarities(k1_model, k1_updates)


def find_max_xi(model, updates, tolerance=1):
    assert tolerance > 0
    min_xi = 0.
    max_xi = 1000.
    max_undefined = True
    print("finding an infeasible xi")
    while max_undefined:
        # FIXME: Why do I get an error why I comment the next two lines, and
        # why does the speed not change?
        print('.')
        try:
            run_fba(model, updates, max_xi)
            max_xi *= 2
        except Infeasible:
            max_undefined = False
    print("\nfinding max xi", end='')
    while max_xi - min_xi > tolerance:
        assert min_xi <= max_xi
        avg_xi = (max_xi + min_xi)/2
        print('min_xi=', min_xi, ', max_xi=', max_xi, ', avg_xi=', avg_xi,
              flush=True)
        try:
            run_fba(model, updates, avg_xi)
            min_xi = avg_xi
        except Infeasible:
            max_xi = avg_xi
    print('Returning ' + str(min_xi))
    return min_xi


def run_fba(model, updates, xi):
    enzyme_mass = float(sys.argv[1])
    release_bounds(model)
    constrain_uptakes(model, xi)
    coef_forward = {}
    coef_backward = {}
    model.reactions.DM_atp_c_.lower_bound = 1
    get_coefficients(updates, coef_forward, coef_backward, 3.6E6)
    flux_constraint(model, coef_forward, coef_backward, enzyme_mass)
    # m_coefs = mitochondrial_coefs(model, coef_forward, coef_backward)
    # FIXME: second call to flux_constraint freezes find_max_xi.
    # flux_constraint(model, *m_coefs, enzyme_mass)
    return fba_and_min_enzyme(model, coef_forward, coef_backward)


def mitochondrial_coefs(model, coef_forward, coef_backward):
    m_coef_forward = {}
    m_coef_backward = {}
    for reaction in model.reactions:
        if 'm' in reaction.compartments:
            m_coef_forward[reaction.id] = coef_forward[reaction.id]
            m_coef_backward[reaction.id] = coef_backward[reaction.id]
    return coef_forward, coef_backward


def subprocess(model, updates, min_xi, max_xi, slices, medium_osmolarity,
               process=4):
    if process == 1:
        bar = Bar('Calculating osmolarities: ', max=slices)
    exchanged = []
    osmolarities = []
    xis = []
    for xi in numpy.linspace(min_xi, max_xi, slices):
        # model = copy.deepcopy(orig_model)
        if process == 1:
            bar.next()
        try:
            sol = run_fba(model, updates, xi)
            osmo = osmolarity(model, sol)
            xis.append(xi)
            to_append = {}
            for mol in osmo:
                if mol not in exchanged:
                    exchanged.append(mol)
                if mol == 'total':
                    to_append[mol] = medium_osmolarity + osmo[mol]*xi
                else:
                    to_append[mol] = osmo[mol]*xi
            osmolarities.append(to_append)
        except Infeasible:
            print('xi value: ' + str(xi) + 'infeasible')
            break
    return {
        'xi': xis,
        'osmolarities': osmolarities,
        'ex': exchanged
    }


def population_osmolarities(model, updates, min_xi=0):
    # 1st infeasible value = 1202.9
    # max_xi = 1200/3
    enzyme_mass = sys.argv[1]
    if len(sys.argv) > 2:
        max_xi = float(sys.argv[2])
    else:
        max_xi = find_max_xi(model, updates)
    slices = (max_xi - min_xi)/50
    xis = []
    osmolarities = []
    exchanged = []
    medium_osmolarity = 280 * 0  # mM/L (estimate)
    q1 = (max_xi-min_xi)/4 + min_xi
    q2 = (max_xi+min_xi)/2
    q3 = 3*(max_xi-min_xi)/4 + min_xi
    out = []
    quarters = [min_xi, q1, q2, q3, max_xi]
    with Pool(processes=4) as pool:
        print('Started pool.')
        for block in range(0, 4):
            print('Starting process ' + str(block + 1))
            out.append(pool.apply_async(subprocess, (model, updates,
                       quarters[block], quarters[block + 1], slices,
                       medium_osmolarity, block,)))
        pool.close()
        pool.join()
    for output in out:
        xis.extend(output.get()['xi'])
        osmolarities.extend(output.get()['osmolarities'])
        exchanged.extend(output.get()['ex'])
    '''
        lac_present = False
        for entry in exchanged:
            if 'lac' in entry:
                lac_present = True
        assert lac_present
    '''
    bar = Bar('Drawing plots: ', max=len(exchanged))
    for mol in exchanged:
        bar.next()
        ex_over_ss = []
        for exchanges in osmolarities:
            if mol not in exchanges:
                exchanges[mol] = 0
            ex_over_ss.append(exchanges[mol])
        fig = figure()
        ax = fig.add_subplot(111)
        ax.set_title(mol + r' vs. $\xi$')
        ax.set_ylabel(mol)
        ax.set_xlabel(r'$\xi$')
        ax.scatter(xis, ex_over_ss)
        if not os.path.exists('Figures_' + enzyme_mass):
            os.mkdir('Figures_' + enzyme_mass)
        fig.savefig('Figures_' + enzyme_mass + '/' + mol)


def openJson(path):
    '''Shortened call to open JSON files.'''
    with open(path) as r:
        return json.load(r)


def flux_constraint(cobra_model, coefficients_forward, coefficients_reverse,
                    bound=1):
    """
    Adds a linear constrain of the form
        0 <= cf[1] vf[1] + cb[1] vb[1] + ... <= bounds
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
    cobra_model.solver.update()


def get_coefficients(model_updates, coef_forward, coef_backward, mult=1):
    for reaction in model_updates:
        coef_forward[reaction] = model_updates[reaction]['forward']*mult
        coef_backward[reaction] = model_updates[reaction]['backward']*mult


def fba_and_min_enzyme(cobra_model, coefficients_forward,
                       coefficients_reverse):
    """
    Performs FBA follows by minimization of enzyme content
    """
    with cobra_model as model:
        cobra.util.fix_objective_as_constraint(model)
        set_enzymatic_objective(model, coefficients_forward,
                                coefficients_reverse)
        return model.optimize()


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
    v_glc = .5  # mmol/gDW/h
    v_aa = 0.05

    # FIXME: Changing these values changes behavior of find_max_xi.
    # Bistable around min_xi = 0 and min_xi = 1200?
    # If v_glc and v_aa are bigger, program fails with xi = 0 infeasible...
    # TODO: On what time interval are uptakes defined in the model.
    IMDM = pandas.read_table('IMDM.txt', comment='#', sep='\s+')
    conc = {}
    for index, row in IMDM.iterrows():
        conc[row['id']] = row['mM']/1000  # mM -> M
    for met in conc:
        if met == 'glc_D_e':
            if xi == 0:
                    model.reactions.get_by_id('EX_glc_e_').lower_bound = -v_glc
            else:
                # TEST:

                model.reactions.get_by_id('EX_glc_e_').lower_bound = \
                    -min(v_glc, conc[met]/xi)
                '''
                model.reactions.get_by_id('EX_glc_e_').lower_bound = \
                    -conc[met]/xi
                '''
        else:
            try:
                name = met.split('_')[0]
                if xi == 0:
                    model.reactions.get_by_id(
                        'EX_' + name + '_e_').lower_bound = -v_aa
                else:
                    # TEST:

                    model.reactions.get_by_id(
                        'EX_' + name + '_e_').lower_bound = \
                        -min(v_aa, conc[met]/xi)
                    '''
                    model.reactions.get_by_id(
                        'EX_' + name + '_e_').lower_bound = -conc[met]/xi
                    '''
            except KeyError:
                if xi == 0:
                    model.reactions.get_by_id(
                        'EX_' + name + '_L_e_').lower_bound = -v_aa
                else:
                    model.reactions.get_by_id(
                        'EX_' + name + '_L_e_').lower_bound \
                        = -min(v_aa, conc[met]/xi)


def get_enzyme_mass(solution, coef_forward, coef_backward):
    enzyme_mass = 0
    for rxn, cf in coef_forward.items():
        enzyme_mass += cf * max(solution.fluxes[rxn], 0)
    for rxn, cb in coef_backward.items():
        enzyme_mass -= cb * min(solution.fluxes[rxn], 0)
    return enzyme_mass


def exchanges_consumption(cobra_model, solution):
    "Subset of exchanges that consume."
    import_list = [rxn for rxn in cobra_model.exchanges if
                   not rxn.products and solution.fluxes[rxn.id] < 0 or
                   not rxn.reactants and solution.fluxes[rxn.id] > 0]
    imports = {}
    for imp in import_list:
        imports[imp.id] = solution.fluxes[imp.id]
    return imports


def exchanges_secretion(cobra_model, solution):
    "Subset of exchanges that secrete"
    export_list = [rxn for rxn in cobra_model.exchanges if
                   not rxn.products and solution.fluxes[rxn.id] > 0 or
                   not rxn.reactants and solution.fluxes[rxn.id] < 0]
    exports = {}
    for exp in export_list:
        exports[exp.id] = solution.fluxes[exp.id]
    return exports


def osmolarity(model, solution):
    '''
    'Osmolarity' for one steady state.

    '''
    osmo = {}
    osmo['total'] = float(0)
    # CHANGED : checking for reaction id's
    permeable = ['EX_o2_e_', 'EX_co2_e_', 'EX_h2o_e_']
    imports = exchanges_consumption(model, solution)
    exports = exchanges_secretion(model, solution)
    for mol in imports:
        if imports[mol] < -0.001:
            osmo[mol] = imports[mol]
            if mol not in permeable:
                osmo['total'] += imports[mol]
    for mol in exports:
        if exports[mol] > 0.001:
            osmo[mol] = exports[mol]
            if mol not in permeable:
                osmo['total'] += exports[mol]
    return osmo


if __name__ == '__main__':
    assert len(sys.argv) > 1
    main()
    print()
