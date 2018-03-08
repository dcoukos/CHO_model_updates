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


def main():
    k1_model = cobra.io.read_sbml_model('iCHOv1_K1_final.xml')
    k1_updates = openJson('JSONs/k1_updates.json')
    population_osmolarities(k1_model, k1_updates)


def find_max_xi(model, coef_forward, coef_backward, tolerance=1):
    """Returns the max ξ value for which the model is still feasible through a
        bissection algorithm.

    Parameters
    ----------
    model : cobra.model
    coef_forward : dict
        Enzyme costs of forward reactions.
    coef_backward: dict
        Enzyme costs of reverse reactions.
    tolerance : tolerance defining when to stop searching
        difference = high value - low value. If difference < tolerance, return
        low value.

    Returns
    -------
    float
        Highest ξ for which the model is feasible.
    """
    assert tolerance > 0
    min_xi = 0.
    max_xi = 1000.
    max_undefined = True
    print("finding an infeasible xi")
    while max_undefined:
        print('.')
        sol = run_fba(model, coef_forward, coef_backward, max_xi)
        if sol.status == 'optimal':
            run_fba(model, coef_forward, coef_backward, max_xi)
            max_xi *= 2
        else:
            assert sol.status == 'infeasible'
            max_undefined = False
    print("\nfinding max xi", end='')
    while max_xi - min_xi > tolerance:
        assert min_xi <= max_xi
        avg_xi = (max_xi + min_xi)/2
        print('min_xi=', min_xi, ', max_xi=', max_xi, ', avg_xi=', avg_xi,
              flush=True)
        try:
            run_fba(model, coef_forward, coef_backward, avg_xi)
            min_xi = avg_xi
        except Infeasible:
            max_xi = avg_xi
    print('Returning ' + str(min_xi))
    return float(min_xi)


def run_fba(model, coef_forward, coef_backward, xi):
    """Adds the constraints to the model, releases bounds, and sets atp
        maintenance before minimization.

    Parameters
    ----------
    model : cobra.model
    coef_forward : dict
        Enzyme costs for the forward direction.
    coef_backward : dict
        Enzyme costs for the backward direction.
    xi : float
        xi = D/X.

    Returns
    -------
    cobra.solution
        Solution object of the minimized model, to which enzymatic constraints
        have been applied.
    """
    enzyme_mass = float(sys.argv[1])
    release_bounds(model)  # Give the model wiggle room.
    constrain_uptakes(model, xi)  # Model competition between cells.
    model.reactions.DM_atp_c_.lower_bound = 1  # atp maintenance.
    flux_constraint(model, coef_forward, coef_backward, enzyme_mass)
    # m_coefs = mitochondrial_coefs(model, coef_forward, coef_backward)
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


def subprocess(model, coef_forward, coef_backward, min_xi, max_xi, slices,
               medium_osmolarity, process=4):
    """Main subprocesses where the model is optimized and osmolarity
        contributions are calculated based on a given sub-range of ξ values.

    Parameters
    ----------
    model : cobra.Model
    coef_forward : dictionary
        Enzyme costs of forward reactions
    coef_backward : dictionary
        Enzyme costs of backwards reactions
    min_xi : type
        Description of parameter `min_xi`.
    max_xi : type
        Description of parameter `max_xi`.
    slices : type
        Description of parameter `slices`.
    medium_osmolarity : type
        Description of parameter `medium_osmolarity`.
    process : type
        Description of parameter `process`.

    Returns
    -------
    type
        Description of returned object.

    """
    if process == 1:
        bar = Bar('Calculating osmolarities: ', max=slices)
    exchanged = []
    osmolarities = []
    xis = []
    for xi in numpy.linspace(min_xi, max_xi, slices):
        if process == 1:
            bar.next()
        try:
            sol = run_fba(model, coef_forward, coef_backward, xi)
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
    """Finds max feasible ξ, optimizes the model within this range, and plots
        osmolarity contributions along these ξ values.
    """
    # 1st infeasible value = 1202.9
    # max_xi = 1200/3
    enzyme_mass = sys.argv[1]
    coef_forward = {}
    coef_backward = {}
    get_coefficients(updates, coef_forward, coef_backward, 3.6E6)
    if len(sys.argv) > 2:
        max_xi = float(sys.argv[2])
    else:
        max_xi = find_max_xi(model, coef_forward, coef_backward)
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
            out.append(pool.apply_async(subprocess, (model, coef_forward,
                       coef_backward, quarters[block], quarters[block + 1],
                       slices, medium_osmolarity, block,)))
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


def get_coefficients(model_updates, coef_forward, coef_backward, mult=1.):
    """Reads specific activities and fills coeficient dictionaries.

    Parameters
    ----------
    model_updates : cobra.Model
    coef_forward : empty dict
        Dictionary of enzyme costs of forward reactions.
    coef_backward : empty dict
        Dictionary of enzyme costs of backward reactions.
    mult : int
        mu

    Returns
    -------
    type
        Description of returned object.

    """
    # TODO: Is this division right?
    for reaction in model_:
        coef_forward[reaction] = mult/model_updates[reaction]['forward']
        coef_backward[reaction] = mult/model_updates[reaction]['backward']


def fba_and_min_enzyme(cobra_model, coefficients_forward,
                       coefficients_reverse):
    """
    Performs FBA follows by minimization of enzyme content
    """
    with cobra_model as model:
        model.objective = model.reactions.biomass_cho_producing
        # TODO: Should the solution object be collected here?
        model.optimize(objective_sense='maximize')
        model.reactions.biomass_cho_producing.lower_bound = \
            model.reactions.biomass_cho_producing.flux
        # cobra.util.fix_objective_as_constraint(model)
        set_enzymatic_objective(model, coefficients_forward,
                                coefficients_reverse)
        return model.optimize(objective_sense='minimize')


def release_bounds(model):
    """The iCHO model has experimental uptakes rigidly set to experimental
        observations. Since our experimental conditions and modeling
        constraints are different, it is important to give the model some
        wiggle room.

    Parameters
    ----------
    model : cobra.Model

    """
    for rxn in model.reactions:
        if rxn.upper_bound == rxn.lower_bound:
            if rxn.upper_bound < 0:
                rxn.upper_bound = 0
            if rxn.upper_bound > 0:
                rxn.lower_bound = 0


def set_enzymatic_objective(cobra_model, coefficients_forward,
                            coefficients_reverse):
    """Set enzyme cost as a constraint for the cobra model, and tells it to
        minimize this cost.

    Parameters
    ----------
    cobra_model : cobra.Model
    coefficients_forward : dict
        Enzyme costs for forward reactions.
    coefficients_reverse : dict
        Enzyme costs for reverse reactions.
    """
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
    """When the cells are not competing (resource excess, high ξ) uptakes are
        restricted by the limits of metabolite importers. When competion occurs
        competiton is modeled by limiting uptakes to metabolite_concentration/ξ
    Parameters
    ----------
    model : cobra.Model
    xi : float
        = D/X

    """

    v_glc = .5  # mmol/gDW/h
    v_aa = 0.05
    IMDM = pandas.read_table('IMDM.txt', comment='#', sep='\s+')
    conc = {}
    for index, row in IMDM.iterrows():
        conc[row['id']] = row['mM']
    for met in conc:
        if met == 'glc_D_e':
            if xi == 0:
                    model.reactions.get_by_id('EX_glc_e_').lower_bound = -v_glc
            else:
                model.reactions.get_by_id('EX_glc_e_').lower_bound = \
                    -min(v_glc, conc[met]/xi)
        else:
            try:
                name = met.split('_')[0]
                if xi == 0:
                    model.reactions.get_by_id(
                        'EX_' + name + '_e_').lower_bound = -v_aa
                else:
                    model.reactions.get_by_id(
                        'EX_' + name + '_e_').lower_bound = \
                        -min(v_aa, conc[met]/xi)
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
    """Get fluxes of reactions that consume metabolites. Only considers
        exchanges

    Parameters
    ----------
    cobra_model : cobra.Model
    solution : cobra.Solution
        Solution of optimized model.

    Returns
    -------
    dictionary
        keys are bigg ids of import reactions, values are fluxes.
    """
    import_list = [rxn for rxn in cobra_model.exchanges if
                   not rxn.products and solution.fluxes[rxn.id] < 0 or
                   not rxn.reactants and solution.fluxes[rxn.id] > 0]
    imports = {}
    for imp in import_list:
        imports[imp.id] = solution.fluxes[imp.id]
    return imports


def exchanges_secretion(cobra_model, solution):
    """Get fluxes of reactions that export metabolites. Only considers
            exchanges

    Parameters
    ----------
    cobra_model : cobra.Model
    solution : cobra.Solution
        Solution of optimized model.

    Returns
    -------
    dictionary
        keys are bigg ids of import reactions, values are fluxes.
    """
    export_list = [rxn for rxn in cobra_model.exchanges if
                   not rxn.products and solution.fluxes[rxn.id] > 0 or
                   not rxn.reactants and solution.fluxes[rxn.id] < 0]
    exports = {}
    for exp in export_list:
        exports[exp.id] = solution.fluxes[exp.id]
    return exports


def osmolarity(model, solution):
    ''' Calculates a cell's contribution to osmolarity for a given steady
        state defined by ξ and for an optmized model (thru solution object)
    '''
    osmo = {}
    osmo['total'] = float(0)
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
