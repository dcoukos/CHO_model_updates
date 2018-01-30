import cobra


model = cobra.io.read_sbml_model('iCHOv1.xml')

flux_expression = 0
for bigg_id in model_updates:
    assert bigg_id in model.reactions
    assert 0 <= coefficients_backward[bigg_id] < float('inf')
    assert 0 <= coefficients_forward[rxn] < float('inf')
    flux_expression += coefficients_backward[bigg_id] * bigg_id.backward_variable
    flux_expression += coefficients_forward[bigg_id] * bigg_id.forward_variable

cobra_model.problem.Constraint(flux_expression, lb=0, ub=1)
