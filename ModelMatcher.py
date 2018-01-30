'''
    This file contains a script which matches the new bigg reaction names in
    the iCHO model downloaded from BiGG to the older names in the iCHOv1_K1
    model from the paper.

    This allows us match data form BRENDA by reaction to the iCHOv1_K1 model.
    This is because the xml for this model did not contain EC codes, while the
    xml from BiGG did.
'''
import cobra
from bs4 import BeautifulSoup as Soup
from progress.bar import Bar


def MatchReaction(bid, reactants, products, reactions):
    ''' PSEUDOCODE:

    If for every metabolite there is a matching link, return true. Else
    return false.'''

    r_names = []
    for reactant in reactants:
        r_names.append('M_' + reactant.id)
    p_names = []
    for product in products:
        p_names.append('M_' + product.id)
    names = r_names + p_names
    for reaction in reactions:
        try:
            reactant_list = reaction.find('listOfReactants')
            products_list = reaction.find('listOfProducts')
            k1_reactants = reactant_list.find_all('speciesReference')
            k1_products = products_list.find_all('speciesReference')
            rids = []
            pids = []
            for reactant in k1_reactants:
                rids.append(reactant['species'])
            for product in k1_products:
                pids.append(product['species'])
            k1_metabolites = rids + pids
            if k1_metabolites == names:
                return reaction['id']
        except AttributeError:
            pass


if __name__ == '__main__':
    print('Loading models....')
    bigg_model = cobra.io.read_sbml_model('iCHOv1.xml')
    k1_model = cobra.io.read_sbml_model('iCHOv1_K1_final.xml')
    print('Parsing xml into a beautiful soup....')
    handler = open('iCHOv1.xml').read()
    bigg_xml = Soup(handler, 'xml')
    handler = open('iCHOv1_K1_final.xml').read()
    k1_xml = Soup(handler, 'xml')
    v1_to_k1 = {}
    # xml_bigg_reactions = bigg_xml.find_all('reaction')
    xml_k1_reactions = k1_xml.find_all('reaction')
    total = len(bigg_model.reactions)
    bar = Bar('Matching models ', max=total)
    for reaction in bigg_model.reactions:
        if reaction.id not in k1_model.reactions:
            reactants_bigg = reaction.reactants
            products_bigg = reaction.products
            v1_to_k1[reaction.id] = MatchReaction(reaction.id, reactants_bigg,
                                                  products_bigg,
                                                  xml_k1_reactions)
        bar.next()
