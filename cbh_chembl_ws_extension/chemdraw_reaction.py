import xmltodict

def get_keys(x):
	dataKeys = {}

	for component in x['CDXML']['page']['stoichiometrygrid']['sgcomponent']:
		if component.get('@ComponentIsHeader', False):
			for datum in component['sgdatum']:
				dataKeys[datum['@SGPropertyType']] = datum['@SGDataValue']
	return dataKeys


def compounds(x, dataKeys,  product_ids, reagent_ids, reactant_ids):
	for component in x['CDXML']['page']['stoichiometrygrid']['sgcomponent']:
		if not component.get('@ComponentIsHeader', False):
			role = None
			if component['@ComponentReferenceID'] in product_ids:
				role = 'product'
			elif component['@ComponentReferenceID'] in reagent_ids:
				role = 'reagent'
			elif component['@ComponentReferenceID'] in reagent_ids:
				role = 'reagent'
			dicttoyield = {	
				dataKeys[datum['@SGPropertyType']]: 
				datum['@SGDataValue'] for datum in component['sgdatum']
			}
			dicttoyield['role'] = role
			yield dicttoyield


def parse(xml_path):
	with open(xml_path) as xfile:
		xml = xfile.read()
	x = xmltodict.parse(xml)
	if x['CDXML']['page'].get('scheme'):
		reactant_ids = x['CDXML']['page']['scheme']['step']['@ReactionStepReactants'].split(' ')
		reagent_ids = x['CDXML']['page']['scheme']['step']['@ReactionStepObjectsAboveArrow'].split(' ')
		product_ids = x['CDXML']['page']['scheme']['step']['@ReactionStepProducts'].split(' ')
		keys = get_keys(x)
		return  [p for p in compounds(x,keys,  product_ids, reagent_ids, reactant_ids)]
	print('not a reaction')
	return []

