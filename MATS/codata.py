#Constants
CODATA = {
    'h' : {
        'value': 6.62607015e-27,
        'units': 'erg s',
        'source' : 'https://physics.nist.gov/cgi-bin/cuu/Value?h|search_for=h',
        'date': '5/21/2020',
        'name' : 'Planks constant',
    },
    'c' : {
        'value': 29979245800,
        'units' : 'cm/s',
        'source': "https://physics.nist.gov/cgi-bin/cuu/Value?c|search_for=c",
        'date': "5/21/2020",
        'name': 'speed of light',
    },
    'k' : {
        'value': 1.380649e-16,
        'units': "erg / K",
        'source': "https://physics.nist.gov/cgi-bin/cuu/Value?k",
        'date': "5/21/2020",
        'name': 'Boltzmann constant',
    },
    'Na' : {
        'value': 6.02214076e23,
        'units': 'mol-1',
        'source' : "https://physics.nist.gov/cgi-bin/cuu/Value?na",
        'date': '5/21/2020',
        'name': 'Avogadro number',
    },
    'cpa_atm' : {
        'value': (10*101325)**-1,
        'units': 'atm/cpa',
        'source': 'https://physics.nist.gov/cgi-bin/cuu/Value?stdatm|search_for=atmosphere',
        'date': '5/21/2020',
        'name': 'convert cpa to atm',
    },
}

CONSTANTS = {
    k : v['value'] for k, v in CODATA.items()
}
CONSTANTS['c2'] = CONSTANTS['h'] * CONSTANTS['c'] / CONSTANTS['k']
