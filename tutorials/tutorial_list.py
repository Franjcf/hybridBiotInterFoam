# tutorials list for python validation script
# adapted from porousMultiphaseFoam source code, Horgue P. et al. (2014)

tutorials = [
             
	     # Elastic Cases
             {'solver': "elasticHBIF", 'category' : "elastic_cases", 'case' : "terzaghi"},
             {'solver': "elasticHBIF", 'category' : "elastic_cases", 'case' : "oscillating_core"}, \
             {'solver': "elasticHBIF", 'category' : "elastic_cases", 'case' : "coastal_barrier"}, \
            
	     # Plastic Cases
             {'solver': "plasticHBIF", 'category' : "plastic_cases", 'case' : "hele_shaw_fracturing"}, \
             {'solver': "plasticHBIF", 'category' : "plastic_cases", 'case' : "fracture_wellbore_pressure"}, \
             {'solver': "plasticHBIF", 'category' : "plastic_cases", 'case' : "surface_deformation"}

	    ]
