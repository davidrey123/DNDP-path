# Created on : Mar 27, 2024, 5:00:04 PM
# Author     : michaellevin

class Params:
 
    def __init__(self):
        self.bush_gap = 0.001
        self.pas_cost_mu = 0.000005
        self.pas_flow_mu = 0.01
        self.flow_epsilon = 0.001
        self.line_search_gap = 0.000001
        self.tapas_equilibrate_iter = 3
    
        self.DEBUG_CHECKS = True

        self.PRINT_PAS_INFO = False
        self.PRINT_TAPAS_INFO = True

        self.INFTY = 1.0e9

        self.printBushEquilibrate = False
        self.printReducedCosts = False
    
    
    