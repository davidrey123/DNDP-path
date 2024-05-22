# Created on : Mar 27, 2024, 5:00:04 PM
# Author     : michaellevin

INFTY = 1.0e9

class Params:   
 
    def __init__(self):
        self.bush_gap = 0.01
        self.pas_cost_mu = 0.001
        
        self.pas_cost_epsilon = 0.000001

        self.pas_flow_mu = 0.01
        self.flow_epsilon = 0.000001
        
        self.line_search_gap = 1E-7
        self.tapas_equilibrate_iter = 3
    
        self.DEBUG_CHECKS = True

        self.PRINT_PAS_INFO = False
        self.PRINT_BRANCH_INFO = False
        self.PRINT_TAPAS_INFO = False
        
        self.PRINT_TAP_ITER = False

        self.printBushEquilibrate = False
        self.printReducedCosts = False
        
        
        self.tapas_max_iter = 100
        self.min_gap = 1E-3
        self.msa_max_iter = 500
    
        self.warmstart = False
        
        self.PRINT_BB_INFO = False
        self.BB_timelimit = 600
        self.BB_tol = 1E-2
