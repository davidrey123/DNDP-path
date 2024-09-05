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
        
        self.tapas_max_iter = 10
        self.min_gap = 1E-3
        self.min_gap_SO_OA_cuts = 1E-1 #---warning should not be used to compute lower bounds on TSTT - only OA cuts
        
        self.msa_max_iter = 500
    
        self.warmstart = False
        
        self.PRINT_BB_INFO = True #---prints detailed BB info
        self.PRINT_BB_BASIC = True #---prints only basic BB info
        
        self.CPLEX_threads = 4
        self.BB_timelimit = 3600
        self.BB_tol = 1E-2

        self.OABPC_tol = 1E-2
        self.KNP = False
        
        self.useAONcuts = False
        
        self.runUEifCGIntegral = True
        self.useInterdictionCuts = True
