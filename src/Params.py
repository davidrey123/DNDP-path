# Created on : Mar 27, 2024, 5:00:04 PM
# Author     : michaellevin

INFTY = 1.0e9

class Params:   
    
    def resetPAS(self):
        # if these are too large, then we avoid flow shifting on PAS
        # if these are too small, then we use PAS that are ineffective
        self.pas_cost_mu = 0.01
        self.pas_flow_mu = 0.01
        
        # if this is too large, then we avoid flow shifting in PAS
        # if this is too small, we waste time shifting in PAS that is useless
        self.pas_cost_epsilon = 0.01

        
        self.line_search_gap = 1E-2
        
        
        
        
    def __init__(self):
        self.bush_gap = 0.01
        self.pas_cost_mu = 0
        
        self.pas_cost_epsilon = 0

        self.pas_flow_mu = 0
        
        self.line_search_gap = 1E-3
        
        self.resetPAS()
        
        self.flow_epsilon = 0.000001
        
        
        self.min_line_search_gap = 1E-6
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
        self.min_gap_SO_OA_cuts = 1E-1 #---warning should not be used to compute lower bounds on TSTT - only OA cuts
        
        self.msa_max_iter = 500
    
        self.warmstart = False
        
        self.PRINT_BB_INFO = True #---prints detailed BB info
        self.PRINT_BB_BASIC = False #---prints only basic BB info
        
        self.CPLEX_threads = 1
        self.BB_timelimit = 3600
        self.BB_tol = 1E-2

        self.OABPC_tol = 1E-2
        self.KNP = False
        
        self.useAONcuts = False
        
        self.runUEifCGIntegral = True
        self.useInterdictionCuts = True
        
        self.OAcut_tol = 0.01
        
        
        
        # used within TAPAS don't change
        
        self.good_pas_cost_mu = 0
        self.good_pas_flow_mu = 0
        self.good_bush_gap = 0
        self.good_pas_cost_epsilon = 0
        
    
