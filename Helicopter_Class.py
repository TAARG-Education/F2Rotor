class Helicopter_HS:
    def __init__(self, tau, b_w, AR, R_LS, Hf, r, gamma, k, 
                 R_w_b, c_w, R_mr, ala_fissa, profilo, turbolentflow):
        

        self.tau = tau
        self.b_w = b_w
        self.AR = AR
        self.R_LS = R_LS
        self. Hf = Hf
        self.r = r
        self.gamma = gamma
        self.k = k
        self.R_w_b = R_w_b
        self.c_w=c_w
        self.R_mr=R_mr
        self.M_inf = None
        self.Re_inf = None
        self.ala_fissa = ala_fissa
        self.profilo = profilo
        self.turbolentflow = turbolentflow

        

    def set_Re_inf(self, v, rho, nu):
        self.Re_inf = (v * self.c_w) / nu

    def set_M_inf(self, v, a):
        self.M_inf = v/a
        

   
        
        

