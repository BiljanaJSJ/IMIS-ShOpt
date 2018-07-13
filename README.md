
This paper proposes a general optimization strategy, which combines results from different optimization or parameter estimation methods to overcome shortcomings of a single method.  
Shotgun optimization is developed as a framework which employs different optimization strategies, criteria, or conditional targets to enable wider likelihood exploration.  The introduced Shotgun optimization approach is embedded into an incremental mixture importance sampling algorithm to produce improved posterior samples for multimodal densities and creates robustness in cases where the likelihood and prior are in disagreement.  Despite using different optimization approaches, the samples are combined into samples from a single target posterior.   The diversity of the framework is demonstrated on parameter estimation from differential equation models employing diverse strategies including numerical solutions and approximations thereof.  Additionally the approach is demonstrated on mixtures of discrete and continuous parameters and is shown to ease estimation from synthetic likelihood models. 


This the code for the IMIS-ShOpt Examples.





Brief descriptions of the files included in this directory:

FhN_fullModel_IMIS_ShOpt               - example in the Section 4.1.1. 
FhN_One_IMIS_ShOpt_IMIS_Opt            - example in the Sections 4.1.2.
FhN_One_IMIS_ShOpt_IMIS_Opt_NLS        - example in the Sections 4.1.2.
FhN_One_IMIS_ShOpt_IMIS_Opt_Profiling  - example in the Sections 4.1.2.
FhN_One_IMIS_ShOpt_IMIS_Opt_TwoStage   - example in the Sections 4.1.2.
SIR_IMIS_Opt                           - example in the Section 4.2.
SIR_IMIS_ShOpt                         - example in the Section 4.2.
theta_Ricker_IMIS_ShOpt_SL             - example in the Section 4.3.