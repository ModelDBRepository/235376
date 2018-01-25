
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                           %
%            			Model of Cerebellar PF-PC 						    %
%	       			Long-Term Depression & Potentiation                     %
%                                                                           %
%                     developed by Andrew R. Gallimore                      %  
%          						December 2017      						    %	
%							for Matlab Simbiology	                        %
%                                                                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%TO SIMULATE THIS MODEL %%%%%%%%%%%%%%%%%%%%%%
%To simulate a standard LTD or LTP protocol, run the RUN_LTD_LTP.m file

global pfpc_model;
% Create model object
pfpc_model = sbiomodel('pfpc_model');
% Add a compartment named cell to model
cytosol = addcompartment(pfpc_model, 'cytosol');
cytosol.Capacity = 1;
cytosol.CapacityUnits = 'microliter';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%CALCIUM MODEL%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Add species
Ca = addspecies(cytosol, 'Ca', 'InitialAmount', 0.06, 'InitialAmountUnits', 'micromolarity');
PMCA = addspecies(cytosol, 'PMCA', 'InitialAmount', 0.2076, 'InitialAmountUnits', 'micromolarity');
CaPMCA = addspecies(cytosol, 'CaPMCA', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');

NCX = addspecies(cytosol, 'NCX', 'InitialAmount', 0.0623, 'InitialAmountUnits', 'micromolarity');
CaNCX = addspecies(cytosol, 'CaNCX', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
Ca2NCX = addspecies(cytosol, 'Ca2NCX', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');

SERCA = addspecies(cytosol, 'SERCA', 'InitialAmount', 2.0757, 'InitialAmountUnits', 'micromolarity');
CaSERCA = addspecies(cytosol, 'CaSERCA', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
Ca2SERCA = addspecies(cytosol, 'Ca2SERCA', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');

PV = addspecies(cytosol, 'PV', 'InitialAmount', 1.15, 'InitialAmountUnits', 'micromolarity');
CaPV = addspecies(cytosol, 'CaPV', 'InitialAmount', 8.4, 'InitialAmountUnits', 'micromolarity');
Ca2PV = addspecies(cytosol, 'Ca2PV', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
MgPV = addspecies(cytosol, 'MgPV', 'InitialAmount', 30.45, 'InitialAmountUnits', 'micromolarity');
Mg2PV = addspecies(cytosol, 'Mg2PV', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');

Mg = addspecies(cytosol, 'Mg', 'InitialAmount', 590, 'InitialAmountUnits', 'micromolarity');
CBs = addspecies(cytosol, 'CBs', 'InitialAmount', 36.25, 'InitialAmountUnits', 'micromolarity');
CaCBs = addspecies(cytosol, 'CaCBs', 'InitialAmount', 3.4, 'InitialAmountUnits', 'micromolarity');
Ca2CBs = addspecies(cytosol, 'Ca2CBs', 'InitialAmount', 0.125, 'InitialAmountUnits', 'micromolarity');

CBf = addspecies(cytosol, 'CBf', 'InitialAmount', 37.77, 'InitialAmountUnits', 'micromolarity');
CaCBf = addspecies(cytosol, 'CaCBf', 'InitialAmount', 2.1, 'InitialAmountUnits', 'micromolarity');
Ca2CBf = addspecies(cytosol, 'Ca2CBf', 'InitialAmount', 0.125, 'InitialAmountUnits', 'micromolarity');

Ca_ER = addspecies(cytosol, 'Ca_ER', 'InitialAmount', 100, 'InitialAmountUnits', 'micromolarity');


%%%%%%Ca model parameters%%%%%%
kon_pmca_ca = addparameter(pfpc_model,'kon_pmca_ca', 'Value', 25000, 'ValueUnits', '1/(micromolarity*second)');
koff_pmca_ca = addparameter(pfpc_model,'koff_pmca_ca', 'Value', 2000, 'ValueUnits', '1/second');
kflux_pmca = addparameter(pfpc_model,'kflux_pmca', 'Value', 500, 'ValueUnits', '1/second');

kon_ncx_ca = addparameter(pfpc_model,'kon_ncx_ca', 'Value', 93.827, 'ValueUnits', '1/(micromolarity*second)');
koff_ncx_ca = addparameter(pfpc_model,'koff_ncx_ca', 'Value', 612.6, 'ValueUnits', '1/second');
kflux_ncx = addparameter(pfpc_model,'kflux_ncx', 'Value', 1000, 'ValueUnits', '1/second');

kon_serca_ca = addparameter(pfpc_model,'kon_serca_ca', 'Value', 17147, 'ValueUnits', '1/(micromolarity*second)');
koff_serca_ca = addparameter(pfpc_model,'koff_serca_ca', 'Value', 8426.3, 'ValueUnits', '1/second');
kflux_serca = addparameter(pfpc_model,'kflux_serca', 'Value', 2.60, 'ValueUnits', '1/second');

kleak = addparameter(pfpc_model,'kleak', 'Value', 39.4378, 'ValueUnits', 'micromole/second');

kon_pv_ca = addparameter(pfpc_model,'kon_pv_ca', 'Value', 107, 'ValueUnits', '1/(micromolarity*second)');
koff_pv_ca = addparameter(pfpc_model,'koff_pv_ca', 'Value', 0.95, 'ValueUnits', '1/second');

kon_pv_mg = addparameter(pfpc_model,'kon_pv_mg', 'Value', 472, 'ValueUnits', '1/(micromolarity*second)');
koff_pv_mg = addparameter(pfpc_model,'koff_pv_mg', 'Value', 25, 'ValueUnits', '1/second');

kon_cbs_ca = addparameter(pfpc_model,'kon_cbs_ca', 'Value', 5.5, 'ValueUnits', '1/(micromolarity*second)');
koff_cbs_ca = addparameter(pfpc_model,'koff_cbs_ca', 'Value', 2.6, 'ValueUnits', '1/second');

kon_cbf_ca = addparameter(pfpc_model,'kon_cbf_ca', 'Value', 43.5, 'ValueUnits', '1/(micromolarity*second)');
koff_cbf_ca = addparameter(pfpc_model,'koff_cbf_ca', 'Value', 35.8, 'ValueUnits', '1/second');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%INPUTS FOR CALCIUM & GLUTAMATE%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
kinflux = addparameter(pfpc_model,'kinflux', 'Value', 0, 'ValueUnits', 'micromole/second', 'ConstantValue', false);
kinflux_pf = addparameter(pfpc_model,'kinflux_pf', 'Value', 0, 'ValueUnits', 'micromole/second', 'ConstantValue', false);
kinflux_glu = addparameter(pfpc_model,'kinflux_glu', 'Value', 0, 'ValueUnits', 'micromole/second', 'ConstantValue', false);

%%%%Calcium Influx and Buffering Reactions%%%%

% Ca + PMCA <->  Ca1PMCA -> PMCA (efflux pump)
ca_reac_1 = addreaction(pfpc_model, 'Ca + PMCA -> CaPMCA')
ca_reac_1_k = addkineticlaw(ca_reac_1,'MassAction');
set(ca_reac_1_k, 'ParameterVariableNames', {'kon_pmca_ca'});

ca_reac_2 = addreaction(pfpc_model, 'CaPMCA -> Ca + PMCA')
ca_reac_2_k = addkineticlaw(ca_reac_2,'MassAction');
set(ca_reac_2_k, 'ParameterVariableNames', {'koff_pmca_ca'});

ca_reac_3 = addreaction(pfpc_model, 'CaPMCA -> PMCA')
ca_reac_3_k = addkineticlaw(ca_reac_3,'MassAction');
set(ca_reac_3_k, 'ParameterVariableNames', {'kflux_pmca'});

% Ca + NCX <->  Ca1NCX + Ca <->  Ca2NCX -> NCX
ca_reac_4 = addreaction(pfpc_model, 'Ca + NCX -> CaNCX')
ca_reac_4_k = addkineticlaw(ca_reac_4,'MassAction');
set(ca_reac_4_k, 'ParameterVariableNames', {'kon_ncx_ca'});

ca_reac_5 = addreaction(pfpc_model, 'CaNCX -> Ca + NCX')
ca_reac_5_k = addkineticlaw(ca_reac_5,'MassAction');
set(ca_reac_5_k, 'ParameterVariableNames', {'koff_ncx_ca'});

ca_reac_6 = addreaction(pfpc_model, 'CaNCX + Ca -> Ca2NCX')
ca_reac_6_k = addkineticlaw(ca_reac_6,'MassAction');
set(ca_reac_6_k, 'ParameterVariableNames', {'kon_ncx_ca'});

ca_reac_7 = addreaction(pfpc_model, 'Ca2NCX -> CaNCX + Ca')
ca_reac_7_k = addkineticlaw(ca_reac_7,'MassAction');
set(ca_reac_7_k, 'ParameterVariableNames', {'koff_ncx_ca'});

ca_reac_8 = addreaction(pfpc_model, 'Ca2NCX -> NCX')
ca_reac_8_k = addkineticlaw(ca_reac_8,'MassAction');
set(ca_reac_8_k, 'ParameterVariableNames', {'kflux_ncx'});


% Ca + SERCA <->  Ca1SERCA +Ca <->  Ca2SERCA  ->  SERCA
ca_reac_9 = addreaction(pfpc_model, 'Ca + SERCA -> CaSERCA')
ca_reac_9_k = addkineticlaw(ca_reac_9,'MassAction');
set(ca_reac_9_k, 'ParameterVariableNames', {'kon_serca_ca'});

ca_reac_10 = addreaction(pfpc_model, 'CaSERCA -> Ca + SERCA')
ca_reac_10_k = addkineticlaw(ca_reac_10,'MassAction');
set(ca_reac_10_k, 'ParameterVariableNames', {'koff_serca_ca'});

ca_reac_11 = addreaction(pfpc_model, 'CaSERCA + Ca -> Ca2SERCA')
ca_reac_11_k = addkineticlaw(ca_reac_11,'MassAction');
set(ca_reac_11_k, 'ParameterVariableNames', {'kon_serca_ca'});

ca_reac_12 = addreaction(pfpc_model, 'Ca2SERCA -> CaSERCA + Ca')
ca_reac_12_k = addkineticlaw(ca_reac_12,'MassAction');
set(ca_reac_12_k, 'ParameterVariableNames', {'koff_serca_ca'});

ca_reac_13 = addreaction(pfpc_model, 'Ca2SERCA -> SERCA + Ca_ER + Ca_ER')
ca_reac_13_k = addkineticlaw(ca_reac_13,'MassAction');
set(ca_reac_13_k, 'ParameterVariableNames', {'kflux_serca'});


% Leak
ca_reac_14 = addreaction(pfpc_model, 'null -> Ca')
ca_reac_14_k = addkineticlaw(ca_reac_14,'MassAction');
set(ca_reac_14_k, 'ParameterVariableNames', {'kleak'});


% PV + Ca <-> CaPV + Ca <-> Ca2PV
ca_reac_15 = addreaction(pfpc_model, 'Ca + PV -> CaPV')
ca_reac_15_k = addkineticlaw(ca_reac_15,'MassAction');
set(ca_reac_15_k, 'ParameterVariableNames', {'kon_pv_ca'});

ca_reac_16 = addreaction(pfpc_model, 'CaPV -> Ca + PV')
ca_reac_16_k = addkineticlaw(ca_reac_16,'MassAction');
set(ca_reac_16_k, 'ParameterVariableNames', {'koff_pv_ca'});

ca_reac_17 = addreaction(pfpc_model, 'CaPV + Ca -> Ca2PV')
ca_reac_17_k = addkineticlaw(ca_reac_17,'MassAction');
set(ca_reac_17_k, 'ParameterVariableNames', {'kon_pv_ca'});

ca_reac_18 = addreaction(pfpc_model, 'Ca2PV -> CaPV + Ca')
ca_reac_18_k = addkineticlaw(ca_reac_18,'MassAction');
set(ca_reac_18_k, 'ParameterVariableNames', {'koff_pv_ca'});


% PV + Mg <-> MgPV <-> Mg2PV, apparent rate constants, [Mg] = 590uM
ca_reac_19 = addreaction(pfpc_model, 'PV + Mg -> MgPV')
ca_reac_19_k = addkineticlaw(ca_reac_19,'MassAction');
set(ca_reac_19_k, 'ParameterVariableNames', {'kon_pv_mg'});

ca_reac_20 = addreaction(pfpc_model, 'MgPV -> PV + Mg')
ca_reac_20_k = addkineticlaw(ca_reac_20,'MassAction');
set(ca_reac_20_k, 'ParameterVariableNames', {'koff_pv_mg'});

ca_reac_21 = addreaction(pfpc_model, 'MgPV + Mg -> Mg2PV')
ca_reac_21_k = addkineticlaw(ca_reac_21,'MassAction');
set(ca_reac_21_k, 'ParameterVariableNames', {'kon_pv_mg'});

ca_reac_22 = addreaction(pfpc_model, 'Mg2PV -> MgPV')
ca_reac_22_k = addkineticlaw(ca_reac_22,'MassAction');
set(ca_reac_22_k, 'ParameterVariableNames', {'koff_pv_mg'});


% CBs + Ca <-> CaCBs + Ca <-> Ca2CBs, CBs - high affinity site
ca_reac_23 = addreaction(pfpc_model, 'Ca + CBs -> CaCBs')
ca_reac_23_k = addkineticlaw(ca_reac_23,'MassAction');
set(ca_reac_23_k, 'ParameterVariableNames', {'kon_cbs_ca'});

ca_reac_24 = addreaction(pfpc_model, 'CaCBs -> Ca + CBs')
ca_reac_24_k = addkineticlaw(ca_reac_24,'MassAction');
set(ca_reac_24_k, 'ParameterVariableNames', {'koff_cbs_ca'});

ca_reac_25 = addreaction(pfpc_model, 'CaCBs + Ca -> Ca2CBs')
ca_reac_25_k = addkineticlaw(ca_reac_25,'MassAction');
set(ca_reac_25_k, 'ParameterVariableNames', {'kon_cbs_ca'});

ca_reac_26 = addreaction(pfpc_model, 'Ca2CBs -> CaCBs + Ca')
ca_reac_26_k = addkineticlaw(ca_reac_26,'MassAction');
set(ca_reac_26_k, 'ParameterVariableNames', {'koff_cbs_ca'});


% CBf + Ca <-> CaCBf + Ca <-> Ca2CBf, CBf - medium affinity site
ca_reac_27 = addreaction(pfpc_model, 'Ca + CBf -> CaCBf')
ca_reac_27_k = addkineticlaw(ca_reac_27,'MassAction');
set(ca_reac_27_k, 'ParameterVariableNames', {'kon_cbf_ca'});

ca_reac_28 = addreaction(pfpc_model, 'CaCBf -> Ca + CBf')
ca_reac_28_k = addkineticlaw(ca_reac_28,'MassAction');
set(ca_reac_28_k, 'ParameterVariableNames', {'koff_cbf_ca'});

ca_reac_29 = addreaction(pfpc_model, 'CaCBf + Ca -> Ca2CBf')
ca_reac_29_k = addkineticlaw(ca_reac_29,'MassAction');
set(ca_reac_29_k, 'ParameterVariableNames', {'kon_cbf_ca'});

ca_reac_30 = addreaction(pfpc_model, 'Ca2CBf -> CaCBf + Ca')
ca_reac_30_k = addkineticlaw(ca_reac_30,'MassAction');
set(ca_reac_30_k, 'ParameterVariableNames', {'koff_cbf_ca'});

%%%mGluR Model%%%
%%%%%mGluR Model%%%%%
%Species (uM)
Ca_deg = addspecies(cytosol, 'Ca_deg', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
glu = addspecies(cytosol, 'glu', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');%250 when released
glu_deg = addspecies(cytosol, 'glu_deg', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');%250 when released
mGLUR = addspecies(cytosol, 'mGLUR', 'InitialAmount', 8.3333, 'InitialAmountUnits', 'micromolarity');
mGLUR_gq = addspecies(cytosol, 'mGLUR_gq', 'InitialAmount', 6.6667, 'InitialAmountUnits', 'micromolarity');
mGLUR_glu = addspecies(cytosol, 'mGLUR_glu', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
mGLUR_gq_glu = addspecies(cytosol, 'mGLUR_gq_glu', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
gq_gdp = addspecies(cytosol, 'gq_gdp', 'InitialAmount', 43.333, 'InitialAmountUnits', 'micromolarity');
ga_gdp = addspecies(cytosol, 'ga_gdp', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
ga_gtp = addspecies(cytosol, 'ga_gtp', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
gbc = addspecies(cytosol, 'gbc', 'InitialAmount', 250, 'InitialAmountUnits', 'micromolarity');
mglur_tot = addspecies(cytosol, 'mglur_tot', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
mglur_tot_Rule = addrule(pfpc_model,'mglur_tot = mGLUR + mGLUR_gq + mGLUR_glu + mGLUR_gq_glu + gq_gdp + ga_gdp + ga_gtp + gbc','repeatedAssignment');

%Parameters
kon_mglur_glu = addparameter(pfpc_model,'kon_mglur_glu', 'Value', 11.11, 'ValueUnits', '1/(micromolarity*second)');
koff_mglur_glu = addparameter(pfpc_model,'koff_mglur_glu', 'Value', 100, 'ValueUnits', '1/second');
kon_mglur_gq = addparameter(pfpc_model,'kon_mglur_gq', 'Value', 2, 'ValueUnits', '1/(micromolarity*second)');
koff_mglur_gq = addparameter(pfpc_model,'koff_mglur_gq', 'Value', 100, 'ValueUnits', '1/second');
kact_gq = addparameter(pfpc_model,'kact_gq', 'Value', 116, 'ValueUnits', '1/second');
kact_gq_basal = addparameter(pfpc_model,'kact_gq_basal', 'Value', 0.0001, 'ValueUnits', '1/second');
kdeact_ga = addparameter(pfpc_model,'kdeact_ga', 'Value', 0.02, 'ValueUnits', '1/second');
ktrimer = addparameter(pfpc_model,'ktrimer', 'Value', 6, 'ValueUnits', '1/(micromolarity*second)');
kdeg_glu = addparameter(pfpc_model,'kdeg_glu', 'Value', 9000, 'ValueUnits', '1/second', 'ConstantValue', false);
kdeg_ca = addparameter(pfpc_model,'kdeg_ca', 'Value', 0, 'ValueUnits', '1/second', 'ConstantValue', false);
%Reactions
%A1
mglur_reac_1 = addreaction(pfpc_model, 'mGLUR + glu -> mGLUR_glu')
mglur_reac_1_k = addkineticlaw(mglur_reac_1,'MassAction');
set(mglur_reac_1_k, 'ParameterVariableNames', {'kon_mglur_glu'});

mglur_reac_2 = addreaction(pfpc_model, 'mGLUR_glu -> mGLUR + glu')
mglur_reac_2_k = addkineticlaw(mglur_reac_2,'MassAction');
set(mglur_reac_2_k, 'ParameterVariableNames', {'koff_mglur_glu'});
%A2
mglur_reac_3 = addreaction(pfpc_model, 'mGLUR_gq + glu -> mGLUR_gq_glu')
mglur_reac_3_k = addkineticlaw(mglur_reac_3,'MassAction');
set(mglur_reac_3_k, 'ParameterVariableNames', {'kon_mglur_glu'});

mglur_reac_4 = addreaction(pfpc_model, 'mGLUR_gq_glu -> mGLUR_gq + glu')
mglur_reac_4_k = addkineticlaw(mglur_reac_4,'MassAction');
set(mglur_reac_4_k, 'ParameterVariableNames', {'koff_mglur_glu'});

%A3
mglur_reac_5 = addreaction(pfpc_model, 'mGLUR + gq_gdp -> mGLUR_gq')
mglur_reac_5_k = addkineticlaw(mglur_reac_5,'MassAction');
set(mglur_reac_5_k, 'ParameterVariableNames', {'kon_mglur_gq'});

mglur_reac_6 = addreaction(pfpc_model, 'mGLUR_gq -> mGLUR + gq_gdp')
mglur_reac_6_k = addkineticlaw(mglur_reac_6,'MassAction');
set(mglur_reac_6_k, 'ParameterVariableNames', {'koff_mglur_gq'});

%A4

mglur_reac_7 = addreaction(pfpc_model, 'mGLUR_glu + gq_gdp -> mGLUR_gq_glu')
mglur_reac_7_k = addkineticlaw(mglur_reac_7,'MassAction');
set(mglur_reac_7_k, 'ParameterVariableNames', {'kon_mglur_gq'});

mglur_reac_8 = addreaction(pfpc_model, 'mGLUR_gq_glu -> mGLUR_glu + gq_gdp')
mglur_reac_8_k = addkineticlaw(mglur_reac_8,'MassAction');
set(mglur_reac_8_k, 'ParameterVariableNames', {'koff_mglur_gq'});

%A5
mglur_reac_9 = addreaction(pfpc_model, 'mGLUR_gq_glu -> mGLUR_glu + ga_gtp + gbc')
mglur_reac_9_k = addkineticlaw(mglur_reac_9,'MassAction');
set(mglur_reac_9_k, 'ParameterVariableNames', {'kact_gq'});


%A7
mglur_reac_11 = addreaction(pfpc_model, 'ga_gtp -> ga_gdp')
mglur_reac_11_k = addkineticlaw(mglur_reac_11,'MassAction');
set(mglur_reac_11_k, 'ParameterVariableNames', {'kdeact_ga'});

%A8
mglur_reac_12 = addreaction(pfpc_model, 'ga_gdp + gbc -> gq_gdp')
mglur_reac_12_k = addkineticlaw(mglur_reac_12,'MassAction');
set(mglur_reac_12_k, 'ParameterVariableNames', {'ktrimer'});


%A9 CLEAR GLUTAMATE
mglur_reac_13 = addreaction(pfpc_model, 'glu -> glu_deg')
mglur_reac_13_k = addkineticlaw(mglur_reac_13,'MassAction');
set(mglur_reac_13_k, 'ParameterVariableNames', {'kdeg_glu'});

%%%%%PLC Model%%%%%
%Species (uM)
IP3 = addspecies(cytosol, 'IP3', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
DAG = addspecies(cytosol, 'DAG', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
PLC = addspecies(cytosol, 'PLC', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
PLC_Ca = addspecies(cytosol, 'PLC_Ca', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
PLC_pip2 = addspecies(cytosol, 'PLC_pip2', 'InitialAmount', 35, 'InitialAmountUnits', 'micromolarity');
PLC_gq_Ca = addspecies(cytosol, 'PLC_gq_Ca', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
pip2 = addspecies(cytosol, 'pip2', 'InitialAmount', 4166.7, 'InitialAmountUnits', 'micromolarity');
PLC_pip2_Ca = addspecies(cytosol, 'PLC_pip2_Ca', 'InitialAmount', 6.25, 'InitialAmountUnits', 'micromolarity');
PLC_gq_pip2 = addspecies(cytosol, 'PLC_gq_pip2', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
PLC_gq_pip2_Ca = addspecies(cytosol, 'PLC_gq_pip2_Ca', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');

%Parameters
kon_plc_ca = addparameter(pfpc_model,'kon_plc_ca', 'Value', 300, 'ValueUnits', '1/(micromolarity*second)');
koff_plc_ca = addparameter(pfpc_model,'koff_plc_ca', 'Value', 100, 'ValueUnits', '1/second');
kon_plc_gq_ca = addparameter(pfpc_model,'kon_plc_gq_ca', 'Value', 900, 'ValueUnits', '1/(micromolarity*second)');
koff_plc_gq_ca = addparameter(pfpc_model,'koff_plc_gq_ca', 'Value', 30, 'ValueUnits', '1/second');
kon_plc_ga = addparameter(pfpc_model,'kon_plc_ga', 'Value', 800, 'ValueUnits', '1/(micromolarity*second)');
koff_plc_ga = addparameter(pfpc_model,'koff_plc_ga', 'Value', 40, 'ValueUnits', '1/second');
kon_plc_ca_ga = addparameter(pfpc_model,'kon_plc_ca_ga', 'Value', 1200, 'ValueUnits', '1/(micromolarity*second)');
koff_plc_ca_ga = addparameter(pfpc_model,'koff_plc_ca_ga', 'Value', 6, 'ValueUnits', '1/second');
kcat_plc_basal = addparameter(pfpc_model,'kcat_plc_basal', 'Value', 0, 'ValueUnits', '1/second');
kcat_plc_ga = addparameter(pfpc_model,'kcat_plc_ga', 'Value', 160, 'ValueUnits', '1/second');
kon_plc_ca_pip2 = addparameter(pfpc_model,'kon_plc_ca_pip2', 'Value', 1, 'ValueUnits', '1/(micromolarity*second)');
koff_plc_ca_pip2 = addparameter(pfpc_model,'koff_plc_ca_pip2', 'Value', 170, 'ValueUnits', '1/second');
kdeact_plc_ga = addparameter(pfpc_model,'kdeact_plc_ga', 'Value', 8, 'ValueUnits', '1/second');

%Reactions

%B1
plc_reac_1 = addreaction(pfpc_model, 'PLC_pip2 + Ca -> PLC_pip2_Ca')
plc_reac_1_k = addkineticlaw(plc_reac_1,'MassAction');
set(plc_reac_1_k, 'ParameterVariableNames', {'kon_plc_ca'});

plc_reac_2 = addreaction(pfpc_model, 'PLC_pip2_Ca -> PLC_pip2 + Ca')
plc_reac_2_k = addkineticlaw(plc_reac_2,'MassAction');
set(plc_reac_2_k, 'ParameterVariableNames', {'koff_plc_ca'});

%B2
plc_reac_3 = addreaction(pfpc_model, 'PLC_gq_pip2 + Ca -> PLC_gq_pip2_Ca')
plc_reac_3_k = addkineticlaw(plc_reac_3,'MassAction');
set(plc_reac_3_k, 'ParameterVariableNames', {'kon_plc_gq_ca'});

plc_reac_4 = addreaction(pfpc_model, 'PLC_gq_pip2_Ca -> PLC_gq_pip2 + Ca')
plc_reac_4_k = addkineticlaw(plc_reac_4,'MassAction');
set(plc_reac_4_k, 'ParameterVariableNames', {'koff_plc_gq_ca'});

%B3
plc_reac_5 = addreaction(pfpc_model, 'PLC_pip2 + ga_gtp -> PLC_gq_pip2')
plc_reac_5_k = addkineticlaw(plc_reac_5,'MassAction');
set(plc_reac_5_k, 'ParameterVariableNames', {'kon_plc_ga'});

plc_reac_6 = addreaction(pfpc_model, 'PLC_gq_pip2 -> PLC_pip2 + ga_gtp')
plc_reac_6_k = addkineticlaw(plc_reac_6,'MassAction');
set(plc_reac_6_k, 'ParameterVariableNames', {'koff_plc_ga'});

%B4
plc_reac_7 = addreaction(pfpc_model, 'PLC_pip2_Ca + ga_gtp -> PLC_gq_pip2_Ca')
plc_reac_7_k = addkineticlaw(plc_reac_7,'MassAction');
set(plc_reac_7_k, 'ParameterVariableNames', {'kon_plc_ca_ga'});

plc_reac_8 = addreaction(pfpc_model, 'PLC_gq_pip2_Ca -> PLC_pip2_Ca + ga_gtp')
plc_reac_8_k = addkineticlaw(plc_reac_8,'MassAction');
set(plc_reac_8_k, 'ParameterVariableNames', {'koff_plc_ca_ga'});

%B5
plc_reac_9 = addreaction(pfpc_model, 'PLC_Ca + ga_gtp -> PLC_gq_Ca')
plc_reac_9_k = addkineticlaw(plc_reac_9,'MassAction');
set(plc_reac_9_k, 'ParameterVariableNames', {'kon_plc_ca_ga'});

plc_reac_10 = addreaction(pfpc_model, 'PLC_gq_Ca -> PLC_Ca + ga_gtp')
plc_reac_10_k = addkineticlaw(plc_reac_10,'MassAction');
set(plc_reac_10_k, 'ParameterVariableNames', {'koff_plc_ca_ga'});

%B6
plc_reac_11 = addreaction(pfpc_model, 'PLC_pip2_Ca -> PLC_Ca + IP3 + DAG')
plc_reac_11_k = addkineticlaw(plc_reac_11,'MassAction');
set(plc_reac_11_k, 'ParameterVariableNames', {'kcat_plc_basal'});

%B7
plc_reac_12 = addreaction(pfpc_model, 'PLC_gq_pip2_Ca -> PLC_Ca + ga_gtp + IP3 + DAG + pip2')
plc_reac_12_k = addkineticlaw(plc_reac_12,'MassAction');
set(plc_reac_12_k, 'ParameterVariableNames', {'kcat_plc_ga'});

%B8
plc_reac_13 = addreaction(pfpc_model, 'PLC_Ca + pip2 -> PLC_pip2_Ca')
plc_reac_13_k = addkineticlaw(plc_reac_13,'MassAction');
set(plc_reac_13_k, 'ParameterVariableNames', {'kon_plc_ca_pip2'});

plc_reac_14 = addreaction(pfpc_model, 'PLC_pip2_Ca -> PLC_Ca + pip2')
plc_reac_14_k = addkineticlaw(plc_reac_14,'MassAction');
set(plc_reac_14_k, 'ParameterVariableNames', {'koff_plc_ca_pip2'});

%B9
plc_reac_15 = addreaction(pfpc_model, 'PLC_gq_Ca + pip2 -> PLC_gq_pip2_Ca')
plc_reac_15_k = addkineticlaw(plc_reac_15,'MassAction');
set(plc_reac_15_k, 'ParameterVariableNames', {'kon_plc_ca_pip2'});

plc_reac_16 = addreaction(pfpc_model, 'PLC_gq_pip2_Ca -> PLC_gq_Ca + pip2')
plc_reac_16_k = addkineticlaw(plc_reac_16,'MassAction');
set(plc_reac_16_k, 'ParameterVariableNames', {'koff_plc_ca_pip2'});

%B10
plc_reac_17 = addreaction(pfpc_model, 'PLC_gq_pip2 -> PLC_pip2 + ga_gdp')
plc_reac_17_k = addkineticlaw(plc_reac_17,'MassAction');
set(plc_reac_17_k, 'ParameterVariableNames', {'kdeact_plc_ga'});

%B11
plc_reac_18 = addreaction(pfpc_model, 'PLC_gq_pip2_Ca -> PLC_pip2_Ca + ga_gdp')
plc_reac_18_k = addkineticlaw(plc_reac_18,'MassAction');
set(plc_reac_18_k, 'ParameterVariableNames', {'kdeact_plc_ga'});

%B12
plc_reac_19 = addreaction(pfpc_model, 'PLC_gq_Ca -> PLC_Ca + ga_gdp')
plc_reac_19_k = addkineticlaw(plc_reac_19,'MassAction');
set(plc_reac_19_k, 'ParameterVariableNames', {'kdeact_plc_ga'});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% IP3R MODEL %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Species
IP3R = addspecies(cytosol, 'IP3R', 'InitialAmount', 0.337, 'InitialAmountUnits', 'micromolarity');
IP3R_IP3 = addspecies(cytosol, 'IP3R_IP3', 'InitialAmount', 0.1, 'InitialAmountUnits', 'micromolarity');
IP3R_open = addspecies(cytosol, 'IP3R_open', 'InitialAmount', 0.00016667, 'InitialAmountUnits', 'micromolarity');
IP3R_Ca = addspecies(cytosol, 'IP3R_Ca', 'InitialAmount', 0.025, 'InitialAmountUnits', 'micromolarity');
IP3R_2Ca = addspecies(cytosol, 'IP3R_2Ca', 'InitialAmount', 0.003, 'InitialAmountUnits', 'micromolarity');
IP3R_3Ca = addspecies(cytosol, 'IP3R_3Ca', 'InitialAmount', 0.0005, 'InitialAmountUnits', 'micromolarity');
IP3R_4Ca = addspecies(cytosol, 'IP3R_4Ca', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');

%Parameters
kon_ip3r_ip3 = addparameter(pfpc_model,'kon_ip3r_ip3', 'Value', 1000, 'ValueUnits', '1/(micromolarity*second)');
koff_ip3r_ip3 = addparameter(pfpc_model,'koff_ip3r_ip3', 'Value', 25800, 'ValueUnits', '1/second');
kopen_ip3r = addparameter(pfpc_model,'kopen_ip3r', 'Value', 8000, 'ValueUnits', '1/(micromolarity*second)');
kclose_ip3r = addparameter(pfpc_model,'kclose_ip3r', 'Value', 2000, 'ValueUnits', '1/second');
kon_ip3r_ca = addparameter(pfpc_model,'kon_ip3r_ca', 'Value', 8.889/4, 'ValueUnits', '1/(micromolarity*second)');
koff_ip3r_ca = addparameter(pfpc_model,'koff_ip3r_ca', 'Value', 5, 'ValueUnits', '1/second');
kon_ip3r_ca2 = addparameter(pfpc_model,'kon_ip3r_ca2', 'Value', 20/4, 'ValueUnits', '1/(micromolarity*second)');
koff_ip3r_ca2 = addparameter(pfpc_model,'koff_ip3r_ca2', 'Value', 10, 'ValueUnits', '1/second');
kon_ip3r_ca3 = addparameter(pfpc_model,'kon_ip3r_ca3', 'Value', 40/4, 'ValueUnits', '1/(micromolarity*second)');
koff_ip3r_ca3 = addparameter(pfpc_model,'koff_ip3r_ca3', 'Value', 15, 'ValueUnits', '1/second');
kon_ip3r_ca4 = addparameter(pfpc_model,'kon_ip3r_ca4', 'Value', 60/4, 'ValueUnits', '1/(micromolarity*second)');
koff_ip3r_ca4 = addparameter(pfpc_model,'koff_ip3r_ca4', 'Value', 20, 'ValueUnits', '1/second');

%Reactions
%D1
ip3r_reac_1 = addreaction(pfpc_model, 'IP3R + IP3 -> IP3R_IP3')
ip3r_reac_1_k = addkineticlaw(ip3r_reac_1,'MassAction');
set(ip3r_reac_1_k, 'ParameterVariableNames', {'kon_ip3r_ip3'});
ip3r_reac_2 = addreaction(pfpc_model, 'IP3R_IP3 -> IP3R + IP3')
ip3r_reac_2_k = addkineticlaw(ip3r_reac_2,'MassAction');
set(ip3r_reac_2_k, 'ParameterVariableNames', {'koff_ip3r_ip3'});
%D2
ip3r_reac_3 = addreaction(pfpc_model, 'IP3R_IP3 + Ca -> IP3R_open')
ip3r_reac_3_k = addkineticlaw(ip3r_reac_3,'MassAction');
set(ip3r_reac_3_k, 'ParameterVariableNames', {'kopen_ip3r'});
ip3r_reac_4 = addreaction(pfpc_model, 'IP3R_open -> IP3R_IP3 + Ca')
ip3r_reac_4_k = addkineticlaw(ip3r_reac_4,'MassAction');
set(ip3r_reac_4_k, 'ParameterVariableNames', {'kclose_ip3r'});
%D3
ip3r_reac_5 = addreaction(pfpc_model, 'IP3R + Ca -> IP3R_Ca')
ip3r_reac_5_k = addkineticlaw(ip3r_reac_5,'MassAction');
set(ip3r_reac_5_k, 'ParameterVariableNames', {'kon_ip3r_ca'});
ip3r_reac_6 = addreaction(pfpc_model, 'IP3R_Ca -> IP3R + Ca')
ip3r_reac_6_k = addkineticlaw(ip3r_reac_6,'MassAction');
set(ip3r_reac_6_k, 'ParameterVariableNames', {'koff_ip3r_ca'});
%D4
ip3r_reac_7 = addreaction(pfpc_model, 'IP3R_Ca + Ca -> IP3R_2Ca')
ip3r_reac_7_k = addkineticlaw(ip3r_reac_7,'MassAction');
set(ip3r_reac_7_k, 'ParameterVariableNames', {'kon_ip3r_ca2'});
ip3r_reac_8 = addreaction(pfpc_model, 'IP3R_2Ca -> IP3R_Ca + Ca')
ip3r_reac_8_k = addkineticlaw(ip3r_reac_8,'MassAction');
set(ip3r_reac_8_k, 'ParameterVariableNames', {'koff_ip3r_ca2'});
%D5
ip3r_reac_9 = addreaction(pfpc_model, 'IP3R_2Ca + Ca -> IP3R_3Ca')
ip3r_reac_9_k = addkineticlaw(ip3r_reac_9,'MassAction');
set(ip3r_reac_9_k, 'ParameterVariableNames', {'kon_ip3r_ca3'});
ip3r_reac_10 = addreaction(pfpc_model, 'IP3R_3Ca -> IP3R_2Ca + Ca')
ip3r_reac_10_k = addkineticlaw(ip3r_reac_10,'MassAction');
set(ip3r_reac_10_k, 'ParameterVariableNames', {'koff_ip3r_ca3'});
%D6
ip3r_reac_11 = addreaction(pfpc_model, 'IP3R_3Ca + Ca -> IP3R_4Ca')
ip3r_reac_11_k = addkineticlaw(ip3r_reac_11,'MassAction');
set(ip3r_reac_11_k, 'ParameterVariableNames', {'kon_ip3r_ca4'});
ip3r_reac_12 = addreaction(pfpc_model, 'IP3R_4Ca -> IP3R_3Ca + Ca')
ip3r_reac_12_k = addkineticlaw(ip3r_reac_12,'MassAction');
set(ip3r_reac_12_k, 'ParameterVariableNames', {'koff_ip3r_ca4'});

%%%%%Calcium Regulation Model%%%%%
%Species (uM)
Ca_EXT = addspecies(cytosol, 'Ca_EXT', 'InitialAmount', 2000, 'InitialAmountUnits', 'micromolarity');


%%%%%%IP3 DEGRADATION MODEL%%%%%%
%Species
IP3K = addspecies(cytosol, 'IP3K', 'InitialAmount', 0.86667, 'InitialAmountUnits', 'micromolarity');
IP3K_Ca = addspecies(cytosol, 'IP3K_Ca', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
IP3K_2Ca = addspecies(cytosol, 'IP3K_2Ca', 'InitialAmount', 0.033333, 'InitialAmountUnits', 'micromolarity');
IP3K_2Ca_IP3 = addspecies(cytosol, 'IP3K_2Ca_IP3', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
IP2 = addspecies(cytosol, 'IP2', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
IP4 = addspecies(cytosol, 'IP4', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
IP3P = addspecies(cytosol, 'IP3P', 'InitialAmount', 0.98, 'InitialAmountUnits', 'micromolarity');
IP3P_IP3 = addspecies(cytosol, 'IP3P_IP3', 'InitialAmount', 0.02, 'InitialAmountUnits', 'micromolarity');

%Parameters
kon_ip3k_ca = addparameter(pfpc_model,'kon_ip3k_ca', 'Value', 1111.1, 'ValueUnits', '1/(micromolarity*second)');
koff_ip3k_ca = addparameter(pfpc_model,'koff_ip3k_ca', 'Value', 100, 'ValueUnits', '1/second');
kon_ip3k_ip3 = addparameter(pfpc_model,'kon_ip3k_ip3', 'Value', 100, 'ValueUnits', '1/(micromolarity*second)');
koff_ip3k_ip3 = addparameter(pfpc_model,'koff_ip3k_ip3', 'Value', 80, 'ValueUnits', '1/second');
kcat_ip3k = addparameter(pfpc_model,'kcat_ip3k', 'Value', 3, 'ValueUnits', '1/second');
kon_ip3p_ip3 = addparameter(pfpc_model,'kon_ip3p_ip3', 'Value', 9, 'ValueUnits', '1/(micromolarity*second)');
koff_ip3p_ip3 = addparameter(pfpc_model,'koff_ip3p_ip3', 'Value', 72, 'ValueUnits', '1/second');
kcat_ip3p = addparameter(pfpc_model,'kcat_ip3p', 'Value', 10, 'ValueUnits', '1/second');

%Reactions

%C1
ip3_deg_reac_1 = addreaction(pfpc_model, 'IP3K + Ca -> IP3K_Ca')
ip3_deg_reac_1_k = addkineticlaw(ip3_deg_reac_1,'MassAction');
set(ip3_deg_reac_1_k, 'ParameterVariableNames', {'kon_ip3k_ca'});

ip3_deg_reac_2 = addreaction(pfpc_model, 'IP3K_Ca -> IP3K + Ca')
ip3_deg_reac_2_k = addkineticlaw(ip3_deg_reac_2,'MassAction');
set(ip3_deg_reac_2_k, 'ParameterVariableNames', {'koff_ip3k_ca'});

ip3_deg_reac_3 = addreaction(pfpc_model, 'IP3K_Ca + Ca -> IP3K_2Ca')
ip3_deg_reac_3_k = addkineticlaw(ip3_deg_reac_3,'MassAction');
set(ip3_deg_reac_3_k, 'ParameterVariableNames', {'kon_ip3k_ca'});

ip3_deg_reac_4 = addreaction(pfpc_model, 'IP3K_2Ca -> IP3K_Ca + Ca')
ip3_deg_reac_4_k = addkineticlaw(ip3_deg_reac_4,'MassAction');
set(ip3_deg_reac_4_k, 'ParameterVariableNames', {'koff_ip3k_ca'});

%C2
ip3_deg_reac_5 = addreaction(pfpc_model, 'IP3K_2Ca + IP3 -> IP3K_2Ca_IP3')
ip3_deg_reac_5_k = addkineticlaw(ip3_deg_reac_5,'MassAction');
set(ip3_deg_reac_5_k, 'ParameterVariableNames', {'kon_ip3k_ip3'});

ip3_deg_reac_6 = addreaction(pfpc_model, 'IP3K_2Ca_IP3 -> IP3K_2Ca + IP3')
ip3_deg_reac_6_k = addkineticlaw(ip3_deg_reac_6,'MassAction');
set(ip3_deg_reac_6_k, 'ParameterVariableNames', {'koff_ip3k_ip3'});


%C3
ip3_deg_reac_7 = addreaction(pfpc_model, 'IP3K_2Ca_IP3 -> IP3K_2Ca + IP4')
ip3_deg_reac_7_k = addkineticlaw(ip3_deg_reac_7,'MassAction');
set(ip3_deg_reac_7_k, 'ParameterVariableNames', {'kcat_ip3k'});

%C4
ip3_deg_reac_8 = addreaction(pfpc_model, 'IP3P + IP3 -> IP3P_IP3')
ip3_deg_reac_8_k = addkineticlaw(ip3_deg_reac_8,'MassAction');
set(ip3_deg_reac_8_k, 'ParameterVariableNames', {'kon_ip3p_ip3'});

ip3_deg_reac_9 = addreaction(pfpc_model, 'IP3P_IP3 -> IP3P + IP3')
ip3_deg_reac_9_k = addkineticlaw(ip3_deg_reac_9,'MassAction');
set(ip3_deg_reac_9_k, 'ParameterVariableNames', {'koff_ip3p_ip3'});

%C5
ip3_deg_reac_10 = addreaction(pfpc_model, 'IP3P_IP3 -> IP3P + IP2')
ip3_deg_reac_10_k = addkineticlaw(ip3_deg_reac_10,'MassAction');
set(ip3_deg_reac_10_k, 'ParameterVariableNames', {'kcat_ip3p'});


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%CAMKII MODEL%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Add species
PP2A = addspecies(cytosol, 'PP2A', 'InitialAmount', 0.5, 'InitialAmountUnits', 'micromolarity');
PP1 = addspecies(cytosol, 'PP1', 'InitialAmount', 0.5, 'InitialAmountUnits', 'micromolarity');

CaM = addspecies(cytosol, 'CaM', 'InitialAmount', 60, 'InitialAmountUnits', 'micromolarity');
Ca2CaM = addspecies(cytosol, 'Ca2CaM', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
Ca3CaM = addspecies(cytosol, 'Ca3CaM', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
Ca4CaM = addspecies(cytosol, 'Ca4CaM', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');

CaMKII_00_00 = addspecies(cytosol, 'CaMKII_00_00', 'InitialAmount', 20, 'InitialAmountUnits', 'micromolarity');
CaMKII_00_10 = addspecies(cytosol, 'CaMKII_00_10', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
CaMKII_00_11 = addspecies(cytosol, 'CaMKII_00_11', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
CaMKII_10_00 = addspecies(cytosol, 'CaMKII_10_00', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
CaMKII_10_00_PP2A = addspecies(cytosol, 'CaMKII_10_00_PP2A', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
CaMKII_11_00 = addspecies(cytosol, 'CaMKII_11_00', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
CaMKII_11_00_PP2A = addspecies(cytosol, 'CaMKII_11_00_PP2A', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
CaMKII_10a_00 = addspecies(cytosol, 'CaMKII_10a_00', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
CaMKII_12_00 = addspecies(cytosol, 'CaMKII_12_00', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
CaMKII_12_00_PP2A = addspecies(cytosol, 'CaMKII_12_00_PP2A', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
CaMKII_20a_00 = addspecies(cytosol, 'CaMKII_20a_00', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
CaMKII_22_00 = addspecies(cytosol, 'CaMKII_22_00', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
CaMKII_22_00_PP2A = addspecies(cytosol, 'CaMKII_22_00_PP2A', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
CaMKII_20_00 = addspecies(cytosol, 'CaMKII_20_00', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
CaMKII_20_00_PP2A = addspecies(cytosol, 'CaMKII_20_00_PP2A', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
CaMKII_10_10 = addspecies(cytosol, 'CaMKII_10_10', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
CaMKII_10_10_PP2A = addspecies(cytosol, 'CaMKII_10_10_PP2A', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
CaMKII_10_01 = addspecies(cytosol, 'CaMKII_10_01', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
CaMKII_10_01_PP2A = addspecies(cytosol, 'CaMKII_10_01_PP2A', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
CaMKII_11_10 = addspecies(cytosol, 'CaMKII_11_10', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
CaMKII_11_10_PP2A = addspecies(cytosol, 'CaMKII_11_10_PP2A', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
CaMKII_10a_10 = addspecies(cytosol, 'CaMKII_10a_10', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
CaMKII_20_01 = addspecies(cytosol, 'CaMKII_20_01', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
CaMKII_20_01_PP2A = addspecies(cytosol, 'CaMKII_20_01_PP2A', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
CaMKII_12_10 = addspecies(cytosol, 'CaMKII_12_10', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
CaMKII_12_10_PP2A = addspecies(cytosol, 'CaMKII_12_10_PP2A', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
CaMKII_20a_01 = addspecies(cytosol, 'CaMKII_20a_01', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
CaMKII_10_11 = addspecies(cytosol, 'CaMKII_10_11', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
CaMKII_10_11_PP2A = addspecies(cytosol, 'CaMKII_10_11_PP2A', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
CaMKII_11_11 = addspecies(cytosol, 'CaMKII_11_11', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
CaMKII_11_11_PP2A = addspecies(cytosol, 'CaMKII_11_11_PP2A', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');

%%%%%PDE1 Activation and Regulation%%%%%
PDE1 = addspecies(cytosol, 'PDE1', 'InitialAmount', 0.6, 'InitialAmountUnits', 'micromolarity');
PDE5 = addspecies(cytosol, 'PDE5', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
PDE1_CaM = addspecies(cytosol, 'PDE1_CaM', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
PDE1p = addspecies(cytosol, 'PDE1p', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
PDE1p_CaM = addspecies(cytosol, 'PDE1p_CaM', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');

PDE_active = addspecies(cytosol, 'PDE_active', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
PDE_Active_Rule = addrule(pfpc_model,'PDE_active = PDE1_CaM + PDE1p_CaM','repeatedAssignment');

%%phosphorylation of PDE1 by CaMKII
CaMKII_00_10_PDE1 = addspecies(cytosol, 'CaMKII_00_10_PDE1', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
CaMKII_00_11_PDE1 = addspecies(cytosol, 'CaMKII_00_11_PDE1', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
CaMKII_10_00_PDE1 = addspecies(cytosol, 'CaMKII_10_00_PDE1', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
CaMKII_20_00_PDE1 = addspecies(cytosol, 'CaMKII_20_00_PDE1', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
CaMKII_10_10_PDE1 = addspecies(cytosol, 'CaMKII_10_10_PDE1', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
CaMKII_10_01_PDE1 = addspecies(cytosol, 'CaMKII_10_01_PDE1', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
CaMKII_20_01_PDE1 = addspecies(cytosol, 'CaMKII_20_01_PDE1', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
CaMKII_10_11_PDE1 = addspecies(cytosol, 'CaMKII_10_11_PDE1', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
CaMKII_11_11_PDE1 = addspecies(cytosol, 'CaMKII_11_11_PDE1', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
CaMKII_11_00_PDE1 = addspecies(cytosol, 'CaMKII_11_00_PDE1', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
CaMKII_12_00_PDE1 = addspecies(cytosol, 'CaMKII_12_00_PDE1', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
CaMKII_22_00_PDE1 = addspecies(cytosol, 'CaMKII_22_00_PDE1', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
CaMKII_11_10_PDE1 = addspecies(cytosol, 'CaMKII_11_10_PDE1', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
CaMKII_12_10_PDE1 = addspecies(cytosol, 'CaMKII_12_10_PDE1', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
CaMKII_00_10_PDE1_CaM = addspecies(cytosol, 'CaMKII_00_10_PDE1_CaM', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
CaMKII_00_11_PDE1_CaM = addspecies(cytosol, 'CaMKII_00_11_PDE1_CaM', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
CaMKII_10_00_PDE1_CaM = addspecies(cytosol, 'CaMKII_10_00_PDE1_CaM', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
CaMKII_20_00_PDE1_CaM = addspecies(cytosol, 'CaMKII_20_00_PDE1_CaM', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
CaMKII_10_10_PDE1_CaM = addspecies(cytosol, 'CaMKII_10_10_PDE1_CaM', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
CaMKII_10_01_PDE1_CaM = addspecies(cytosol, 'CaMKII_10_01_PDE1_CaM', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
CaMKII_20_01_PDE1_CaM = addspecies(cytosol, 'CaMKII_20_01_PDE1_CaM', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
CaMKII_10_11_PDE1_CaM = addspecies(cytosol, 'CaMKII_10_11_PDE1_CaM', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
CaMKII_11_11_PDE1_CaM = addspecies(cytosol, 'CaMKII_11_11_PDE1_CaM', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
CaMKII_11_00_PDE1_CaM = addspecies(cytosol, 'CaMKII_11_00_PDE1_CaM', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
CaMKII_12_00_PDE1_CaM = addspecies(cytosol, 'CaMKII_12_00_PDE1_CaM', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
CaMKII_22_00_PDE1_CaM = addspecies(cytosol, 'CaMKII_22_00_PDE1_CaM', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
CaMKII_11_10_PDE1_CaM = addspecies(cytosol, 'CaMKII_11_10_PDE1_CaM', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
CaMKII_12_10_PDE1_CaM = addspecies(cytosol, 'CaMKII_12_10_PDE1_CaM', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');

%%%%cGMP Production by GC%%%%
GC = addspecies(cytosol, 'GC', 'InitialAmount', 0.04, 'InitialAmountUnits', 'micromolarity');
NO = addspecies(cytosol, 'NO', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity', 'ConstantAmount', false);
GC_NO = addspecies(cytosol, 'GC_NO', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');

GTP = addspecies(cytosol, 'GTP', 'InitialAmount', 100, 'InitialAmountUnits', 'micromolarity');
GC_NO_GTP = addspecies(cytosol, 'GC_NO_GTP', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
cGMP = addspecies(cytosol, 'cGMP', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');

%%%%Hydrolysis of cGMP by PDE1%%%%
%Basal hydrolysis
PDE1_cGMP = addspecies(cytosol, 'PDE1_cGMP', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
PDE1p_cGMP = addspecies(cytosol, 'PDE1p_cGMP', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');

%Active hydrolysis
PDE1_CaM_cGMP = addspecies(cytosol, 'PDE1_CaM_cGMP', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
PDE1p_CaM_cGMP = addspecies(cytosol, 'PDE1p_CaM_cGMP', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');

PDE5_cGMP = addspecies(cytosol, 'PDE5_cGMP', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');

%%%%Activation of PKG by cGMP%%%%
PKG = addspecies(cytosol, 'PKG', 'InitialAmount', 0.2, 'InitialAmountUnits', 'micromolarity');
PKG_cGMP = addspecies(cytosol, 'PKG_cGMP', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');

%%%%Phosphorylation of Gsubstrate by PKG
Gsub = addspecies(cytosol, 'Gsub', 'InitialAmount', 3.5, 'InitialAmountUnits', 'micromolarity');
PKG_cGMP_Gsub = addspecies(cytosol, 'PKG_cGMP_Gsub', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
pGsub = addspecies(cytosol, 'pGsub', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');

%%%%Inhibition of PP2A by pGsub%%%%
PP2A_pGsub = addspecies(cytosol, 'PP2A_pGsub', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');

%%%%Inhibition of PP1 by pGsub%%%%
PP1_pGsub = addspecies(cytosol, 'PP1_pGsub', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');

%%%%Dephosphorylation of Gsubstrate by PP2A
pGsub_PP2A = addspecies(cytosol, 'pGsub_PP2A', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');

%%%%Calcineurin binding to CaM
CaN = addspecies(cytosol, 'CaN', 'InitialAmount', 1, 'InitialAmountUnits', 'micromolarity');
CaN_CaM = addspecies(cytosol, 'CaN_CaM', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');

%%%%Dephosphorylation of Gsubstrate by CaN
CaN_pGsub = addspecies(cytosol, 'CaN_pGsub', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');

%%%%Autonomous CaMKII Percentage%%%%%
CaMKII_Auton = addspecies(cytosol, 'CaMKII_Auton', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');

%Each CaMKII has two subunits, each with a single T286 site. So this percentage is of the total number of phosphorylated sites (i.e. 2*[CaMKII]_t0)
CaMKII_Auton_Rule = addrule(pfpc_model, 'CaMKII_Auton = 100*((CaMKII_10_00 + CaMKII_20_00 + CaMKII_10_10 + CaMKII_10_01 + CaMKII_20_01  + CaMKII_10_11 + CaMKII_10a_00 + CaMKII_20a_00 + CaMKII_10a_10 + CaMKII_20a_01 + CaMKII_10_10_PP2A + CaMKII_10_01_PP2A + CaMKII_20_00_PP2A + CaMKII_20_01_PP2A + CaMKII_10_11_PP2A + CaMKII_10_10_PP1 + CaMKII_10_01_PP1 + CaMKII_20_00_PP1 + CaMKII_20_01_PP1 + CaMKII_10_11_PP1 + (2*(CaMKII_11_11 + CaMKII_11_00 + CaMKII_12_00 + CaMKII_22_00 + CaMKII_11_10 + CaMKII_12_10 + CaMKII_11_00_PP2A + CaMKII_12_00_PP2A + CaMKII_22_00_PP2A + CaMKII_11_10_PP2A + CaMKII_12_10_PP2A + CaMKII_11_11_PP2A + CaMKII_11_00_PP1 + CaMKII_12_00_PP1 + CaMKII_22_00_PP1 + CaMKII_11_10_PP1 + CaMKII_12_10_PP1 + CaMKII_11_11_PP1)))/40)', 'repeatedAssignment');

CaMKII_305p = addspecies(cytosol, 'CaMKII_305p', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');

CaMKII_305p_Rule = addrule(pfpc_model, 'CaMKII_305p = CaMKII_20a_00 + CaMKII_22_00 + CaMKII_22_00_PP2A + CaMKII_20_00 + CaMKII_20_00_PP2A + CaMKII_20_01 + CaMKII_20_01_PP2A + CaMKII_20a_01 + CaMKII_20_00_PDE1 + CaMKII_20_01_PDE1 + CaMKII_22_00_PDE1 + CaMKII_20_00_PDE1_CaM + CaMKII_20_01_PDE1_CaM + CaMKII_22_00_PDE1_CaM', 'repeatedAssignment');

%%%%Active PP2A%%%%%
PP2A_active = addspecies(cytosol, 'PP2A_active', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');

PP2A_active_rule = addrule(pfpc_model, 'PP2A_active = PP2A+CaMKII_10_00_PP2A+CaMKII_11_00_PP2A+CaMKII_12_00_PP2A+CaMKII_22_00_PP2A+CaMKII_20_00_PP2A+CaMKII_10_10_PP2A+CaMKII_10_01_PP2A+CaMKII_11_10_PP2A+CaMKII_20_01_PP2A+CaMKII_12_10_PP2A+CaMKII_10_11_PP2A+CaMKII_11_11_PP2A+pGsub_PP2A', 'repeatedAssignment');


%%%%%%CaM Activation model parameters%%%%%%
kon_cam_ca1 = addparameter(pfpc_model,'kon_cam_ca1', 'Value', 16.7788, 'ValueUnits', '1/(micromolarity*micromolarity*second)');
koff_cam_ca1 = addparameter(pfpc_model,'koff_cam_ca1', 'Value', 34.9426, 'ValueUnits', '1/second');

kon_cam_ca3 = addparameter(pfpc_model,'kon_cam_ca3', 'Value', 13.5158, 'ValueUnits', '1/(micromolarity*second)');
koff_cam_ca3 = addparameter(pfpc_model,'koff_cam_ca3', 'Value', 228.3048, 'ValueUnits', '1/second');

kon_cam_ca4 = addparameter(pfpc_model,'kon_cam_ca4', 'Value', 26.2878, 'ValueUnits', '1/(micromolarity*second)');
koff_cam_ca4 = addparameter(pfpc_model,'koff_cam_ca4', 'Value', 64.0459, 'ValueUnits', '1/second');


%%%%%%%%%%CaM Activation Model Reactions%%%%%%%%%%%%%

%CaM Reactions
cam_reac_1 = addreaction(pfpc_model, 'CaM + Ca + Ca -> Ca2CaM');
cam_reac_1_k = addkineticlaw(cam_reac_1,'MassAction');
set(cam_reac_1_k, 'ParameterVariableNames', {'kon_cam_ca1'});

cam_reac_2 = addreaction(pfpc_model, 'Ca2CaM -> CaM + Ca + Ca');
cam_reac_2_k = addkineticlaw(cam_reac_2,'MassAction');
set(cam_reac_2_k, 'ParameterVariableNames', {'koff_cam_ca1'});

cam_reac_3 = addreaction(pfpc_model, 'Ca2CaM + Ca -> Ca3CaM');
cam_reac_3_k = addkineticlaw(cam_reac_3,'MassAction');
set(cam_reac_3_k, 'ParameterVariableNames', {'kon_cam_ca3'});

cam_reac_4 = addreaction(pfpc_model, 'Ca3CaM -> Ca2CaM + Ca');
cam_reac_4_k = addkineticlaw(cam_reac_4,'MassAction');
set(cam_reac_4_k, 'ParameterVariableNames', {'koff_cam_ca3'});

cam_reac_5 = addreaction(pfpc_model, 'Ca3CaM + Ca -> Ca4CaM');
cam_reac_5_k = addkineticlaw(cam_reac_5,'MassAction');
set(cam_reac_5_k, 'ParameterVariableNames', {'kon_cam_ca4'});

cam_reac_6 = addreaction(pfpc_model, 'Ca4CaM -> Ca3CaM + Ca');
cam_reac_6_k = addkineticlaw(cam_reac_6,'MassAction');
set(cam_reac_6_k, 'ParameterVariableNames', {'koff_cam_ca4'});

%%%%%CaMKII Activation Model Parameters%%%%
%binding of CaM to CamKII
k1 = addparameter(pfpc_model,'k1', 'Value', 0.0590, 'ValueUnits', '1/(micromolarity*second)');
k1_2 = addparameter(pfpc_model,'k1_2', 'Value', 1.3150, 'ValueUnits', '1/(micromolarity*second)');
%unbinding of CaM to CamKII
k1r = addparameter(pfpc_model,'k1r', 'Value', 2.848, 'ValueUnits', '1/second');
k1r_2 = addparameter(pfpc_model,'k1r_2', 'Value', 5.696, 'ValueUnits', '1/second');
%autophosphorylation at T286
k2 = addparameter(pfpc_model,'k2', 'Value', 4.5176, 'ValueUnits', '1/second');
k2_2 = addparameter(pfpc_model,'k2_2', 'Value', 16.3875, 'ValueUnits', '1/second');
%binding of CaM to CamKIIpT286
k10 = addparameter(pfpc_model,'k10', 'Value', 38.4349, 'ValueUnits', '1/(micromolarity*second)');
k10_2 = addparameter(pfpc_model,'k10_2', 'Value', 8.4000, 'ValueUnits', '1/(micromolarity*second)');
%unbinding of CaM to CamKIIpT286
k10r = addparameter(pfpc_model,'k10r', 'Value', 0.33, 'ValueUnits', '1/second');
k10r_2 = addparameter(pfpc_model,'k10r_2', 'Value', 0.62, 'ValueUnits', '1/second');
%autophosphorylation at T305
k3 = addparameter(pfpc_model,'k3', 'Value', 0.1, 'ValueUnits', '1/second');
k3_2 = addparameter(pfpc_model,'k3_2', 'Value', 0.2, 'ValueUnits', '1/second');
%Deactivation of CaMKII after dephosphorylation
kdeact_cam = addparameter(pfpc_model,'kdeact_cam', 'Value', 2, 'ValueUnits', '1/second');
%dephosporylation by PP2A
kon_camkii_pp2a = addparameter(pfpc_model,'kon_camkii_pp2a', 'Value', 2, 'ValueUnits', '1/(micromolarity*second)', 'ConstantValue', false);
kon_camkii_pp2a_2 = addparameter(pfpc_model,'kon_camkii_pp2a_2', 'Value', 4, 'ValueUnits', '1/(micromolarity*second)', 'ConstantValue', false);
koff_camkii_pp2a = addparameter(pfpc_model,'koff_camkii_pp2a', 'Value', 0.8, 'ValueUnits', '1/second');
kcat_pp2a_cam = addparameter(pfpc_model,'kcat_pp2a_cam', 'Value', 0.2, 'ValueUnits', '1/second');

%dephosporylation by PP1
kon_camkii_pp1 = addparameter(pfpc_model,'kon_camkii_pp1', 'Value', 2, 'ValueUnits', '1/(micromolarity*second)', 'ConstantValue', false);
kon_camkii_pp1_2 = addparameter(pfpc_model,'kon_camkii_pp1_2', 'Value', 4, 'ValueUnits', '1/(micromolarity*second)', 'ConstantValue', false);
koff_camkii_pp1 = addparameter(pfpc_model,'koff_camkii_pp1', 'Value', 0.8, 'ValueUnits', '1/second');
kcat_pp1_cam = addparameter(pfpc_model,'kcat_pp1_cam', 'Value', 0.2, 'ValueUnits', '1/second');

%%%%%CaMKII Activation Model Reactions%%%%
%1
camkii_reac_1 = addreaction(pfpc_model, 'CaMKII_00_00 + Ca4CaM -> CaMKII_00_10');
camkii_reac_1_k = addkineticlaw(camkii_reac_1,'MassAction');
set(camkii_reac_1_k, 'ParameterVariableNames', {'k1_2'});

%2
camkii_reac_2a = addreaction(pfpc_model, 'CaMKII_00_10 -> CaMKII_00_00 + Ca4CaM');
camkii_reac_2a_k = addkineticlaw(camkii_reac_2a,'MassAction');
set(camkii_reac_2a_k, 'ParameterVariableNames', {'k1r'});

camkii_reac_2b = addreaction(pfpc_model, 'CaMKII_00_10 + Ca4CaM -> CaMKII_00_11');
camkii_reac_2b_k = addkineticlaw(camkii_reac_2b,'MassAction');
set(camkii_reac_2b_k, 'ParameterVariableNames', {'k1'});

%3
camkii_reac_3a = addreaction(pfpc_model, 'CaMKII_00_11 -> CaMKII_00_10 + Ca4CaM');
camkii_reac_3a_k = addkineticlaw(camkii_reac_3a,'MassAction');
set(camkii_reac_3a_k, 'ParameterVariableNames', {'k1r_2'});

camkii_reac_3b = addreaction(pfpc_model, 'CaMKII_00_11 -> CaMKII_10_11');
camkii_reac_3b_k = addkineticlaw(camkii_reac_3b,'MassAction');
set(camkii_reac_3b_k, 'ParameterVariableNames', {'k2_2'});

%4
camkii_reac_4a = addreaction(pfpc_model, 'CaMKII_10_00 + Ca4CaM -> CaMKII_10_10');
camkii_reac_4a_k = addkineticlaw(camkii_reac_4a,'MassAction');
set(camkii_reac_4a_k, 'ParameterVariableNames', {'k10'});

camkii_reac_4b = addreaction(pfpc_model, 'CaMKII_10_00 + Ca4CaM -> CaMKII_10_01');
camkii_reac_4b_k = addkineticlaw(camkii_reac_4b,'MassAction');
set(camkii_reac_4b_k, 'ParameterVariableNames', {'k1'});

camkii_reac_4c = addreaction(pfpc_model, 'CaMKII_10_00 + PP2A -> CaMKII_10_00_PP2A');
camkii_reac_4c_k = addkineticlaw(camkii_reac_4c,'MassAction');
set(camkii_reac_4c_k, 'ParameterVariableNames', {'kon_camkii_pp2a'});

camkii_reac_4d = addreaction(pfpc_model, 'CaMKII_10_00_PP2A -> CaMKII_10_00 + PP2A');
camkii_reac_4d_k = addkineticlaw(camkii_reac_4d,'MassAction');
set(camkii_reac_4d_k, 'ParameterVariableNames', {'koff_camkii_pp2a'});

camkii_reac_4e = addreaction(pfpc_model, 'CaMKII_10_00_PP2A -> CaMKII_00_00 + PP2A');
camkii_reac_4e_k = addkineticlaw(camkii_reac_4e,'MassAction');
set(camkii_reac_4e_k, 'ParameterVariableNames', {'kcat_pp2a_cam'});

camkii_reac_4f = addreaction(pfpc_model, 'CaMKII_10_00 -> CaMKII_20_00');
camkii_reac_4f_k = addkineticlaw(camkii_reac_4f,'MassAction');
set(camkii_reac_4f_k, 'ParameterVariableNames', {'k3'});

%5
camkii_reac_5a = addreaction(pfpc_model, 'CaMKII_11_00 + Ca4CaM -> CaMKII_11_10');
camkii_reac_5a_k = addkineticlaw(camkii_reac_5a,'MassAction');
set(camkii_reac_5a_k, 'ParameterVariableNames', {'k10_2'});

camkii_reac_5b = addreaction(pfpc_model, 'CaMKII_11_00 + PP2A -> CaMKII_11_00_PP2A');
camkii_reac_5b_k = addkineticlaw(camkii_reac_5b,'MassAction');
set(camkii_reac_5b_k, 'ParameterVariableNames', {'kon_camkii_pp2a_2'});

camkii_reac_5c = addreaction(pfpc_model, 'CaMKII_11_00_PP2A -> CaMKII_11_00 + PP2A');
camkii_reac_5c_k = addkineticlaw(camkii_reac_5c,'MassAction');
set(camkii_reac_5c_k, 'ParameterVariableNames', {'koff_camkii_pp2a'});

camkii_reac_5d = addreaction(pfpc_model, 'CaMKII_11_00_PP2A -> CaMKII_10a_00 + PP2A');
camkii_reac_5d_k = addkineticlaw(camkii_reac_5d,'MassAction');
set(camkii_reac_5d_k, 'ParameterVariableNames', {'kcat_pp2a_cam'});

camkii_reac_5e = addreaction(pfpc_model, 'CaMKII_10a_00 -> CaMKII_11_00');
camkii_reac_5e_k = addkineticlaw(camkii_reac_5e,'MassAction');
set(camkii_reac_5e_k, 'ParameterVariableNames', {'k2'});

camkii_reac_5f = addreaction(pfpc_model, 'CaMKII_10a_00 -> CaMKII_10_00');
camkii_reac_5f_k = addkineticlaw(camkii_reac_5f,'MassAction');
set(camkii_reac_5f_k, 'ParameterVariableNames', {'kdeact_cam'});

camkii_reac_5g = addreaction(pfpc_model, 'CaMKII_11_00 -> CaMKII_12_00');
camkii_reac_5g_k = addkineticlaw(camkii_reac_5g,'MassAction');
set(camkii_reac_5g_k, 'ParameterVariableNames', {'k3_2'});

%6
camkii_reac_6a = addreaction(pfpc_model, 'CaMKII_12_00 + Ca4CaM -> CaMKII_12_10');
camkii_reac_6a_k = addkineticlaw(camkii_reac_6a,'MassAction');
set(camkii_reac_6a_k, 'ParameterVariableNames', {'k10'});

camkii_reac_6b = addreaction(pfpc_model, 'CaMKII_12_00 + PP2A -> CaMKII_12_00_PP2A');
camkii_reac_6b_k = addkineticlaw(camkii_reac_6b,'MassAction');
set(camkii_reac_6b_k, 'ParameterVariableNames', {'kon_camkii_pp2a'});

camkii_reac_6c = addreaction(pfpc_model, 'CaMKII_12_00_PP2A -> CaMKII_12_00 + PP2A');
camkii_reac_6c_k = addkineticlaw(camkii_reac_6c,'MassAction');
set(camkii_reac_6c_k, 'ParameterVariableNames', {'koff_camkii_pp2a'});

camkii_reac_6d = addreaction(pfpc_model, 'CaMKII_12_00_PP2A -> CaMKII_10a_00 + PP2A');
camkii_reac_6d_k = addkineticlaw(camkii_reac_6d,'MassAction');
set(camkii_reac_6d_k, 'ParameterVariableNames', {'kcat_pp2a_cam'});

camkii_reac_6e = addreaction(pfpc_model, 'CaMKII_12_00 + PP2A -> CaMKII_12_00_PP2A');
camkii_reac_6e_k = addkineticlaw(camkii_reac_6e,'MassAction');
set(camkii_reac_6e_k, 'ParameterVariableNames', {'kon_camkii_pp2a'});

camkii_reac_6f = addreaction(pfpc_model, 'CaMKII_12_00_PP2A -> CaMKII_12_00 + PP2A');
camkii_reac_6f_k = addkineticlaw(camkii_reac_6f,'MassAction');
set(camkii_reac_6f_k, 'ParameterVariableNames', {'koff_camkii_pp2a'});

camkii_reac_6g = addreaction(pfpc_model, 'CaMKII_12_00_PP2A -> CaMKII_20a_00 + PP2A');
camkii_reac_6g_k = addkineticlaw(camkii_reac_6g,'MassAction');
set(camkii_reac_6g_k, 'ParameterVariableNames', {'kcat_pp2a_cam'});

camkii_reac_6h = addreaction(pfpc_model, 'CaMKII_20a_00 -> CaMKII_12_00');
camkii_reac_6h_k = addkineticlaw(camkii_reac_6h,'MassAction');
set(camkii_reac_6h_k, 'ParameterVariableNames', {'k2'});

camkii_reac_6i = addreaction(pfpc_model, 'CaMKII_20a_00 -> CaMKII_20_00');
camkii_reac_6i_k = addkineticlaw(camkii_reac_6i,'MassAction');
set(camkii_reac_6i_k, 'ParameterVariableNames', {'kdeact_cam'});

camkii_reac_6j = addreaction(pfpc_model, 'CaMKII_12_00 -> CaMKII_22_00');
camkii_reac_6j_k = addkineticlaw(camkii_reac_6j,'MassAction');
set(camkii_reac_6j_k, 'ParameterVariableNames', {'k3'});

%7
camkii_reac_7a = addreaction(pfpc_model, 'CaMKII_22_00 + PP2A -> CaMKII_22_00_PP2A');
camkii_reac_7a_k = addkineticlaw(camkii_reac_7a,'MassAction');
set(camkii_reac_7a_k, 'ParameterVariableNames', {'kon_camkii_pp2a_2'});

camkii_reac_7b = addreaction(pfpc_model, 'CaMKII_22_00_PP2A -> CaMKII_22_00 + PP2A');
camkii_reac_7b_k = addkineticlaw(camkii_reac_7b,'MassAction');
set(camkii_reac_7b_k, 'ParameterVariableNames', {'koff_camkii_pp2a'});

camkii_reac_7c = addreaction(pfpc_model, 'CaMKII_22_00_PP2A -> CaMKII_20a_00 + PP2A');
camkii_reac_7c_k = addkineticlaw(camkii_reac_7c,'MassAction');
set(camkii_reac_7c_k, 'ParameterVariableNames', {'kcat_pp2a_cam'});

%8
camkii_reac_8a = addreaction(pfpc_model, 'CaMKII_20_00 + Ca4CaM -> CaMKII_20_01');
camkii_reac_8a_k = addkineticlaw(camkii_reac_8a,'MassAction');
set(camkii_reac_8a_k, 'ParameterVariableNames', {'k1'});

camkii_reac_8b = addreaction(pfpc_model, 'CaMKII_20_00 + PP2A -> CaMKII_20_00_PP2A');
camkii_reac_8b_k = addkineticlaw(camkii_reac_8b,'MassAction');
set(camkii_reac_8b_k, 'ParameterVariableNames', {'kon_camkii_pp2a'});

camkii_reac_8c = addreaction(pfpc_model, 'CaMKII_20_00_PP2A -> CaMKII_20_00 + PP2A');
camkii_reac_8c_k = addkineticlaw(camkii_reac_8c,'MassAction');
set(camkii_reac_8c_k, 'ParameterVariableNames', {'koff_camkii_pp2a'});

camkii_reac_8d = addreaction(pfpc_model, 'CaMKII_20_00_PP2A -> CaMKII_00_00 + PP2A');
camkii_reac_8d_k = addkineticlaw(camkii_reac_8d,'MassAction');
set(camkii_reac_8d_k, 'ParameterVariableNames', {'kcat_pp2a_cam'});

%9
camkii_reac_9a = addreaction(pfpc_model, 'CaMKII_10_10 -> CaMKII_10_00 + Ca4CaM');
camkii_reac_9a_k = addkineticlaw(camkii_reac_9a,'MassAction');
set(camkii_reac_9a_k, 'ParameterVariableNames', {'k10r'});

camkii_reac_9b = addreaction(pfpc_model, 'CaMKII_10_10 + Ca4CaM -> CaMKII_10_11');
camkii_reac_9b_k = addkineticlaw(camkii_reac_9b,'MassAction');
set(camkii_reac_9b_k, 'ParameterVariableNames', {'k1'});

camkii_reac_9c = addreaction(pfpc_model, 'CaMKII_10_10 + PP2A -> CaMKII_10_10_PP2A');
camkii_reac_9c_k = addkineticlaw(camkii_reac_9c,'MassAction');
set(camkii_reac_9c_k, 'ParameterVariableNames', {'kon_camkii_pp2a'});

camkii_reac_9d = addreaction(pfpc_model, 'CaMKII_10_10_PP2A -> CaMKII_10_10 + PP2A');
camkii_reac_9d_k = addkineticlaw(camkii_reac_9d,'MassAction');
set(camkii_reac_9d_k, 'ParameterVariableNames', {'koff_camkii_pp2a'});

camkii_reac_9e = addreaction(pfpc_model, 'CaMKII_10_10_PP2A -> CaMKII_00_10 + PP2A');
camkii_reac_9e_k = addkineticlaw(camkii_reac_9e,'MassAction');
set(camkii_reac_9e_k, 'ParameterVariableNames', {'kcat_pp2a_cam'});


%10
camkii_reac_10a = addreaction(pfpc_model, 'CaMKII_10_01 -> CaMKII_10_00 + Ca4CaM');
camkii_reac_10a_k = addkineticlaw(camkii_reac_10a,'MassAction');
set(camkii_reac_10a_k, 'ParameterVariableNames', {'k1r'});

camkii_reac_10b = addreaction(pfpc_model, 'CaMKII_10_01 + Ca4CaM -> CaMKII_10_11');
camkii_reac_10b_k = addkineticlaw(camkii_reac_10b,'MassAction');
set(camkii_reac_10b_k, 'ParameterVariableNames', {'k10'});

camkii_reac_10c = addreaction(pfpc_model, 'CaMKII_10_01 + PP2A -> CaMKII_10_01_PP2A');
camkii_reac_10c_k = addkineticlaw(camkii_reac_10c,'MassAction');
set(camkii_reac_10c_k, 'ParameterVariableNames', {'kon_camkii_pp2a'});

camkii_reac_10d = addreaction(pfpc_model, 'CaMKII_10_01_PP2A -> CaMKII_10_01 + PP2A');
camkii_reac_10d_k = addkineticlaw(camkii_reac_10d,'MassAction');
set(camkii_reac_10d_k, 'ParameterVariableNames', {'koff_camkii_pp2a'});

camkii_reac_10e = addreaction(pfpc_model, 'CaMKII_10_01_PP2A -> CaMKII_00_10 + PP2A');
camkii_reac_10e_k = addkineticlaw(camkii_reac_10e,'MassAction');
set(camkii_reac_10e_k, 'ParameterVariableNames', {'kcat_pp2a_cam'});

camkii_reac_10f = addreaction(pfpc_model, 'CaMKII_10_01 -> CaMKII_20_01');
camkii_reac_10f_k = addkineticlaw(camkii_reac_10f,'MassAction');
set(camkii_reac_10f_k, 'ParameterVariableNames', {'k3'});

%%%%%%%%%%%%NEW%%%%%%%%%%%%
camkii_reac_10g = addreaction(pfpc_model, 'CaMKII_10_01 -> CaMKII_11_10');
camkii_reac_10g_k = addkineticlaw(camkii_reac_10g,'MassAction');
set(camkii_reac_10g_k, 'ParameterVariableNames', {'k2'});

%11
camkii_reac_11a = addreaction(pfpc_model, 'CaMKII_11_10 -> CaMKII_11_00 + Ca4CaM');
camkii_reac_11a_k = addkineticlaw(camkii_reac_11a,'MassAction');
set(camkii_reac_11a_k, 'ParameterVariableNames', {'k10r'});

camkii_reac_11b = addreaction(pfpc_model, 'CaMKII_11_10 + Ca4CaM -> CaMKII_11_11');
camkii_reac_11b_k = addkineticlaw(camkii_reac_11b,'MassAction');
set(camkii_reac_11b_k, 'ParameterVariableNames', {'k10'});

camkii_reac_11c = addreaction(pfpc_model, 'CaMKII_11_10 + PP2A -> CaMKII_11_10_PP2A');
camkii_reac_11c_k = addkineticlaw(camkii_reac_11c,'MassAction');
set(camkii_reac_11c_k, 'ParameterVariableNames', {'kon_camkii_pp2a'});

camkii_reac_11d = addreaction(pfpc_model, 'CaMKII_11_10_PP2A -> CaMKII_11_10 + PP2A');
camkii_reac_11d_k = addkineticlaw(camkii_reac_11d,'MassAction');
set(camkii_reac_11d_k, 'ParameterVariableNames', {'koff_camkii_pp2a'});

camkii_reac_11e = addreaction(pfpc_model, 'CaMKII_11_10_PP2A -> CaMKII_10a_10 + PP2A');
camkii_reac_11e_k = addkineticlaw(camkii_reac_11e,'MassAction');
set(camkii_reac_11e_k, 'ParameterVariableNames', {'kcat_pp2a_cam'});

camkii_reac_11f = addreaction(pfpc_model, 'CaMKII_10a_10 -> CaMKII_11_10');
camkii_reac_11f_k = addkineticlaw(camkii_reac_11f,'MassAction');
set(camkii_reac_11f_k, 'ParameterVariableNames', {'k2'});

camkii_reac_11g = addreaction(pfpc_model, 'CaMKII_10a_10 -> CaMKII_10_10');
camkii_reac_11g_k = addkineticlaw(camkii_reac_11g,'MassAction');
set(camkii_reac_11g_k, 'ParameterVariableNames', {'kdeact_cam'});

camkii_reac_11h = addreaction(pfpc_model, 'CaMKII_11_10 + PP2A -> CaMKII_11_10_PP2A');
camkii_reac_11h_k = addkineticlaw(camkii_reac_11h,'MassAction');
set(camkii_reac_11h_k, 'ParameterVariableNames', {'kon_camkii_pp2a'});

camkii_reac_11i = addreaction(pfpc_model, 'CaMKII_11_10_PP2A -> CaMKII_11_10 + PP2A');
camkii_reac_11i_k = addkineticlaw(camkii_reac_11i,'MassAction');
set(camkii_reac_11i_k, 'ParameterVariableNames', {'koff_camkii_pp2a'});

camkii_reac_11j = addreaction(pfpc_model, 'CaMKII_11_10_PP2A -> CaMKII_10_01 + PP2A');
camkii_reac_11j_k = addkineticlaw(camkii_reac_11j,'MassAction');
set(camkii_reac_11j_k, 'ParameterVariableNames', {'kcat_pp2a_cam'});

camkii_reac_11k = addreaction(pfpc_model, 'CaMKII_11_10 -> CaMKII_12_10');
camkii_reac_11k_k = addkineticlaw(camkii_reac_11k,'MassAction');
set(camkii_reac_11k_k, 'ParameterVariableNames', {'k3'});


%12
camkii_reac_12a = addreaction(pfpc_model, 'CaMKII_20_01 -> CaMKII_20_00 + Ca4CaM');
camkii_reac_12a_k = addkineticlaw(camkii_reac_12a,'MassAction');
set(camkii_reac_12a_k, 'ParameterVariableNames', {'k1r'});

camkii_reac_12b = addreaction(pfpc_model, 'CaMKII_20_01 + PP2A -> CaMKII_20_01_PP2A');
camkii_reac_12b_k = addkineticlaw(camkii_reac_12b,'MassAction');
set(camkii_reac_12b_k, 'ParameterVariableNames', {'kon_camkii_pp2a'});

camkii_reac_12c = addreaction(pfpc_model, 'CaMKII_20_01_PP2A -> CaMKII_20_01 + PP2A');
camkii_reac_12c_k = addkineticlaw(camkii_reac_12c,'MassAction');
set(camkii_reac_12c_k, 'ParameterVariableNames', {'koff_camkii_pp2a'});

camkii_reac_12d = addreaction(pfpc_model, 'CaMKII_20_01_PP2A -> CaMKII_00_10 + PP2A');
camkii_reac_12d_k = addkineticlaw(camkii_reac_12d,'MassAction');
set(camkii_reac_12d_k, 'ParameterVariableNames', {'kcat_pp2a_cam'});

%13
camkii_reac_13a = addreaction(pfpc_model, 'CaMKII_12_10 -> CaMKII_12_00 + Ca4CaM');
camkii_reac_13a_k = addkineticlaw(camkii_reac_13a,'MassAction');
set(camkii_reac_13a_k, 'ParameterVariableNames', {'k10r'});

camkii_reac_13b = addreaction(pfpc_model, 'CaMKII_12_10 + PP2A -> CaMKII_12_10_PP2A');
camkii_reac_13b_k = addkineticlaw(camkii_reac_13b,'MassAction');
set(camkii_reac_13b_k, 'ParameterVariableNames', {'kon_camkii_pp2a'});

camkii_reac_13c = addreaction(pfpc_model, 'CaMKII_12_10_PP2A -> CaMKII_12_10 + PP2A');
camkii_reac_13c_k = addkineticlaw(camkii_reac_13c,'MassAction');
set(camkii_reac_13c_k, 'ParameterVariableNames', {'koff_camkii_pp2a'});

camkii_reac_13d = addreaction(pfpc_model, 'CaMKII_12_10_PP2A -> CaMKII_10a_10 + PP2A');
camkii_reac_13d_k = addkineticlaw(camkii_reac_13d,'MassAction');
set(camkii_reac_13d_k, 'ParameterVariableNames', {'kcat_pp2a_cam'});

camkii_reac_13e = addreaction(pfpc_model, 'CaMKII_12_10 + PP2A -> CaMKII_12_10_PP2A');
camkii_reac_13e_k = addkineticlaw(camkii_reac_13e,'MassAction');
set(camkii_reac_13e_k, 'ParameterVariableNames', {'kon_camkii_pp2a'});

camkii_reac_13f = addreaction(pfpc_model, 'CaMKII_12_10_PP2A -> CaMKII_12_10 + PP2A');
camkii_reac_13f_k = addkineticlaw(camkii_reac_13f,'MassAction');
set(camkii_reac_13f_k, 'ParameterVariableNames', {'koff_camkii_pp2a'});

camkii_reac_13g = addreaction(pfpc_model, 'CaMKII_12_10_PP2A -> CaMKII_20a_01 + PP2A');
camkii_reac_13g_k = addkineticlaw(camkii_reac_13g,'MassAction');
set(camkii_reac_13g_k, 'ParameterVariableNames', {'kcat_pp2a_cam'});

camkii_reac_13h = addreaction(pfpc_model, 'CaMKII_20a_01 -> CaMKII_12_10');
camkii_reac_13h_k = addkineticlaw(camkii_reac_13h,'MassAction');
set(camkii_reac_13h_k, 'ParameterVariableNames', {'k2'});

camkii_reac_13i = addreaction(pfpc_model, 'CaMKII_20a_01 -> CaMKII_20_01');
camkii_reac_13i_k = addkineticlaw(camkii_reac_13i,'MassAction');
set(camkii_reac_13i_k, 'ParameterVariableNames', {'kdeact_cam'});

%%%%%%%%%%NEW%%%%%%%%%%
camkii_reac_13j = addreaction(pfpc_model, 'CaMKII_20_01 -> CaMKII_12_10');
camkii_reac_13j_k = addkineticlaw(camkii_reac_13j,'MassAction');
set(camkii_reac_13j_k, 'ParameterVariableNames', {'k2'});

%14
camkii_reac_14a = addreaction(pfpc_model, 'CaMKII_10_11 -> CaMKII_10_10 + Ca4CaM');
camkii_reac_14a_k = addkineticlaw(camkii_reac_14a,'MassAction');
set(camkii_reac_14a_k, 'ParameterVariableNames', {'k1r'});

camkii_reac_14b = addreaction(pfpc_model, 'CaMKII_10_11 -> CaMKII_10_01 + Ca4CaM');
camkii_reac_14b_k = addkineticlaw(camkii_reac_14b,'MassAction');
set(camkii_reac_14b_k, 'ParameterVariableNames', {'k10r'});

camkii_reac_14c = addreaction(pfpc_model, 'CaMKII_10_11 + PP2A -> CaMKII_10_11_PP2A');
camkii_reac_14c_k = addkineticlaw(camkii_reac_14c,'MassAction');
set(camkii_reac_14c_k, 'ParameterVariableNames', {'kon_camkii_pp2a'});

camkii_reac_14d = addreaction(pfpc_model, 'CaMKII_10_11_PP2A -> CaMKII_10_11 + PP2A');
camkii_reac_14d_k = addkineticlaw(camkii_reac_14d,'MassAction');
set(camkii_reac_14d_k, 'ParameterVariableNames', {'koff_camkii_pp2a'});

camkii_reac_14e = addreaction(pfpc_model, 'CaMKII_10_11_PP2A -> CaMKII_00_11 + PP2A');
camkii_reac_14e_k = addkineticlaw(camkii_reac_14e,'MassAction');
set(camkii_reac_14e_k, 'ParameterVariableNames', {'kcat_pp2a_cam'});

camkii_reac_14f = addreaction(pfpc_model, 'CaMKII_10_11 -> CaMKII_11_11');
camkii_reac_14f_k = addkineticlaw(camkii_reac_14f,'MassAction');
set(camkii_reac_14f_k, 'ParameterVariableNames', {'k2'});

%15
camkii_reac_15a = addreaction(pfpc_model, 'CaMKII_11_11 -> CaMKII_11_10 + Ca4CaM');
camkii_reac_15a_k = addkineticlaw(camkii_reac_15a,'MassAction');
set(camkii_reac_15a_k, 'ParameterVariableNames', {'k10r_2'});

camkii_reac_15b = addreaction(pfpc_model, 'CaMKII_11_11 + PP2A -> CaMKII_11_11_PP2A');
camkii_reac_15b_k = addkineticlaw(camkii_reac_15b,'MassAction');
set(camkii_reac_15b_k, 'ParameterVariableNames', {'kon_camkii_pp2a_2'});

camkii_reac_15c = addreaction(pfpc_model, 'CaMKII_11_11_PP2A -> CaMKII_11_11 + PP2A');
camkii_reac_15c_k = addkineticlaw(camkii_reac_15c,'MassAction');
set(camkii_reac_15c_k, 'ParameterVariableNames', {'koff_camkii_pp2a'});

camkii_reac_15d = addreaction(pfpc_model, 'CaMKII_11_11_PP2A -> CaMKII_10_11 + PP2A');
camkii_reac_15d_k = addkineticlaw(camkii_reac_15d,'MassAction');
set(camkii_reac_15d_k, 'ParameterVariableNames', {'kcat_pp2a_cam'});

%%%%%%%CaMKII PP1 Dephosphorylations%%%%%%%%%

camkii_pp1_reac_1a = addreaction(pfpc_model, 'CaMKII_10_00 + PP1 -> CaMKII_10_00_PP1');
camkii_pp1_reac_1a_k = addkineticlaw(camkii_pp1_reac_1a,'MassAction');
set(camkii_pp1_reac_1a_k, 'ParameterVariableNames', {'kon_camkii_pp1'});

camkii_pp1_reac_1b = addreaction(pfpc_model, 'CaMKII_10_00_PP1 -> CaMKII_10_00 + PP1');
camkii_pp1_reac_1b_k = addkineticlaw(camkii_pp1_reac_1b,'MassAction');
set(camkii_pp1_reac_1b_k, 'ParameterVariableNames', {'koff_camkii_pp1'});

camkii_pp1_reac_1c = addreaction(pfpc_model, 'CaMKII_10_00_PP1 -> CaMKII_00_00 + PP1');
camkii_pp1_reac_1c_k = addkineticlaw(camkii_pp1_reac_1c,'MassAction');
set(camkii_pp1_reac_1c_k, 'ParameterVariableNames', {'kcat_pp1_cam'});


camkii_pp1_reac_2a = addreaction(pfpc_model, 'CaMKII_11_00 + PP1 -> CaMKII_11_00_PP1');
camkii_pp1_reac_2a_k = addkineticlaw(camkii_pp1_reac_2a,'MassAction');
set(camkii_pp1_reac_2a_k, 'ParameterVariableNames', {'kon_camkii_pp1_2'});

camkii_pp1_reac_2b = addreaction(pfpc_model, 'CaMKII_11_00_PP1 -> CaMKII_11_00 + PP1');
camkii_pp1_reac_2b_k = addkineticlaw(camkii_pp1_reac_2b,'MassAction');
set(camkii_pp1_reac_2b_k, 'ParameterVariableNames', {'koff_camkii_pp1'});

camkii_pp1_reac_2c = addreaction(pfpc_model, 'CaMKII_11_00_PP1 -> CaMKII_10a_00 + PP1');
camkii_pp1_reac_2c_k = addkineticlaw(camkii_pp1_reac_2c,'MassAction');
set(camkii_pp1_reac_2c_k, 'ParameterVariableNames', {'kcat_pp1_cam'});


camkii_pp1_reac_3a = addreaction(pfpc_model, 'CaMKII_12_00 + PP1 -> CaMKII_12_00_PP1');
camkii_pp1_reac_3a_k = addkineticlaw(camkii_pp1_reac_3a,'MassAction');
set(camkii_pp1_reac_3a_k, 'ParameterVariableNames', {'kon_camkii_pp1'});

camkii_pp1_reac_3b = addreaction(pfpc_model, 'CaMKII_12_00_PP1 -> CaMKII_12_00 + PP1');
camkii_pp1_reac_3b_k = addkineticlaw(camkii_pp1_reac_3b,'MassAction');
set(camkii_pp1_reac_3b_k, 'ParameterVariableNames', {'koff_camkii_pp1'});

camkii_pp1_reac_3c = addreaction(pfpc_model, 'CaMKII_12_00_PP1 -> CaMKII_10a_00 + PP1');
camkii_pp1_reac_3c_k = addkineticlaw(camkii_pp1_reac_3c,'MassAction');
set(camkii_pp1_reac_3c_k, 'ParameterVariableNames', {'kcat_pp1_cam'});


camkii_pp1_reac_4a = addreaction(pfpc_model, 'CaMKII_12_00 + PP1 -> CaMKII_12_00_PP1');
camkii_pp1_reac_4a_k = addkineticlaw(camkii_pp1_reac_4a,'MassAction');
set(camkii_pp1_reac_4a_k, 'ParameterVariableNames', {'kon_camkii_pp1'});

camkii_pp1_reac_4b = addreaction(pfpc_model, 'CaMKII_12_00_PP1 -> CaMKII_12_00 + PP1');
camkii_pp1_reac_4b_k = addkineticlaw(camkii_pp1_reac_4b,'MassAction');
set(camkii_pp1_reac_4b_k, 'ParameterVariableNames', {'koff_camkii_pp1'});

camkii_pp1_reac_4c = addreaction(pfpc_model, 'CaMKII_12_00_PP1 -> CaMKII_20a_00 + PP1');
camkii_pp1_reac_4c_k = addkineticlaw(camkii_pp1_reac_4c,'MassAction');
set(camkii_pp1_reac_4c_k, 'ParameterVariableNames', {'kcat_pp1_cam'});


camkii_pp1_reac_5a = addreaction(pfpc_model, 'CaMKII_22_00 + PP1 -> CaMKII_22_00_PP1');
camkii_pp1_reac_5a_k = addkineticlaw(camkii_pp1_reac_5a,'MassAction');
set(camkii_pp1_reac_5a_k, 'ParameterVariableNames', {'kon_camkii_pp1_2'});

camkii_pp1_reac_5b = addreaction(pfpc_model, 'CaMKII_22_00_PP1 -> CaMKII_22_00 + PP1');
camkii_pp1_reac_5b_k = addkineticlaw(camkii_pp1_reac_5b,'MassAction');
set(camkii_pp1_reac_5b_k, 'ParameterVariableNames', {'koff_camkii_pp1'});

camkii_pp1_reac_5c = addreaction(pfpc_model, 'CaMKII_22_00_PP1 -> CaMKII_20a_00 + PP1');
camkii_pp1_reac_5c_k = addkineticlaw(camkii_pp1_reac_5c,'MassAction');
set(camkii_pp1_reac_5c_k, 'ParameterVariableNames', {'kcat_pp1_cam'});


camkii_pp1_reac_6a = addreaction(pfpc_model, 'CaMKII_20_00 + PP1 -> CaMKII_20_00_PP1');
camkii_pp1_reac_6a_k = addkineticlaw(camkii_pp1_reac_6a,'MassAction');
set(camkii_pp1_reac_6a_k, 'ParameterVariableNames', {'kon_camkii_pp1'});

camkii_pp1_reac_6b = addreaction(pfpc_model, 'CaMKII_20_00_PP1 -> CaMKII_20_00 + PP1');
camkii_pp1_reac_6b_k = addkineticlaw(camkii_pp1_reac_6b,'MassAction');
set(camkii_pp1_reac_6b_k, 'ParameterVariableNames', {'koff_camkii_pp1'});

camkii_pp1_reac_6c = addreaction(pfpc_model, 'CaMKII_20_00_PP1 -> CaMKII_00_00 + PP1');
camkii_pp1_reac_6c_k = addkineticlaw(camkii_pp1_reac_6c,'MassAction');
set(camkii_pp1_reac_6c_k, 'ParameterVariableNames', {'kcat_pp1_cam'});


camkii_pp1_reac_7a = addreaction(pfpc_model, 'CaMKII_10_10 + PP1 -> CaMKII_10_10_PP1');
camkii_pp1_reac_7a_k = addkineticlaw(camkii_pp1_reac_7a,'MassAction');
set(camkii_pp1_reac_7a_k, 'ParameterVariableNames', {'kon_camkii_pp1'});

camkii_pp1_reac_7b = addreaction(pfpc_model, 'CaMKII_10_10_PP1 -> CaMKII_10_10 + PP1');
camkii_pp1_reac_7b_k = addkineticlaw(camkii_pp1_reac_7b,'MassAction');
set(camkii_pp1_reac_7b_k, 'ParameterVariableNames', {'koff_camkii_pp1'});

camkii_pp1_reac_7c = addreaction(pfpc_model, 'CaMKII_10_10_PP1 -> CaMKII_00_10 + PP1');
camkii_pp1_reac_7c_k = addkineticlaw(camkii_pp1_reac_7c,'MassAction');
set(camkii_pp1_reac_7c_k, 'ParameterVariableNames', {'kcat_pp1_cam'});


camkii_pp1_reac_8a = addreaction(pfpc_model, 'CaMKII_10_01 + PP1 -> CaMKII_10_01_PP1');
camkii_pp1_reac_8a_k = addkineticlaw(camkii_pp1_reac_8a,'MassAction');
set(camkii_pp1_reac_8a_k, 'ParameterVariableNames', {'kon_camkii_pp1'});

camkii_pp1_reac_8b = addreaction(pfpc_model, 'CaMKII_10_01_PP1 -> CaMKII_10_01 + PP1');
camkii_pp1_reac_8b_k = addkineticlaw(camkii_pp1_reac_8b,'MassAction');
set(camkii_pp1_reac_8b_k, 'ParameterVariableNames', {'koff_camkii_pp1'});

camkii_pp1_reac_8c = addreaction(pfpc_model, 'CaMKII_10_01_PP1 -> CaMKII_00_10 + PP1');
camkii_pp1_reac_8c_k = addkineticlaw(camkii_pp1_reac_8c,'MassAction');
set(camkii_pp1_reac_8c_k, 'ParameterVariableNames', {'kcat_pp1_cam'});


camkii_pp1_reac_9a = addreaction(pfpc_model, 'CaMKII_11_10 + PP1 -> CaMKII_11_10_PP1');
camkii_pp1_reac_9a_k = addkineticlaw(camkii_pp1_reac_9a,'MassAction');
set(camkii_pp1_reac_9a_k, 'ParameterVariableNames', {'kon_camkii_pp1'});

camkii_pp1_reac_9b = addreaction(pfpc_model, 'CaMKII_11_10_PP1 -> CaMKII_11_10 + PP1');
camkii_pp1_reac_9b_k = addkineticlaw(camkii_pp1_reac_9b,'MassAction');
set(camkii_pp1_reac_9b_k, 'ParameterVariableNames', {'koff_camkii_pp1'});

camkii_pp1_reac_9c = addreaction(pfpc_model, 'CaMKII_11_10_PP1 -> CaMKII_10a_10 + PP1');
camkii_pp1_reac_9c_k = addkineticlaw(camkii_pp1_reac_9c,'MassAction');
set(camkii_pp1_reac_9c_k, 'ParameterVariableNames', {'kcat_pp1_cam'});


camkii_pp1_reac_10a = addreaction(pfpc_model, 'CaMKII_11_10 + PP1 -> CaMKII_11_10_PP1');
camkii_pp1_reac_10a_k = addkineticlaw(camkii_pp1_reac_10a,'MassAction');
set(camkii_pp1_reac_10a_k, 'ParameterVariableNames', {'kon_camkii_pp1'});

camkii_pp1_reac_10b = addreaction(pfpc_model, 'CaMKII_11_10_PP1 -> CaMKII_11_10 + PP1');
camkii_pp1_reac_10b_k = addkineticlaw(camkii_pp1_reac_10b,'MassAction');
set(camkii_pp1_reac_10b_k, 'ParameterVariableNames', {'koff_camkii_pp1'});

camkii_pp1_reac_10c = addreaction(pfpc_model, 'CaMKII_11_10_PP1 -> CaMKII_10_01 + PP1');
camkii_pp1_reac_10c_k = addkineticlaw(camkii_pp1_reac_10c,'MassAction');
set(camkii_pp1_reac_10c_k, 'ParameterVariableNames', {'kcat_pp1_cam'});


camkii_pp1_reac_11a = addreaction(pfpc_model, 'CaMKII_20_01 + PP1 -> CaMKII_20_01_PP1');
camkii_pp1_reac_11a_k = addkineticlaw(camkii_pp1_reac_11a,'MassAction');
set(camkii_pp1_reac_11a_k, 'ParameterVariableNames', {'kon_camkii_pp1'});

camkii_pp1_reac_11b = addreaction(pfpc_model, 'CaMKII_20_01_PP1 -> CaMKII_20_01 + PP1');
camkii_pp1_reac_11b_k = addkineticlaw(camkii_pp1_reac_11b,'MassAction');
set(camkii_pp1_reac_11b_k, 'ParameterVariableNames', {'koff_camkii_pp1'});

camkii_pp1_reac_11c = addreaction(pfpc_model, 'CaMKII_20_01_PP1 -> CaMKII_00_10 + PP1');
camkii_pp1_reac_11c_k = addkineticlaw(camkii_pp1_reac_11c,'MassAction');
set(camkii_pp1_reac_11c_k, 'ParameterVariableNames', {'kcat_pp1_cam'});


camkii_pp1_reac_12a = addreaction(pfpc_model, 'CaMKII_12_10 + PP1 -> CaMKII_12_10_PP1');
camkii_pp1_reac_12a_k = addkineticlaw(camkii_pp1_reac_12a,'MassAction');
set(camkii_pp1_reac_12a_k, 'ParameterVariableNames', {'kon_camkii_pp1'});

camkii_pp1_reac_12b = addreaction(pfpc_model, 'CaMKII_12_10_PP1 -> CaMKII_12_10 + PP1');
camkii_pp1_reac_12b_k = addkineticlaw(camkii_pp1_reac_12b,'MassAction');
set(camkii_pp1_reac_12b_k, 'ParameterVariableNames', {'koff_camkii_pp1'});

camkii_pp1_reac_12c = addreaction(pfpc_model, 'CaMKII_12_10_PP1 -> CaMKII_10a_10 + PP1');
camkii_pp1_reac_12c_k = addkineticlaw(camkii_pp1_reac_12c,'MassAction');
set(camkii_pp1_reac_12c_k, 'ParameterVariableNames', {'kcat_pp1_cam'});


camkii_pp1_reac_13a = addreaction(pfpc_model, 'CaMKII_12_10 + PP1 -> CaMKII_12_10_PP1');
camkii_pp1_reac_13a_k = addkineticlaw(camkii_pp1_reac_13a,'MassAction');
set(camkii_pp1_reac_13a_k, 'ParameterVariableNames', {'kon_camkii_pp1'});

camkii_pp1_reac_13b = addreaction(pfpc_model, 'CaMKII_12_10_PP1 -> CaMKII_12_10 + PP1');
camkii_pp1_reac_13b_k = addkineticlaw(camkii_pp1_reac_13b,'MassAction');
set(camkii_pp1_reac_13b_k, 'ParameterVariableNames', {'koff_camkii_pp1'});

camkii_pp1_reac_13c = addreaction(pfpc_model, 'CaMKII_12_10_PP1 -> CaMKII_20a_01 + PP1');
camkii_pp1_reac_13c_k = addkineticlaw(camkii_pp1_reac_13c,'MassAction');
set(camkii_pp1_reac_13c_k, 'ParameterVariableNames', {'kcat_pp1_cam'});


camkii_pp1_reac_14a = addreaction(pfpc_model, 'CaMKII_10_11 + PP1 -> CaMKII_10_11_PP1');
camkii_pp1_reac_14a_k = addkineticlaw(camkii_pp1_reac_14a,'MassAction');
set(camkii_pp1_reac_14a_k, 'ParameterVariableNames', {'kon_camkii_pp1'});

camkii_pp1_reac_14b = addreaction(pfpc_model, 'CaMKII_10_11_PP1 -> CaMKII_10_11 + PP1');
camkii_pp1_reac_14b_k = addkineticlaw(camkii_pp1_reac_14b,'MassAction');
set(camkii_pp1_reac_14b_k, 'ParameterVariableNames', {'koff_camkii_pp1'});

camkii_pp1_reac_14c = addreaction(pfpc_model, 'CaMKII_10_11_PP1 -> CaMKII_00_11 + PP1');
camkii_pp1_reac_14c_k = addkineticlaw(camkii_pp1_reac_14c,'MassAction');
set(camkii_pp1_reac_14c_k, 'ParameterVariableNames', {'kcat_pp1_cam'});


camkii_pp1_reac_15a = addreaction(pfpc_model, 'CaMKII_11_11 + PP1 -> CaMKII_11_11_PP1');
camkii_pp1_reac_15a_k = addkineticlaw(camkii_pp1_reac_15a,'MassAction');
set(camkii_pp1_reac_15a_k, 'ParameterVariableNames', {'kon_camkii_pp1_2'});

camkii_pp1_reac_15b = addreaction(pfpc_model, 'CaMKII_11_11_PP1 -> CaMKII_11_11 + PP1');
camkii_pp1_reac_15b_k = addkineticlaw(camkii_pp1_reac_15b,'MassAction');
set(camkii_pp1_reac_15b_k, 'ParameterVariableNames', {'koff_camkii_pp1'});

camkii_pp1_reac_15c = addreaction(pfpc_model, 'CaMKII_11_11_PP1 -> CaMKII_10_11 + PP1');
camkii_pp1_reac_15c_k = addkineticlaw(camkii_pp1_reac_15c,'MassAction');
set(camkii_pp1_reac_15c_k, 'ParameterVariableNames', {'kcat_pp1_cam'});





%%%%%%%PDE1 Activation and Regulation%%%%%%%
%Add species
%PDE-CaM Parameters
kon_cam_pde = addparameter(pfpc_model,'kon_cam_pde', 'Value', 1200, 'ValueUnits', '1/(micromolarity*second)');
koff_cam_pde = addparameter(pfpc_model,'koff_cam_pde', 'Value', 0.01, 'ValueUnits', '1/second');
koff_cam_pdep = addparameter(pfpc_model,'koff_cam_pdep', 'Value', 0.36, 'ValueUnits', '1/second');

%Regulation of PDE1 by CaMKII Parameters
kon_camkii_pde = addparameter(pfpc_model,'kon_camkii_pde', 'Value', 0.45, 'ValueUnits', '1/(micromolarity*second)', 'ConstantValue', false);
koff_camkii_pde = addparameter(pfpc_model,'koff_camkii_pde', 'Value', 4, 'ValueUnits', '1/second');
kcat_camkii = addparameter(pfpc_model,'kcat_camkii', 'Value', 2, 'ValueUnits', '1/second');

kdephos_pde = addparameter(pfpc_model,'kdephos_pde ', 'Value', 2, 'ValueUnits', '1/second');

%binding and activation by Ca4CaM
campde_reac_1 = addreaction(pfpc_model, 'PDE1 + Ca4CaM -> PDE1_CaM');
campde_reac_1_k = addkineticlaw(campde_reac_1,'MassAction');
set(campde_reac_1_k, 'ParameterVariableNames', {'kon_cam_pde'});

campde_reac_2 = addreaction(pfpc_model, 'PDE1_CaM -> PDE1 + Ca4CaM');
campde_reac_2_k = addkineticlaw(campde_reac_2,'MassAction');
set(campde_reac_2_k, 'ParameterVariableNames', {'koff_cam_pde'});

campde_reac_3 = addreaction(pfpc_model, 'PDE1p + Ca4CaM -> PDE1p_CaM');
campde_reac_3_k = addkineticlaw(campde_reac_3,'MassAction');
set(campde_reac_3_k, 'ParameterVariableNames', {'kon_cam_pde'});

campde_reac_4 = addreaction(pfpc_model, 'PDE1p_CaM -> PDE1p + Ca4CaM');
campde_reac_4_k = addkineticlaw(campde_reac_4,'MassAction');
set(campde_reac_4_k, 'ParameterVariableNames', {'koff_cam_pdep'});

%Regulation of PDE1 by CaMKII Reactions

%Free PDE1
%CaMKII_00_10
camkiipde_reac_1a = addreaction(pfpc_model, 'CaMKII_00_10 + PDE1 -> CaMKII_00_10_PDE1');
camkiipde_reac_1a_k = addkineticlaw(camkiipde_reac_1a,'MassAction');
set(camkiipde_reac_1a_k, 'ParameterVariableNames', {'kon_camkii_pde'});

camkiipde_reac_1b = addreaction(pfpc_model, 'CaMKII_00_10_PDE1 -> CaMKII_00_10 + PDE1');
camkiipde_reac_1b_k = addkineticlaw(camkiipde_reac_1b,'MassAction');
set(camkiipde_reac_1b_k, 'ParameterVariableNames', {'koff_camkii_pde'});

camkiipde_reac_1c = addreaction(pfpc_model, 'CaMKII_00_10_PDE1 -> CaMKII_00_10 + PDE1p');
camkiipde_reac_1c_k = addkineticlaw(camkiipde_reac_1c,'MassAction');
set(camkiipde_reac_1c_k, 'ParameterVariableNames', {'kcat_camkii'});

%CaMKII_00_11
camkiipde_reac_2a = addreaction(pfpc_model, 'CaMKII_00_11 + PDE1 -> CaMKII_00_11_PDE1');
camkiipde_reac_2a_k = addkineticlaw(camkiipde_reac_2a,'MassAction');
set(camkiipde_reac_2a_k, 'ParameterVariableNames', {'kon_camkii_pde'});

camkiipde_reac_2b = addreaction(pfpc_model, 'CaMKII_00_11_PDE1 -> CaMKII_00_11 + PDE1');
camkiipde_reac_2b_k = addkineticlaw(camkiipde_reac_2b,'MassAction');
set(camkiipde_reac_2b_k, 'ParameterVariableNames', {'koff_camkii_pde'});

camkiipde_reac_2c = addreaction(pfpc_model, 'CaMKII_00_11_PDE1 -> CaMKII_00_11 + PDE1p');
camkiipde_reac_2c_k = addkineticlaw(camkiipde_reac_2c,'MassAction');
set(camkiipde_reac_2c_k, 'ParameterVariableNames', {'kcat_camkii'});

%CaMKII_10_00
camkiipde_reac_3a = addreaction(pfpc_model, 'CaMKII_10_00 + PDE1 -> CaMKII_10_00_PDE1');
camkiipde_reac_3a_k = addkineticlaw(camkiipde_reac_3a,'MassAction');
set(camkiipde_reac_3a_k, 'ParameterVariableNames', {'kon_camkii_pde'});

camkiipde_reac_3b = addreaction(pfpc_model, 'CaMKII_10_00_PDE1 -> CaMKII_10_00 + PDE1');
camkiipde_reac_3b_k = addkineticlaw(camkiipde_reac_3b,'MassAction');
set(camkiipde_reac_3b_k, 'ParameterVariableNames', {'koff_camkii_pde'});

camkiipde_reac_3c = addreaction(pfpc_model, 'CaMKII_10_00_PDE1 -> CaMKII_10_00 + PDE1p');
camkiipde_reac_3c_k = addkineticlaw(camkiipde_reac_3c,'MassAction');
set(camkiipde_reac_3c_k, 'ParameterVariableNames', {'kcat_camkii'});

%CaMKII_20_00
camkiipde_reac_4a = addreaction(pfpc_model, 'CaMKII_20_00 + PDE1 -> CaMKII_20_00_PDE1');
camkiipde_reac_4a_k = addkineticlaw(camkiipde_reac_4a,'MassAction');
set(camkiipde_reac_4a_k, 'ParameterVariableNames', {'kon_camkii_pde'});

camkiipde_reac_4b = addreaction(pfpc_model, 'CaMKII_20_00_PDE1 -> CaMKII_20_00 + PDE1');
camkiipde_reac_4b_k = addkineticlaw(camkiipde_reac_4b,'MassAction');
set(camkiipde_reac_4b_k, 'ParameterVariableNames', {'koff_camkii_pde'});

camkiipde_reac_4c = addreaction(pfpc_model, 'CaMKII_20_00_PDE1 -> CaMKII_20_00 + PDE1p');
camkiipde_reac_4c_k = addkineticlaw(camkiipde_reac_4c,'MassAction');
set(camkiipde_reac_4c_k, 'ParameterVariableNames', {'kcat_camkii'});

%CaMKII_10_10
camkiipde_reac_5a = addreaction(pfpc_model, 'CaMKII_10_10 + PDE1 -> CaMKII_10_10_PDE1');
camkiipde_reac_5a_k = addkineticlaw(camkiipde_reac_5a,'MassAction');
set(camkiipde_reac_5a_k, 'ParameterVariableNames', {'kon_camkii_pde'});

camkiipde_reac_5b = addreaction(pfpc_model, 'CaMKII_10_10_PDE1 -> CaMKII_10_10 + PDE1');
camkiipde_reac_5b_k = addkineticlaw(camkiipde_reac_5b,'MassAction');
set(camkiipde_reac_5b_k, 'ParameterVariableNames', {'koff_camkii_pde'});

camkiipde_reac_5c = addreaction(pfpc_model, 'CaMKII_10_10_PDE1 -> CaMKII_10_10 + PDE1p');
camkiipde_reac_5c_k = addkineticlaw(camkiipde_reac_5c,'MassAction');
set(camkiipde_reac_5c_k, 'ParameterVariableNames', {'kcat_camkii'});

%CaMKII_10_01
camkiipde_reac_6a = addreaction(pfpc_model, 'CaMKII_10_01 + PDE1 -> CaMKII_10_01_PDE1');
camkiipde_reac_6a_k = addkineticlaw(camkiipde_reac_6a,'MassAction');
set(camkiipde_reac_6a_k, 'ParameterVariableNames', {'kon_camkii_pde'});

camkiipde_reac_6b = addreaction(pfpc_model, 'CaMKII_10_01_PDE1 -> CaMKII_10_01 + PDE1');
camkiipde_reac_6b_k = addkineticlaw(camkiipde_reac_6b,'MassAction');
set(camkiipde_reac_6b_k, 'ParameterVariableNames', {'koff_camkii_pde'});

camkiipde_reac_6c = addreaction(pfpc_model, 'CaMKII_10_01_PDE1 -> CaMKII_10_01 + PDE1p');
camkiipde_reac_6c_k = addkineticlaw(camkiipde_reac_6c,'MassAction');
set(camkiipde_reac_6c_k, 'ParameterVariableNames', {'kcat_camkii'});

%CaMKII_20_01
camkiipde_reac_7a = addreaction(pfpc_model, 'CaMKII_20_01 + PDE1 -> CaMKII_20_01_PDE1');
camkiipde_reac_7a_k = addkineticlaw(camkiipde_reac_7a,'MassAction');
set(camkiipde_reac_7a_k, 'ParameterVariableNames', {'kon_camkii_pde'});

camkiipde_reac_7b = addreaction(pfpc_model, 'CaMKII_20_01_PDE1 -> CaMKII_20_01 + PDE1');
camkiipde_reac_7b_k = addkineticlaw(camkiipde_reac_7b,'MassAction');
set(camkiipde_reac_7b_k, 'ParameterVariableNames', {'koff_camkii_pde'});

camkiipde_reac_7c = addreaction(pfpc_model, 'CaMKII_20_01_PDE1 -> CaMKII_20_01 + PDE1p');
camkiipde_reac_7c_k = addkineticlaw(camkiipde_reac_7c,'MassAction');
set(camkiipde_reac_7c_k, 'ParameterVariableNames', {'kcat_camkii'});

%CaMKII_10_11
camkiipde_reac_8a = addreaction(pfpc_model, 'CaMKII_10_11 + PDE1 -> CaMKII_10_11_PDE1');
camkiipde_reac_8a_k = addkineticlaw(camkiipde_reac_8a,'MassAction');
set(camkiipde_reac_8a_k, 'ParameterVariableNames', {'kon_camkii_pde'});

camkiipde_reac_8b = addreaction(pfpc_model, 'CaMKII_10_11_PDE1 -> CaMKII_10_11 + PDE1');
camkiipde_reac_8b_k = addkineticlaw(camkiipde_reac_8b,'MassAction');
set(camkiipde_reac_8b_k, 'ParameterVariableNames', {'koff_camkii_pde'});

camkiipde_reac_8c = addreaction(pfpc_model, 'CaMKII_10_11_PDE1 -> CaMKII_10_11 + PDE1p');
camkiipde_reac_8c_k = addkineticlaw(camkiipde_reac_8c,'MassAction');
set(camkiipde_reac_8c_k, 'ParameterVariableNames', {'kcat_camkii'});

%CaMKII_11_11
camkiipde_reac_9a = addreaction(pfpc_model, 'CaMKII_11_11 + PDE1 -> CaMKII_11_11_PDE1');
camkiipde_reac_9a_k = addkineticlaw(camkiipde_reac_9a,'MassAction');
set(camkiipde_reac_9a_k, 'ParameterVariableNames', {'kon_camkii_pde'});

camkiipde_reac_9b = addreaction(pfpc_model, 'CaMKII_11_11_PDE1 -> CaMKII_11_11 + PDE1');
camkiipde_reac_9b_k = addkineticlaw(camkiipde_reac_9b,'MassAction');
set(camkiipde_reac_9b_k, 'ParameterVariableNames', {'koff_camkii_pde'});

camkiipde_reac_9c = addreaction(pfpc_model, 'CaMKII_11_11_PDE1 -> CaMKII_11_11 + PDE1p');
camkiipde_reac_9c_k = addkineticlaw(camkiipde_reac_9c,'MassAction');
set(camkiipde_reac_9c_k, 'ParameterVariableNames', {'kcat_camkii'});

%CaMKII_11_00
camkiipde_reac_10a = addreaction(pfpc_model, 'CaMKII_11_00 + PDE1 -> CaMKII_11_00_PDE1');
camkiipde_reac_10a_k = addkineticlaw(camkiipde_reac_10a,'MassAction');
set(camkiipde_reac_10a_k, 'ParameterVariableNames', {'kon_camkii_pde'});

camkiipde_reac_10b = addreaction(pfpc_model, 'CaMKII_11_00_PDE1 -> CaMKII_11_00 + PDE1');
camkiipde_reac_10b_k = addkineticlaw(camkiipde_reac_10b,'MassAction');
set(camkiipde_reac_10b_k, 'ParameterVariableNames', {'koff_camkii_pde'});

camkiipde_reac_10c = addreaction(pfpc_model, 'CaMKII_11_00_PDE1 -> CaMKII_11_00 + PDE1p');
camkiipde_reac_10c_k = addkineticlaw(camkiipde_reac_10c,'MassAction');
set(camkiipde_reac_10c_k, 'ParameterVariableNames', {'kcat_camkii'});

%CaMKII_12_00
camkiipde_reac_11a = addreaction(pfpc_model, 'CaMKII_12_00 + PDE1 -> CaMKII_12_00_PDE1');
camkiipde_reac_11a_k = addkineticlaw(camkiipde_reac_11a,'MassAction');
set(camkiipde_reac_11a_k, 'ParameterVariableNames', {'kon_camkii_pde'});

camkiipde_reac_11b = addreaction(pfpc_model, 'CaMKII_12_00_PDE1 -> CaMKII_12_00 + PDE1');
camkiipde_reac_11b_k = addkineticlaw(camkiipde_reac_11b,'MassAction');
set(camkiipde_reac_11b_k, 'ParameterVariableNames', {'koff_camkii_pde'});

camkiipde_reac_11c = addreaction(pfpc_model, 'CaMKII_12_00_PDE1 -> CaMKII_12_00 + PDE1p');
camkiipde_reac_11c_k = addkineticlaw(camkiipde_reac_11c,'MassAction');
set(camkiipde_reac_11c_k, 'ParameterVariableNames', {'kcat_camkii'});

%CaMKII_22_00
camkiipde_reac_12a = addreaction(pfpc_model, 'CaMKII_22_00 + PDE1 -> CaMKII_22_00_PDE1');
camkiipde_reac_12a_k = addkineticlaw(camkiipde_reac_12a,'MassAction');
set(camkiipde_reac_12a_k, 'ParameterVariableNames', {'kon_camkii_pde'});

camkiipde_reac_12b = addreaction(pfpc_model, 'CaMKII_22_00_PDE1 -> CaMKII_22_00 + PDE1');
camkiipde_reac_12b_k = addkineticlaw(camkiipde_reac_12b,'MassAction');
set(camkiipde_reac_12b_k, 'ParameterVariableNames', {'koff_camkii_pde'});

camkiipde_reac_12c = addreaction(pfpc_model, 'CaMKII_22_00_PDE1 -> CaMKII_22_00 + PDE1p');
camkiipde_reac_12c_k = addkineticlaw(camkiipde_reac_12c,'MassAction');
set(camkiipde_reac_12c_k, 'ParameterVariableNames', {'kcat_camkii'});

%CaMKII_11_10
camkiipde_reac_13a = addreaction(pfpc_model, 'CaMKII_11_10 + PDE1 -> CaMKII_11_10_PDE1');
camkiipde_reac_13a_k = addkineticlaw(camkiipde_reac_13a,'MassAction');
set(camkiipde_reac_13a_k, 'ParameterVariableNames', {'kon_camkii_pde'});

camkiipde_reac_13b = addreaction(pfpc_model, 'CaMKII_11_10_PDE1 -> CaMKII_11_10 + PDE1');
camkiipde_reac_13b_k = addkineticlaw(camkiipde_reac_13b,'MassAction');
set(camkiipde_reac_13b_k, 'ParameterVariableNames', {'koff_camkii_pde'});

camkiipde_reac_13c = addreaction(pfpc_model, 'CaMKII_11_10_PDE1 -> CaMKII_11_10 + PDE1p');
camkiipde_reac_13c_k = addkineticlaw(camkiipde_reac_13c,'MassAction');
set(camkiipde_reac_13c_k, 'ParameterVariableNames', {'kcat_camkii'});

%CaMKII_12_10
camkiipde_reac_14a = addreaction(pfpc_model, 'CaMKII_12_10 + PDE1 -> CaMKII_12_10_PDE1');
camkiipde_reac_14a_k = addkineticlaw(camkiipde_reac_14a,'MassAction');
set(camkiipde_reac_14a_k, 'ParameterVariableNames', {'kon_camkii_pde'});

camkiipde_reac_14b = addreaction(pfpc_model, 'CaMKII_12_10_PDE1 -> CaMKII_12_10 + PDE1');
camkiipde_reac_14b_k = addkineticlaw(camkiipde_reac_14b,'MassAction');
set(camkiipde_reac_14b_k, 'ParameterVariableNames', {'koff_camkii_pde'});

camkiipde_reac_14c = addreaction(pfpc_model, 'CaMKII_12_10_PDE1 -> CaMKII_12_10 + PDE1p');
camkiipde_reac_14c_k = addkineticlaw(camkiipde_reac_14c,'MassAction');
set(camkiipde_reac_14c_k, 'ParameterVariableNames', {'kcat_camkii'});

%PDE1_CaM
%CaMKII_00_10
camkiipde_reac_15a = addreaction(pfpc_model, 'CaMKII_00_10 + PDE1_CaM -> CaMKII_00_10_PDE1_CaM');
camkiipde_reac_15a_k = addkineticlaw(camkiipde_reac_15a,'MassAction');
set(camkiipde_reac_15a_k, 'ParameterVariableNames', {'kon_camkii_pde'});

camkiipde_reac_15b = addreaction(pfpc_model, 'CaMKII_00_10_PDE1_CaM -> CaMKII_00_10 + PDE1_CaM');
camkiipde_reac_15b_k = addkineticlaw(camkiipde_reac_15b,'MassAction');
set(camkiipde_reac_15b_k, 'ParameterVariableNames', {'koff_camkii_pde'});

camkiipde_reac_15c = addreaction(pfpc_model, 'CaMKII_00_10_PDE1_CaM -> CaMKII_00_10 + PDE1p_CaM');
camkiipde_reac_15c_k = addkineticlaw(camkiipde_reac_15c,'MassAction');
set(camkiipde_reac_15c_k, 'ParameterVariableNames', {'kcat_camkii'});

%CaMKII_00_11
camkiipde_reac_16a = addreaction(pfpc_model, 'CaMKII_00_11 + PDE1_CaM -> CaMKII_00_11_PDE1_CaM');
camkiipde_reac_16a_k = addkineticlaw(camkiipde_reac_16a,'MassAction');
set(camkiipde_reac_16a_k, 'ParameterVariableNames', {'kon_camkii_pde'});

camkiipde_reac_16b = addreaction(pfpc_model, 'CaMKII_00_11_PDE1_CaM -> CaMKII_00_11 + PDE1_CaM');
camkiipde_reac_16b_k = addkineticlaw(camkiipde_reac_16b,'MassAction');
set(camkiipde_reac_16b_k, 'ParameterVariableNames', {'koff_camkii_pde'});

camkiipde_reac_16c = addreaction(pfpc_model, 'CaMKII_00_11_PDE1_CaM -> CaMKII_00_11 + PDE1p_CaM');
camkiipde_reac_16c_k = addkineticlaw(camkiipde_reac_16c,'MassAction');
set(camkiipde_reac_16c_k, 'ParameterVariableNames', {'kcat_camkii'});

%CaMKII_10_00
camkiipde_reac_17a = addreaction(pfpc_model, 'CaMKII_10_00 + PDE1_CaM -> CaMKII_10_00_PDE1_CaM');
camkiipde_reac_17a_k = addkineticlaw(camkiipde_reac_17a,'MassAction');
set(camkiipde_reac_17a_k, 'ParameterVariableNames', {'kon_camkii_pde'});

camkiipde_reac_17b = addreaction(pfpc_model, 'CaMKII_10_00_PDE1_CaM -> CaMKII_10_00 + PDE1_CaM');
camkiipde_reac_17b_k = addkineticlaw(camkiipde_reac_17b,'MassAction');
set(camkiipde_reac_17b_k, 'ParameterVariableNames', {'koff_camkii_pde'});

camkiipde_reac_17c = addreaction(pfpc_model, 'CaMKII_10_00_PDE1_CaM -> CaMKII_10_00 + PDE1p_CaM');
camkiipde_reac_17c_k = addkineticlaw(camkiipde_reac_17c,'MassAction');
set(camkiipde_reac_17c_k, 'ParameterVariableNames', {'kcat_camkii'});

%CaMKII_20_00
camkiipde_reac_18a = addreaction(pfpc_model, 'CaMKII_20_00 + PDE1_CaM -> CaMKII_20_00_PDE1_CaM');
camkiipde_reac_18a_k = addkineticlaw(camkiipde_reac_18a,'MassAction');
set(camkiipde_reac_18a_k, 'ParameterVariableNames', {'kon_camkii_pde'});

camkiipde_reac_18b = addreaction(pfpc_model, 'CaMKII_20_00_PDE1_CaM -> CaMKII_20_00 + PDE1_CaM');
camkiipde_reac_18b_k = addkineticlaw(camkiipde_reac_18b,'MassAction');
set(camkiipde_reac_18b_k, 'ParameterVariableNames', {'koff_camkii_pde'});

camkiipde_reac_18c = addreaction(pfpc_model, 'CaMKII_20_00_PDE1_CaM -> CaMKII_20_00 + PDE1p_CaM');
camkiipde_reac_18c_k = addkineticlaw(camkiipde_reac_18c,'MassAction');
set(camkiipde_reac_18c_k, 'ParameterVariableNames', {'kcat_camkii'});

%CaMKII_10_10
camkiipde_reac_19a = addreaction(pfpc_model, 'CaMKII_10_10 + PDE1_CaM -> CaMKII_10_10_PDE1_CaM');
camkiipde_reac_19a_k = addkineticlaw(camkiipde_reac_19a,'MassAction');
set(camkiipde_reac_19a_k, 'ParameterVariableNames', {'kon_camkii_pde'});

camkiipde_reac_19b = addreaction(pfpc_model, 'CaMKII_10_10_PDE1_CaM -> CaMKII_10_10 + PDE1_CaM');
camkiipde_reac_19b_k = addkineticlaw(camkiipde_reac_19b,'MassAction');
set(camkiipde_reac_19b_k, 'ParameterVariableNames', {'koff_camkii_pde'});

camkiipde_reac_19c = addreaction(pfpc_model, 'CaMKII_10_10_PDE1_CaM -> CaMKII_10_10 + PDE1p_CaM');
camkiipde_reac_19c_k = addkineticlaw(camkiipde_reac_19c,'MassAction');
set(camkiipde_reac_19c_k, 'ParameterVariableNames', {'kcat_camkii'});

%CaMKII_10_01
camkiipde_reac_20a = addreaction(pfpc_model, 'CaMKII_10_01 + PDE1_CaM -> CaMKII_10_01_PDE1_CaM');
camkiipde_reac_20a_k = addkineticlaw(camkiipde_reac_20a,'MassAction');
set(camkiipde_reac_20a_k, 'ParameterVariableNames', {'kon_camkii_pde'});

camkiipde_reac_20b = addreaction(pfpc_model, 'CaMKII_10_01_PDE1_CaM -> CaMKII_10_01 + PDE1_CaM');
camkiipde_reac_20b_k = addkineticlaw(camkiipde_reac_20b,'MassAction');
set(camkiipde_reac_20b_k, 'ParameterVariableNames', {'koff_camkii_pde'});

camkiipde_reac_20c = addreaction(pfpc_model, 'CaMKII_10_01_PDE1_CaM -> CaMKII_10_01 + PDE1p_CaM');
camkiipde_reac_20c_k = addkineticlaw(camkiipde_reac_20c,'MassAction');
set(camkiipde_reac_20c_k, 'ParameterVariableNames', {'kcat_camkii'});

%CaMKII_20_01
camkiipde_reac_21a = addreaction(pfpc_model, 'CaMKII_20_01 + PDE1_CaM -> CaMKII_20_01_PDE1_CaM');
camkiipde_reac_21a_k = addkineticlaw(camkiipde_reac_21a,'MassAction');
set(camkiipde_reac_21a_k, 'ParameterVariableNames', {'kon_camkii_pde'});

camkiipde_reac_21b = addreaction(pfpc_model, 'CaMKII_20_01_PDE1_CaM -> CaMKII_20_01 + PDE1_CaM');
camkiipde_reac_21b_k = addkineticlaw(camkiipde_reac_21b,'MassAction');
set(camkiipde_reac_21b_k, 'ParameterVariableNames', {'koff_camkii_pde'});

camkiipde_reac_21c = addreaction(pfpc_model, 'CaMKII_20_01_PDE1_CaM -> CaMKII_20_01 + PDE1p_CaM');
camkiipde_reac_21c_k = addkineticlaw(camkiipde_reac_21c,'MassAction');
set(camkiipde_reac_21c_k, 'ParameterVariableNames', {'kcat_camkii'});

%CaMKII_10_11
camkiipde_reac_22a = addreaction(pfpc_model, 'CaMKII_10_11 + PDE1_CaM -> CaMKII_10_11_PDE1_CaM');
camkiipde_reac_22a_k = addkineticlaw(camkiipde_reac_22a,'MassAction');
set(camkiipde_reac_22a_k, 'ParameterVariableNames', {'kon_camkii_pde'});

camkiipde_reac_22b = addreaction(pfpc_model, 'CaMKII_10_11_PDE1_CaM -> CaMKII_10_11 + PDE1_CaM');
camkiipde_reac_22b_k = addkineticlaw(camkiipde_reac_22b,'MassAction');
set(camkiipde_reac_22b_k, 'ParameterVariableNames', {'koff_camkii_pde'});

camkiipde_reac_22c = addreaction(pfpc_model, 'CaMKII_10_11_PDE1_CaM -> CaMKII_10_11 + PDE1p_CaM');
camkiipde_reac_22c_k = addkineticlaw(camkiipde_reac_22c,'MassAction');
set(camkiipde_reac_22c_k, 'ParameterVariableNames', {'kcat_camkii'});

%CaMKII_11_11
camkiipde_reac_23a = addreaction(pfpc_model, 'CaMKII_11_11 + PDE1_CaM -> CaMKII_11_11_PDE1_CaM');
camkiipde_reac_23a_k = addkineticlaw(camkiipde_reac_23a,'MassAction');
set(camkiipde_reac_23a_k, 'ParameterVariableNames', {'kon_camkii_pde'});

camkiipde_reac_23b = addreaction(pfpc_model, 'CaMKII_11_11_PDE1_CaM -> CaMKII_11_11 + PDE1_CaM');
camkiipde_reac_23b_k = addkineticlaw(camkiipde_reac_23b,'MassAction');
set(camkiipde_reac_23b_k, 'ParameterVariableNames', {'koff_camkii_pde'});

camkiipde_reac_23c = addreaction(pfpc_model, 'CaMKII_11_11_PDE1_CaM -> CaMKII_11_11 + PDE1p_CaM');
camkiipde_reac_23c_k = addkineticlaw(camkiipde_reac_23c,'MassAction');
set(camkiipde_reac_23c_k, 'ParameterVariableNames', {'kcat_camkii'});

%CaMKII_11_00
camkiipde_reac_24a = addreaction(pfpc_model, 'CaMKII_11_00 + PDE1_CaM -> CaMKII_11_00_PDE1_CaM');
camkiipde_reac_24a_k = addkineticlaw(camkiipde_reac_24a,'MassAction');
set(camkiipde_reac_24a_k, 'ParameterVariableNames', {'kon_camkii_pde'});

camkiipde_reac_24b = addreaction(pfpc_model, 'CaMKII_11_00_PDE1_CaM -> CaMKII_11_00 + PDE1_CaM');
camkiipde_reac_24b_k = addkineticlaw(camkiipde_reac_24b,'MassAction');
set(camkiipde_reac_24b_k, 'ParameterVariableNames', {'koff_camkii_pde'});

camkiipde_reac_24c = addreaction(pfpc_model, 'CaMKII_11_00_PDE1_CaM -> CaMKII_11_00 + PDE1p_CaM');
camkiipde_reac_24c_k = addkineticlaw(camkiipde_reac_24c,'MassAction');
set(camkiipde_reac_24c_k, 'ParameterVariableNames', {'kcat_camkii'});

%CaMKII_12_00
camkiipde_reac_25a = addreaction(pfpc_model, 'CaMKII_12_00 + PDE1_CaM -> CaMKII_12_00_PDE1_CaM');
camkiipde_reac_25a_k = addkineticlaw(camkiipde_reac_25a,'MassAction');
set(camkiipde_reac_25a_k, 'ParameterVariableNames', {'kon_camkii_pde'});

camkiipde_reac_25b = addreaction(pfpc_model, 'CaMKII_12_00_PDE1_CaM -> CaMKII_12_00 + PDE1_CaM');
camkiipde_reac_25b_k = addkineticlaw(camkiipde_reac_25b,'MassAction');
set(camkiipde_reac_25b_k, 'ParameterVariableNames', {'koff_camkii_pde'});

camkiipde_reac_25c = addreaction(pfpc_model, 'CaMKII_12_00_PDE1_CaM -> CaMKII_12_00 + PDE1p_CaM');
camkiipde_reac_25c_k = addkineticlaw(camkiipde_reac_25c,'MassAction');
set(camkiipde_reac_25c_k, 'ParameterVariableNames', {'kcat_camkii'});

%CaMKII_22_00
camkiipde_reac_26a = addreaction(pfpc_model, 'CaMKII_22_00 + PDE1_CaM -> CaMKII_22_00_PDE1_CaM');
camkiipde_reac_26a_k = addkineticlaw(camkiipde_reac_26a,'MassAction');
set(camkiipde_reac_26a_k, 'ParameterVariableNames', {'kon_camkii_pde'});

camkiipde_reac_26b = addreaction(pfpc_model, 'CaMKII_22_00_PDE1_CaM -> CaMKII_22_00 + PDE1_CaM');
camkiipde_reac_26b_k = addkineticlaw(camkiipde_reac_26b,'MassAction');
set(camkiipde_reac_26b_k, 'ParameterVariableNames', {'koff_camkii_pde'});

camkiipde_reac_26c = addreaction(pfpc_model, 'CaMKII_22_00_PDE1_CaM -> CaMKII_22_00 + PDE1p_CaM');
camkiipde_reac_26c_k = addkineticlaw(camkiipde_reac_26c,'MassAction');
set(camkiipde_reac_26c_k, 'ParameterVariableNames', {'kcat_camkii'});

%CaMKII_11_10
camkiipde_reac_27a = addreaction(pfpc_model, 'CaMKII_11_10 + PDE1_CaM -> CaMKII_11_10_PDE1_CaM');
camkiipde_reac_27a_k = addkineticlaw(camkiipde_reac_27a,'MassAction');
set(camkiipde_reac_27a_k, 'ParameterVariableNames', {'kon_camkii_pde'});

camkiipde_reac_27b = addreaction(pfpc_model, 'CaMKII_11_10_PDE1_CaM -> CaMKII_11_10 + PDE1_CaM');
camkiipde_reac_27b_k = addkineticlaw(camkiipde_reac_27b,'MassAction');
set(camkiipde_reac_27b_k, 'ParameterVariableNames', {'koff_camkii_pde'});

camkiipde_reac_27c = addreaction(pfpc_model, 'CaMKII_11_10_PDE1_CaM -> CaMKII_11_10 + PDE1p_CaM');
camkiipde_reac_27c_k = addkineticlaw(camkiipde_reac_27c,'MassAction');
set(camkiipde_reac_27c_k, 'ParameterVariableNames', {'kcat_camkii'});

%CaMKII_12_10
camkiipde_reac_28a = addreaction(pfpc_model, 'CaMKII_12_10 + PDE1_CaM -> CaMKII_12_10_PDE1_CaM');
camkiipde_reac_28a_k = addkineticlaw(camkiipde_reac_28a,'MassAction');
set(camkiipde_reac_28a_k, 'ParameterVariableNames', {'kon_camkii_pde'});

camkiipde_reac_28b = addreaction(pfpc_model, 'CaMKII_12_10_PDE1_CaM -> CaMKII_12_10 + PDE1_CaM');
camkiipde_reac_28b_k = addkineticlaw(camkiipde_reac_28b,'MassAction');
set(camkiipde_reac_28b_k, 'ParameterVariableNames', {'koff_camkii_pde'});

camkiipde_reac_28c = addreaction(pfpc_model, 'CaMKII_12_10_PDE1_CaM -> CaMKII_12_10 + PDE1p_CaM');
camkiipde_reac_28c_k = addkineticlaw(camkiipde_reac_28c,'MassAction');
set(camkiipde_reac_28c_k, 'ParameterVariableNames', {'kcat_camkii'});

%dephosphorylation of PDE1
camkiipde_reac_29 = addreaction(pfpc_model, 'PDE1p -> PDE1');
camkiipde_reac_29_k = addkineticlaw(camkiipde_reac_29,'MassAction');
set(camkiipde_reac_29_k, 'ParameterVariableNames', {'kdephos_pde'});

%%%%cGMP Production by GC%%%%
%Parameters
kon_gc_no = addparameter(pfpc_model,'kon_gc_no', 'Value', 120, 'ValueUnits', '1/(micromolarity*second)');
koff_gc_no = addparameter(pfpc_model,'koff_gc_no', 'Value', 40, 'ValueUnits', '1/second');
kdeg_no = addparameter(pfpc_model,'kdeg_no', 'Value', 0.0, 'ValueUnits', '1/second');

kon_gc_gtp = addparameter(pfpc_model,'kon_gc_gtp', 'Value', 0.735, 'ValueUnits', '1/(micromolarity*second)');
koff_gc_gtp = addparameter(pfpc_model,'koff_gc_gtp', 'Value', 29.4, 'ValueUnits', '1/second');
kcat_gc = addparameter(pfpc_model,'kcat_gc', 'Value', 7.35, 'ValueUnits', '1/second');

kinflux_no = addparameter(pfpc_model,'kinflux_no', 'Value', 0, 'ValueUnits', 'micromole/second', 'ConstantValue', false);

%Reactions
gc_reac_1 = addreaction(pfpc_model, 'GC + NO -> GC_NO');
gc_reac_1_k = addkineticlaw(gc_reac_1,'MassAction');
set(gc_reac_1_k, 'ParameterVariableNames', {'kon_gc_no'});

gc_reac_2 = addreaction(pfpc_model, 'GC_NO -> GC + NO');
gc_reac_2_k = addkineticlaw(gc_reac_2,'MassAction');
set(gc_reac_2_k, 'ParameterVariableNames', {'koff_gc_no'});

gc_reac_3 = addreaction(pfpc_model, 'GC_NO -> GC');
gc_reac_3_k = addkineticlaw(gc_reac_3,'MassAction');
set(gc_reac_3_k, 'ParameterVariableNames', {'kdeg_no'});

gc_reac_4 = addreaction(pfpc_model, 'GC_NO + GTP -> GC_NO_GTP');
gc_reac_4_k = addkineticlaw(gc_reac_4,'MassAction');
set(gc_reac_4_k, 'ParameterVariableNames', {'kon_gc_gtp'});

gc_reac_5 = addreaction(pfpc_model, 'GC_NO_GTP -> GC_NO + GTP');
gc_reac_5_k = addkineticlaw(gc_reac_5,'MassAction');
set(gc_reac_5_k, 'ParameterVariableNames', {'koff_gc_gtp'});

gc_reac_6 = addreaction(pfpc_model, 'GC_NO_GTP -> GC_NO + cGMP');
gc_reac_6_k = addkineticlaw(gc_reac_6,'MassAction');
set(gc_reac_6_k, 'ParameterVariableNames', {'kcat_gc'});

gc_reac_7 = addreaction(pfpc_model, 'null -> NO');
gc_reac_7_k = addkineticlaw(gc_reac_7,'MassAction');
set(gc_reac_7_k, 'ParameterVariableNames', {'kinflux_no'});

%%%%Hydrolysis of cGMP by PDE1%%%%
%Parameters
kon_pde_cgmp = addparameter(pfpc_model,'kon_pde_cgmp', 'Value', 2.5, 'ValueUnits', '1/(micromolarity*second)');
koff_pde_cgmp = addparameter(pfpc_model,'koff_pde_cgmp', 'Value', 0.4, 'ValueUnits', '1/second');
kcat_pde = addparameter(pfpc_model,'kcat_pde', 'Value', 0.1, 'ValueUnits', '1/second');
kcat_pdecam = addparameter(pfpc_model,'kcat_pdecam', 'Value', 15, 'ValueUnits', '1/second');

%Basal hydrolysis
cgmphydro_reac_1a = addreaction(pfpc_model, 'PDE1 + cGMP -> PDE1_cGMP');
cgmphydro_reac_1a_k = addkineticlaw(cgmphydro_reac_1a,'MassAction');
set(cgmphydro_reac_1a_k, 'ParameterVariableNames', {'kon_pde_cgmp'});

cgmphydro_reac_1b = addreaction(pfpc_model, 'PDE1_cGMP -> PDE1 + cGMP');
cgmphydro_reac_1b_k = addkineticlaw(cgmphydro_reac_1b,'MassAction');
set(cgmphydro_reac_1b_k, 'ParameterVariableNames', {'koff_pde_cgmp'});

cgmphydro_reac_1c = addreaction(pfpc_model, 'PDE1_cGMP -> PDE1 + GMP');
cgmphydro_reac_1c_k = addkineticlaw(cgmphydro_reac_1c,'MassAction');
set(cgmphydro_reac_1c_k, 'ParameterVariableNames', {'kcat_pde'});

cgmphydro_reac_2a = addreaction(pfpc_model, 'PDE1p + cGMP -> PDE1p_cGMP');
cgmphydro_reac_2a_k = addkineticlaw(cgmphydro_reac_2a,'MassAction');
set(cgmphydro_reac_2a_k, 'ParameterVariableNames', {'kon_pde_cgmp'});

cgmphydro_reac_2b = addreaction(pfpc_model, 'PDE1p_cGMP -> PDE1p + cGMP');
cgmphydro_reac_2b_k = addkineticlaw(cgmphydro_reac_2b,'MassAction');
set(cgmphydro_reac_2b_k, 'ParameterVariableNames', {'koff_pde_cgmp'});

cgmphydro_reac_2c = addreaction(pfpc_model, 'PDE1p_cGMP -> PDE1p + GMP');
cgmphydro_reac_2c_k = addkineticlaw(cgmphydro_reac_2c,'MassAction');
set(cgmphydro_reac_2c_k, 'ParameterVariableNames', {'kcat_pde'});


%Active hydrolysis
cgmphydro_reac_3a = addreaction(pfpc_model, 'PDE1_CaM + cGMP -> PDE1_CaM_cGMP');
cgmphydro_reac_3a_k = addkineticlaw(cgmphydro_reac_3a,'MassAction');
set(cgmphydro_reac_3a_k, 'ParameterVariableNames', {'kon_pde_cgmp'});

cgmphydro_reac_3b = addreaction(pfpc_model, 'PDE1_CaM_cGMP -> PDE1_CaM + cGMP');
cgmphydro_reac_3b_k = addkineticlaw(cgmphydro_reac_3b,'MassAction');
set(cgmphydro_reac_3b_k, 'ParameterVariableNames', {'koff_pde_cgmp'});

cgmphydro_reac_3c = addreaction(pfpc_model, 'PDE1_CaM_cGMP -> PDE1_CaM + GMP');
cgmphydro_reac_3c_k = addkineticlaw(cgmphydro_reac_3c,'MassAction');
set(cgmphydro_reac_3c_k, 'ParameterVariableNames', {'kcat_pdecam'});


cgmphydro_reac_4a = addreaction(pfpc_model, 'PDE1p_CaM + cGMP -> PDE1p_CaM_cGMP');
cgmphydro_reac_4a_k = addkineticlaw(cgmphydro_reac_4a,'MassAction');
set(cgmphydro_reac_4a_k, 'ParameterVariableNames', {'kon_pde_cgmp'});

cgmphydro_reac_4b = addreaction(pfpc_model, 'PDE1p_CaM_cGMP -> PDE1p_CaM + cGMP');
cgmphydro_reac_4b_k = addkineticlaw(cgmphydro_reac_4b,'MassAction');
set(cgmphydro_reac_4b_k, 'ParameterVariableNames', {'koff_pde_cgmp'});

cgmphydro_reac_4c = addreaction(pfpc_model, 'PDE1p_CaM_cGMP -> PDE1p_CaM + GMP');
cgmphydro_reac_4c_k = addkineticlaw(cgmphydro_reac_4c,'MassAction');
set(cgmphydro_reac_4c_k, 'ParameterVariableNames', {'kcat_pdecam'});

%%%%%Hydrolysis of cGMP by PDE5%%%%%
kon_pde5_cgmp = addparameter(pfpc_model,'kon_pde5_cgmp', 'Value', 7.5, 'ValueUnits', '1/(micromolarity*second)');
koff_pde5_cgmp = addparameter(pfpc_model,'koff_pde5_cgmp', 'Value', 12, 'ValueUnits', '1/second');
kcat_pde5 = addparameter(pfpc_model,'kcat_pde5', 'Value', 3, 'ValueUnits', '1/second');


cgmphydro_reac_5a = addreaction(pfpc_model, 'PDE5 + cGMP -> PDE5_cGMP');
cgmphydro_reac_5a_k = addkineticlaw(cgmphydro_reac_5a,'MassAction');
set(cgmphydro_reac_5a_k, 'ParameterVariableNames', {'kon_pde5_cgmp'});

cgmphydro_reac_5b = addreaction(pfpc_model, 'PDE5_cGMP -> PDE5 + cGMP');
cgmphydro_reac_5b_k = addkineticlaw(cgmphydro_reac_5b,'MassAction');
set(cgmphydro_reac_5b_k, 'ParameterVariableNames', {'koff_pde5_cgmp'});

cgmphydro_reac_5c = addreaction(pfpc_model, 'PDE5_cGMP -> PDE5 + GMP');
cgmphydro_reac_5c_k = addkineticlaw(cgmphydro_reac_5c,'MassAction');
set(cgmphydro_reac_5c_k, 'ParameterVariableNames', {'kcat_pde5'});


%%%%Activation of PKG by cGMP%%%%
%Parameters
kon_pkg_cgmp = addparameter(pfpc_model,'kon_pkg_cgmp', 'Value', 0.5, 'ValueUnits', '1/(micromolarity*second)');
koff_pkg_cgmp = addparameter(pfpc_model,'koff_pkg_cgmp', 'Value', 0.65, 'ValueUnits', '1/second');

kon_pkg_gsub = addparameter(pfpc_model,'kon_pkg_gsub', 'Value', 0.5, 'ValueUnits', '1/(micromolarity*second)');
koff_pkg_gsub = addparameter(pfpc_model,'koff_pkg_gsub', 'Value', 8, 'ValueUnits', '1/second');
kcat_pkg = addparameter(pfpc_model,'kcat_pkg', 'Value', 0.5, 'ValueUnits', '1/second');

%Activation reactions
pkg_reac_1 = addreaction(pfpc_model, 'PKG + cGMP -> PKG_cGMP');
pkg_reac_1_k = addkineticlaw(pkg_reac_1,'MassAction');
set(pkg_reac_1_k, 'ParameterVariableNames', {'kon_pkg_cgmp'});

pkg_reac_2 = addreaction(pfpc_model, 'PKG_cGMP -> PKG + cGMP');
pkg_reac_2_k = addkineticlaw(pkg_reac_2,'MassAction');
set(pkg_reac_2_k, 'ParameterVariableNames', {'koff_pkg_cgmp'});

%Phosphorylation of Gsubstrate by PKG
pkg_reac_3a = addreaction(pfpc_model, 'PKG_cGMP + Gsub -> PKG_cGMP_Gsub');
pkg_reac_3a_k = addkineticlaw(pkg_reac_3a,'MassAction');
set(pkg_reac_3a_k, 'ParameterVariableNames', {'kon_pkg_gsub'});

pkg_reac_3b = addreaction(pfpc_model, 'PKG_cGMP_Gsub -> PKG_cGMP + Gsub');
pkg_reac_3b_k = addkineticlaw(pkg_reac_3b,'MassAction');
set(pkg_reac_3b_k, 'ParameterVariableNames', {'koff_pkg_gsub'});

pkg_reac_3c = addreaction(pfpc_model, 'PKG_cGMP_Gsub -> PKG_cGMP + pGsub');
pkg_reac_3c_k = addkineticlaw(pkg_reac_3c,'MassAction');
set(pkg_reac_3c_k, 'ParameterVariableNames', {'kcat_pkg'});

%%%%Reactions of Gsubstrate
%Parameters
kon_pp2a_gsub_inh = addparameter(pfpc_model,'kon_pp2a_gsub_inh', 'Value', 10, 'ValueUnits', '1/(micromolarity*second)');
koff_pp2a_gsub_inh = addparameter(pfpc_model,'koff_pp2a_gsub_inh', 'Value', 0.06, 'ValueUnits', '1/second');

kon_pp1_gsub_inh = addparameter(pfpc_model,'kon_pp1_gsub_inh', 'Value', 10, 'ValueUnits', '1/(micromolarity*second)');
koff_pp1_gsub_inh = addparameter(pfpc_model,'koff_pp1_gsub_inh', 'Value', 0.06, 'ValueUnits', '1/second');

kon_pp2a_gsub_dephos = addparameter(pfpc_model,'kon_pp2a_gsub_dephos', 'Value', 2, 'ValueUnits', '1/(micromolarity*second)');
koff_pp2a_gsub_dephos = addparameter(pfpc_model,'koff_pp2a_gsub_dephos', 'Value', 8, 'ValueUnits', '1/second');
kcat_pp2a_gsub = addparameter(pfpc_model,'kcat_pp2a_gsub', 'Value', 2, 'ValueUnits', '1/second');

%Inhibition of PP2A by pGsub
gsub_reac_1 = addreaction(pfpc_model, 'PP2A + pGsub -> PP2A_pGsub');
gsub_reac_1_k = addkineticlaw(gsub_reac_1,'MassAction');
set(gsub_reac_1_k, 'ParameterVariableNames', {'kon_pp2a_gsub_inh'});

gsub_reac_2 = addreaction(pfpc_model, 'PP2A_pGsub -> PP2A + pGsub');
gsub_reac_2_k = addkineticlaw(gsub_reac_2,'MassAction');
set(gsub_reac_2_k, 'ParameterVariableNames', {'koff_pp2a_gsub_inh'});

%Inhibition of PP1 by pGsub
gsub_reac_1a = addreaction(pfpc_model, 'PP1 + pGsub -> PP1_pGsub');
gsub_reac_1a_k = addkineticlaw(gsub_reac_1a,'MassAction');
set(gsub_reac_1a_k, 'ParameterVariableNames', {'kon_pp1_gsub_inh'});

gsub_reac_2a = addreaction(pfpc_model, 'PP1_pGsub -> PP1 + pGsub');
gsub_reac_2a_k = addkineticlaw(gsub_reac_2a,'MassAction');
set(gsub_reac_2a_k, 'ParameterVariableNames', {'koff_pp1_gsub_inh'});

%Dephosphorylation of Gsubstrate by PP2A
gsub_reac_3a = addreaction(pfpc_model, 'pGsub + PP2A -> pGsub_PP2A');
gsub_reac_3a_k = addkineticlaw(gsub_reac_3a,'MassAction');
set(gsub_reac_3a_k, 'ParameterVariableNames', {'kon_pp2a_gsub_dephos'});

gsub_reac_3b = addreaction(pfpc_model, 'pGsub_PP2A -> pGsub + PP2A');
gsub_reac_3b_k = addkineticlaw(gsub_reac_3b,'MassAction');
set(gsub_reac_3b_k, 'ParameterVariableNames', {'koff_pp2a_gsub_dephos'});

gsub_reac_3c = addreaction(pfpc_model, 'pGsub_PP2A -> Gsub + PP2A');
gsub_reac_3c_k = addkineticlaw(gsub_reac_3c,'MassAction');
set(gsub_reac_3c_k, 'ParameterVariableNames', {'kcat_pp2a_gsub'});

%%%%Calcineurin Reactions%%%%
%Parameters
kon_can_cam = addparameter(pfpc_model,'kon_can_cam', 'Value', 100, 'ValueUnits', '1/(micromolarity*second)', 'ConstantValue', false);
koff_can_cam = addparameter(pfpc_model,'koff_can_cam', 'Value', 3, 'ValueUnits', '1/second');

kon_can_pgsub = addparameter(pfpc_model,'kon_can_pgsub', 'Value', 1.2, 'ValueUnits', '1/(micromolarity*second)', 'ConstantValue', false);
koff_can_pgsub = addparameter(pfpc_model,'koff_can_pgsub', 'Value', 2, 'ValueUnits', '1/second');
kcat_can = addparameter(pfpc_model,'kcat_can', 'Value', 0.5, 'ValueUnits', '1/second');

%Calcineurin binding to CaM
can_reac_1 = addreaction(pfpc_model, 'CaN + Ca4CaM -> CaN_CaM');
can_reac_1_k = addkineticlaw(can_reac_1,'MassAction');
set(can_reac_1_k, 'ParameterVariableNames', {'kon_can_cam'});

can_reac_2 = addreaction(pfpc_model, 'CaN_CaM -> CaN + Ca4CaM');
can_reac_2_k = addkineticlaw(can_reac_2,'MassAction');
set(can_reac_2_k, 'ParameterVariableNames', {'koff_can_cam'});

%Dephosphorylation of Gsubstrate by CaN
can_reac_3a = addreaction(pfpc_model, 'CaN_CaM + pGsub -> CaN_pGsub');
can_reac_3a_k = addkineticlaw(can_reac_3a,'MassAction');
set(can_reac_3a_k, 'ParameterVariableNames', {'kon_can_pgsub'});

can_reac_3b = addreaction(pfpc_model, 'CaN_pGsub -> CaN_CaM + pGsub');
can_reac_3b_k = addkineticlaw(can_reac_3b,'MassAction');
set(can_reac_3b_k, 'ParameterVariableNames', {'koff_can_pgsub'});

can_reac_3c = addreaction(pfpc_model, 'CaN_pGsub -> CaN_CaM + Gsub');
can_reac_3c_k = addkineticlaw(can_reac_3c,'MassAction');
set(can_reac_3c_k, 'ParameterVariableNames', {'kcat_can'});

%%%%%%LOOP MODEL%%%%%%%%%
%Ras-RAF-MEK-ERK
ras_gdp = addspecies(cytosol, 'ras_gdp', 'InitialAmount', 0.4, 'InitialAmountUnits', 'micromolarity');
ras_gtp = addspecies(cytosol, 'ras_gtp', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
RAF = addspecies(cytosol, 'RAF', 'InitialAmount', 0.013, 'InitialAmountUnits', 'micromolarity');
MEK = addspecies(cytosol, 'MEK', 'InitialAmount', 1.4, 'InitialAmountUnits', 'micromolarity');
ERK = addspecies(cytosol, 'ERK', 'InitialAmount', 0.96, 'InitialAmountUnits', 'micromolarity');
PP5 = addspecies(cytosol, 'PP5', 'InitialAmount', 0.1, 'InitialAmountUnits', 'micromolarity');
MKP = addspecies(cytosol, 'MKP', 'InitialAmount', 0.05, 'InitialAmountUnits', 'micromolarity');
RAFK = addspecies(cytosol, 'RAFK', 'InitialAmount', 0.1, 'InitialAmountUnits', 'micromolarity');
ras_gdp_pGEF = addspecies(cytosol, 'ras_gdp_pGEF', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
ras_gtp_pGAP = addspecies(cytosol, 'ras_gtp_pGAP', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
ras_RAF = addspecies(cytosol, 'ras_RAF', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
ras_RAF_act = addspecies(cytosol, 'ras_RAF_act', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
ras_RAF_RAFK = addspecies(cytosol, 'ras_RAF_RAFK', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
ras_pRAF_act = addspecies(cytosol, 'ras_pRAF_act', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
pRAF_act = addspecies(cytosol, 'pRAF_act', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
pRAF_PP5 = addspecies(cytosol, 'pRAF_PP5', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
RAF_RKIP = addspecies(cytosol, 'RAF_RKIP', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
RAF_RKIP_PKC = addspecies(cytosol, 'RAF_RKIP_PKC', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
pRKIP = addspecies(cytosol, 'pRKIP', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
pRKIP_RP = addspecies(cytosol, 'pRKIP_RP', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
ras_RAF_MEK = addspecies(cytosol, 'ras_RAF_MEK', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
pMEK = addspecies(cytosol, 'pMEK', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
ras_RAF_pMEK = addspecies(cytosol, 'ras_RAF_pMEK', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
ppMEK = addspecies(cytosol, 'ppMEK', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
ras_pRAF_MEK = addspecies(cytosol, 'ras_pRAF_MEK', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
ras_pRAF_pMEK = addspecies(cytosol, 'ras_pRAF_pMEK', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
pRAF_MEK = addspecies(cytosol, 'pRAF_MEK', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
pRAF_pMEK = addspecies(cytosol, 'pRAF_pMEK', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
ppMEK_ERK = addspecies(cytosol, 'ppMEK_ERK', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
pERK = addspecies(cytosol, 'pERK', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
ppMEK_pERK = addspecies(cytosol, 'ppMEK_pERK', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
ppERK = addspecies(cytosol, 'ppERK', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
ppMEK_PP2A = addspecies(cytosol, 'ppMEK_PP2A', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
pMEK_PP2A = addspecies(cytosol, 'pMEK_PP2A', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
ppERK_MKP = addspecies(cytosol, 'ppERK_MKP', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
pERK_MKP = addspecies(cytosol, 'pERK_MKP', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
Ras_act = addspecies(cytosol, 'Ras_act', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
RAF_act = addspecies(cytosol, 'RAF_act', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
ERK_act = addspecies(cytosol, 'ERK_act', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
MEK_phospho = addspecies(cytosol, 'MEK_phospho', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');

%PLA2-AA
PLA2c  = addspecies(cytosol, 'PLA2c ', 'InitialAmount', 1, 'InitialAmountUnits', 'micromolarity');
AA = addspecies(cytosol, 'AA', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
PLA2c_ppERK = addspecies(cytosol, 'PLA2c_ppERK', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
pPLA2c = addspecies(cytosol, 'pPLA2c', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
CaPLA2c_ppERK = addspecies(cytosol, 'CaPLA2c_ppERK', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
CapPLA2c = addspecies(cytosol, 'CapPLA2c', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
Ca2PLA2c_ppERK = addspecies(cytosol, 'Ca2PLA2c_ppERK', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
Ca2pPLA2c = addspecies(cytosol, 'Ca2pPLA2c', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
pPLA2c_PP2A = addspecies(cytosol, 'pPLA2c_PP2A', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
CapPLA2c_PP2A = addspecies(cytosol, 'CapPLA2c_PP2A', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
Ca2pPLA2c_PP2A = addspecies(cytosol, 'Ca2pPLA2c_PP2A', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
pPLA2c_PP1 = addspecies(cytosol, 'pPLA2c_PP1', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
CapPLA2c_PP1 = addspecies(cytosol, 'CapPLA2c_PP1', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
Ca2pPLA2c_PP1 = addspecies(cytosol, 'Ca2pPLA2c_PP1', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
Ca2pPLA2 = addspecies(cytosol, 'Ca2pPLA2', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
Ca2pPLA2_act = addspecies(cytosol, 'Ca2pPLA2_act', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
CapPLA2 = addspecies(cytosol, 'CapPLA2', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
DAG_Ca2pPLA2_act = addspecies(cytosol, 'DAG_Ca2pPLA2_act', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
CaPLA2c = addspecies(cytosol, 'CaPLA2c', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
Ca2PLA2c = addspecies(cytosol, 'Ca2PLA2c', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
Ca2PLA2 = addspecies(cytosol, 'Ca2PLA2', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
Ca2PLA2_act = addspecies(cytosol, 'Ca2PLA2_act', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
CaPLA2 = addspecies(cytosol, 'CaPLA2', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
DAG_Ca2PLA2_act = addspecies(cytosol, 'DAG_Ca2PLA2_act', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
Ca2PLA2_act_APC = addspecies(cytosol, 'Ca2PLA2_act_APC', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
DAG_Ca2PLA2_act_APC = addspecies(cytosol, 'DAG_Ca2PLA2_act_APC', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
Ca2pPLA2_act_APC = addspecies(cytosol, 'Ca2pPLA2_act_APC', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
DAG_Ca2pPLA2_act_APC = addspecies(cytosol, 'DAG_Ca2pPLA2_act_APC', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
PLA2_active = addspecies(cytosol, 'PLA2_active', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');

%PKC
GEF = addspecies(cytosol, 'GEF', 'InitialAmount', 0.1, 'InitialAmountUnits', 'micromolarity');
pGAP = addspecies(cytosol, 'pGAP', 'InitialAmount', 0.002, 'InitialAmountUnits', 'micromolarity');
GAP = addspecies(cytosol, 'GAP', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
PKCc = addspecies(cytosol, 'PKCc', 'InitialAmount', 1, 'InitialAmountUnits', 'micromolarity');
CaPKC_act = addspecies(cytosol, 'CaPKC_act', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
GEF_CaPKC_act = addspecies(cytosol, 'GEF_CaPKC_act', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
pGEF = addspecies(cytosol, 'pGEF', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
Ca2PKC_act = addspecies(cytosol, 'Ca2PKC_act', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
GEF_Ca2PKC_act = addspecies(cytosol, 'GEF_Ca2PKC_act', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
CaPKC_AA_act = addspecies(cytosol, 'CaPKC_AA_act', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
GEF_CaPKC_AA_act = addspecies(cytosol, 'GEF_CaPKC_AA_act', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
Ca2PKC_AA_act = addspecies(cytosol, 'Ca2PKC_AA_act', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
GEF_Ca2PKC_AA_act = addspecies(cytosol, 'GEF_Ca2PKC_AA_act', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
DAG_CaPKC_act = addspecies(cytosol, 'DAG_CaPKC_act', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
GEF_DAG_CaPKC_act = addspecies(cytosol, 'GEF_DAG_CaPKC_act', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
DAG_Ca2PKC_act = addspecies(cytosol, 'DAG_Ca2PKC_act', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
GEF_DAG_Ca2PKC_act = addspecies(cytosol, 'GEF_DAG_Ca2PKC_act', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
DAG_CaPKC_AA_act = addspecies(cytosol, 'DAG_CaPKC_AA_act', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
GEF_DAG_CaPKC_AA_act = addspecies(cytosol, 'GEF_DAG_CaPKC_AA_act', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
DAG_Ca2PKC_AA_act = addspecies(cytosol, 'DAG_Ca2PKC_AA_act', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
GEF_DAG_Ca2PKC_AA_act = addspecies(cytosol, 'GEF_DAG_Ca2PKC_AA_act', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
CaPKCc = addspecies(cytosol, 'CaPKCc', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
CaPKC = addspecies(cytosol, 'CaPKC', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
Ca2PKC = addspecies(cytosol, 'Ca2PKC', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
CaPKC_AA = addspecies(cytosol, 'CaPKC_AA', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
Ca2PKC_AA = addspecies(cytosol, 'Ca2PKC_AA', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
RKIP = addspecies(cytosol, 'RKIP', 'InitialAmount', 0.2, 'InitialAmountUnits', 'micromolarity');
RP = addspecies(cytosol, 'RP', 'InitialAmount', 0.01, 'InitialAmountUnits', 'micromolarity');
RAF_RKIP_CaPKC_act = addspecies(cytosol, 'RAF_RKIP_CaPKC_act', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
RAF_RKIP_Ca2PKC_act = addspecies(cytosol, 'RAF_RKIP_Ca2PKC_act', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
RAF_RKIP_CaPKC_AA_act = addspecies(cytosol, 'RAF_RKIP_CaPKC_AA_act', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
RAF_RKIP_Ca2PKC_AA_act = addspecies(cytosol, 'RAF_RKIP_Ca2PKC_AA_act', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
RAF_RKIP_DAG_CaPKC_act = addspecies(cytosol, 'RAF_RKIP_DAG_CaPKC_act', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
RAF_RKIP_DAG_Ca2PKC_act = addspecies(cytosol, 'RAF_RKIP_DAG_Ca2PKC_act', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
RAF_RKIP_DAG_CaPKC_AA_act = addspecies(cytosol, 'RAF_RKIP_DAG_CaPKC_AA_act', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
RAF_RKIP_DAG_Ca2PKC_AA_act = addspecies(cytosol, 'RAF_RKIP_DAG_Ca2PKC_AA_act', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
PKC_active = addspecies(cytosol, 'PKC_active', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');

DGK = addspecies(cytosol, 'DGK', 'InitialAmount', 0.5, 'InitialAmountUnits', 'micromolarity');
pDGK = addspecies(cytosol, 'pDGK', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
PP2A_pDGK = addspecies(cytosol, 'PP2A_pDGK', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');

CaPKC_act_DGK = addspecies(cytosol, 'CaPKC_act_DGK', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
Ca2PKC_act_DGK = addspecies(cytosol, 'Ca2PKC_act_DGK', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
CaPKC_AA_act_DGK = addspecies(cytosol, 'CaPKC_AA_act_DGK', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
Ca2PKC_AA_act_DGK = addspecies(cytosol, 'Ca2PKC_AA_act_DGK', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
DAG_CaPKC_act_DGK = addspecies(cytosol, 'DAG_CaPKC_act_DGK', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
DAG_Ca2PKC_act_DGK = addspecies(cytosol, 'DAG_Ca2PKC_act_DGK', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
DAG_CaPKC_AA_act_DGK = addspecies(cytosol, 'DAG_CaPKC_AA_act_DGK', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
DAG_Ca2PKC_AA_act_DGK = addspecies(cytosol, 'DAG_Ca2PKC_AA_act_DGK', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
DGK_DAG = addspecies(cytosol, 'DGK_DAG', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');

%Parameters
kon_raf_mek = addparameter(pfpc_model,'kon_raf_mek', 'Value', 0.65, 'ValueUnits', '1/(micromolarity*second)');
kon_ras_raf = addparameter(pfpc_model,'kon_ras_raf', 'Value', 0.49, 'ValueUnits', '1/(micromolarity*second)');
koff_ras_raf = addparameter(pfpc_model,'koff_ras_raf', 'Value', 0.049, 'ValueUnits', '1/second');
koff_raf_mek = addparameter(pfpc_model,'koff_raf_mek', 'Value', 0.065, 'ValueUnits', '1/second');
kcat_raf = addparameter(pfpc_model,'kcat_raf', 'Value', 0.18, 'ValueUnits', '1/second');
kon_mek_erk = addparameter(pfpc_model,'kon_mek_erk', 'Value', 12.8, 'ValueUnits', '1/(micromolarity*second)');
koff_mek_erk = addparameter(pfpc_model,'koff_mek_erk', 'Value', 0.88, 'ValueUnits', '1/second');
kcat_mek = addparameter(pfpc_model,'kcat_mek', 'Value', 0.32, 'ValueUnits', '1/second');
kon_mek_pp2a = addparameter(pfpc_model,'kon_mek_pp2a', 'Value', 0.38, 'ValueUnits', '1/(micromolarity*second)', 'ConstantValue', false);
koff_mek_pp2a = addparameter(pfpc_model,'koff_mek_pp2a', 'Value', 2, 'ValueUnits', '1/second');

kcat_pp2a = addparameter(pfpc_model,'kcat_pp2a', 'Value', 0.8, 'ValueUnits', '1/second', 'ConstantValue', false);
kon_erk_mkp = addparameter(pfpc_model,'kon_erk_mkp', 'Value', 2, 'ValueUnits', '1/(micromolarity*second)');
koff_erk_mkp = addparameter(pfpc_model,'koff_erk_mkp', 'Value', 0.5, 'ValueUnits', '1/second');
kcat_mkp = addparameter(pfpc_model,'kcat_mkp', 'Value', 0.5, 'ValueUnits', '1/second');
kras_deact = addparameter(pfpc_model,'kras_deact', 'Value', 0, 'ValueUnits', '1/second');
kon_gef_pkc = addparameter(pfpc_model,'kon_gef_pkc', 'Value', 6, 'ValueUnits', '1/(micromolarity*second)');
kcat_pkc = addparameter(pfpc_model,'kcat_pkc', 'Value', 4, 'ValueUnits', '1/second');
koff_gef_pkc = addparameter(pfpc_model,'koff_gef_pkc', 'Value', 16, 'ValueUnits', '1/second');
k_gef_deact = addparameter(pfpc_model,'k_gef_deact', 'Value', 1, 'ValueUnits', '1/second');
kon_ras_gef = addparameter(pfpc_model,'kon_ras_gef', 'Value', 2, 'ValueUnits', '1/(micromolarity*second)');

koff_ras_gef = addparameter(pfpc_model,'koff_ras_gef', 'Value', 0.08, 'ValueUnits', '1/second');
kcat_gef = addparameter(pfpc_model,'kcat_gef', 'Value', 0.5, 'ValueUnits', '1/second');
kon_ras_gap = addparameter(pfpc_model,'kon_ras_gap', 'Value', 50, 'ValueUnits', '1/(micromolarity*second)');
kcat_gap = addparameter(pfpc_model,'kcat_gap', 'Value', 1.8, 'ValueUnits', '1/second');
koff_ras_gap = addparameter(pfpc_model,'koff_ras_gap', 'Value', 40, 'ValueUnits', '1/second');
%k_pkc_deact = addparameter(pfpc_model,'k_pkc_deact', 'Value', 0, 'ValueUnits', '1/second');
k_raf_act = addparameter(pfpc_model,'k_raf_act', 'Value', 1, 'ValueUnits', '1/second');
k_raf_deact = addparameter(pfpc_model,'k_raf_deact', 'Value', 3, 'ValueUnits', '1/second');
kon_raf_rkip = addparameter(pfpc_model,'kon_raf_rkip', 'Value', 3.5, 'ValueUnits', '1/(micromolarity*second)');
koff_raf_rkip = addparameter(pfpc_model,'koff_raf_rkip', 'Value', 0.0072, 'ValueUnits', '1/second');

kon_rkip_pkc = addparameter(pfpc_model,'kon_rkip_pkc', 'Value', 8, 'ValueUnits', '1/(micromolarity*second)');
koff_rkip_pkc = addparameter(pfpc_model,'koff_rkip_pkc', 'Value', 2, 'ValueUnits', '1/second');
kon_raf_rafk = addparameter(pfpc_model,'kon_raf_rafk', 'Value', 8, 'ValueUnits', '1/(micromolarity*second)');
koff_raf_rafk = addparameter(pfpc_model,'koff_raf_rafk', 'Value', 3.5, 'ValueUnits', '1/second');
kcat_rafk = addparameter(pfpc_model,'kcat_rafk', 'Value', 0.9, 'ValueUnits', '1/second');
kon_pla2_ca1 = addparameter(pfpc_model,'kon_pla2_ca1', 'Value', 17.7199, 'ValueUnits', '1/(micromolarity*second)');
koff_pla2_ca1 = addparameter(pfpc_model,'koff_pla2_ca1', 'Value', 154.563, 'ValueUnits', '1/second');
kon_pla2_ca2 = addparameter(pfpc_model,'kon_pla2_ca2', 'Value', 16.3127, 'ValueUnits', '1/(micromolarity*second)');
koff_pla2_ca2 = addparameter(pfpc_model,'koff_pla2_ca2', 'Value', 101.625, 'ValueUnits', '1/second');
kon_pla2_pm = addparameter(pfpc_model,'kon_pla2_pm', 'Value', 26, 'ValueUnits', '1/second');

kon_pla2_dag = addparameter(pfpc_model,'kon_pla2_dag', 'Value', 2.9, 'ValueUnits', '1/(micromolarity*second)');
koff_pla2_dag = addparameter(pfpc_model,'koff_pla2_dag', 'Value', 3.3, 'ValueUnits', '1/second');
koff1_pla2_pm = addparameter(pfpc_model,'koff1_pla2_pm', 'Value', 0.0035, 'ValueUnits', '1/second');
kact_pla2 = addparameter(pfpc_model,'kact_pla2', 'Value', 1.5, 'ValueUnits', '1/second');
kdeact_pla2 = addparameter(pfpc_model,'kdeact_pla2', 'Value', 8, 'ValueUnits', '1/second');
kon_rkip_rp = addparameter(pfpc_model,'kon_rkip_rp', 'Value', 4.6, 'ValueUnits', '1/(micromolarity*second)');
koff_rkip_rp = addparameter(pfpc_model,'koff_rkip_rp', 'Value', 3.5, 'ValueUnits', '1/second');
kcat_rp = addparameter(pfpc_model,'kcat_rp', 'Value', 0.87, 'ValueUnits', '1/second');
kon_raf_pp5 = addparameter(pfpc_model,'kon_raf_pp5', 'Value', 1, 'ValueUnits', '1/(micromolarity*second)');
koff_raf_pp5 = addparameter(pfpc_model,'koff_raf_pp5', 'Value', 1, 'ValueUnits', '1/second');

kcat_pp5 = addparameter(pfpc_model,'kcat_pp5', 'Value', 1, 'ValueUnits', '1/second');
kdeact_pla2_pS505 = addparameter(pfpc_model,'kdeact_pla2_pS505', 'Value', 8, 'ValueUnits', '1/second');
kon_pla2_erk = addparameter(pfpc_model,'kon_pla2_erk', 'Value', 1, 'ValueUnits', '1/(micromolarity*second)');
koff_pla2_erk = addparameter(pfpc_model,'koff_pla2_erk', 'Value', 4, 'ValueUnits', '1/second');
kcat_erk = addparameter(pfpc_model,'kcat_erk', 'Value', 14, 'ValueUnits', '1/second');
kon_pla2_pp2a = addparameter(pfpc_model,'kon_pla2_pp2a', 'Value', 1, 'ValueUnits', '1/(micromolarity*second)', 'ConstantValue', false);
koff_pla2_pp2a = addparameter(pfpc_model,'koff_pla2_pp2a', 'Value', 1.4, 'ValueUnits', '1/second');
kon_pla2_apc = addparameter(pfpc_model,'kon_pla2_apc', 'Value', 30, 'ValueUnits', '1/second');
koff_pla2_apc = addparameter(pfpc_model,'koff_pla2_apc', 'Value', 20, 'ValueUnits', '1/second');
kcat_pla2 = addparameter(pfpc_model,'kcat_pla2', 'Value', 30, 'ValueUnits', '1/second');
kdeg_aa = addparameter(pfpc_model,'kdeg_aa', 'Value', 15, 'ValueUnits', '1/second');

kon_pkc_ca1 = addparameter(pfpc_model,'kon_pkc_ca1', 'Value', 13.3, 'ValueUnits', '1/(micromolarity*second)');
koff_pkc_ca1 = addparameter(pfpc_model,'koff_pkc_ca1', 'Value', 12, 'ValueUnits', '1/second');
kon_pkc_pm = addparameter(pfpc_model,'kon_pkc_pm', 'Value', 420, 'ValueUnits', '1/second');
koff_pkc_pm = addparameter(pfpc_model,'koff_pkc_pm', 'Value', 17.3, 'ValueUnits', '1/second');
kon_pkc_ca2 = addparameter(pfpc_model,'kon_pkc_ca2', 'Value', 1, 'ValueUnits', '1/(micromolarity*second)');
kon_pkc_ca2_aa = addparameter(pfpc_model,'kon_pkc_ca2_aa', 'Value', 20, 'ValueUnits', '1/(micromolarity*second)');
koff_pkc_ca2 = addparameter(pfpc_model,'koff_pkc_ca2', 'Value', 12, 'ValueUnits', '1/second');
kon_pkc_aa = addparameter(pfpc_model,'kon_pkc_aa', 'Value', 4, 'ValueUnits', '1/(micromolarity*second)');
koff_pkc_aa = addparameter(pfpc_model,'koff_pkc_aa', 'Value', 1, 'ValueUnits', '1/second');
kon_pkc_dag = addparameter(pfpc_model,'kon_pkc_dag', 'Value', 0.8, 'ValueUnits', '1/(micromolarity*second)');

koff_pkc_dag = addparameter(pfpc_model,'koff_pkc_dag', 'Value', 1, 'ValueUnits', '1/second');
kon_pkc_gef = addparameter(pfpc_model,'kon_pkc_gef', 'Value', 6, 'ValueUnits', '1/(micromolarity*second)');
koff_pkc_gef = addparameter(pfpc_model,'koff_pkc_gef', 'Value', 16, 'ValueUnits', '1/second');
kcat_pkc_gef = addparameter(pfpc_model,'kcat_pkc_gef', 'Value', 0.1, 'ValueUnits', '1/second');
kon_pla2_pp1 = addparameter(pfpc_model,'kon_pla2_pp1', 'Value', 1.4, 'ValueUnits', '1/(micromolarity*second)', 'ConstantValue', false);
koff_pla2_pp1 = addparameter(pfpc_model,'koff_pla2_pp1', 'Value', 1.5, 'ValueUnits', '1/second');
kcat_pp1 = addparameter(pfpc_model,'kcat_pp1', 'Value', 2.5, 'ValueUnits', '1/second');
kact_pla2_pS505 = addparameter(pfpc_model,'kact_pla2_pS505', 'Value', 29, 'ValueUnits', '1/second');
kact_pkc = addparameter(pfpc_model,'kact_pkc', 'Value', 3, 'ValueUnits', '1/second');
kdeact_pkc = addparameter(pfpc_model,'kdeact_pkc', 'Value', 5, 'ValueUnits', '1/second');
kdeact_pkc_aa = addparameter(pfpc_model,'kdeact_pkc_aa', 'Value', 5, 'ValueUnits', '1/second');


%same as for PKC-GEF (apart from kcat)
kon_pkc_dgk = addparameter(pfpc_model,'kon_pkc_dgk', 'Value', 6, 'ValueUnits', '1/(micromolarity*second)');
koff_pkc_dgk = addparameter(pfpc_model,'koff_pkc_dgk', 'Value', 16, 'ValueUnits', '1/second');
kcat_pkc_dgk = addparameter(pfpc_model,'kcat_pkc_dgk', 'Value', 4, 'ValueUnits', '1/second');
%assumed
kon_dgk_dag = addparameter(pfpc_model,'kon_dgk_dag', 'Value', 10, 'ValueUnits', '1/(micromolarity*second)');
koff_dgk_dag = addparameter(pfpc_model,'koff_dgk_dag', 'Value', 0.5, 'ValueUnits', '1/second');
kcat_dgk = addparameter(pfpc_model,'kcat_dgk', 'Value', 2, 'ValueUnits', '1/second');
%same as for Gsubstrate
kon_pp2a_dgk = addparameter(pfpc_model,'kon_pp2a_dgk', 'Value', 1, 'ValueUnits', '1/(micromolarity*second)');
koff_pp2a_dgk = addparameter(pfpc_model,'koff_pp2a_dgk', 'Value', 8, 'ValueUnits', '1/second');
kcat_pp2a_dgk = addparameter(pfpc_model,'kcat_pp2a_dgk', 'Value', 0.4, 'ValueUnits', '1/second');

%%%%%%Ras-RAF-MEK-ERK Reactions%%%%%%%%%%

rrme_reac_1 = addreaction(pfpc_model, 'ras_gtp + RAF -> ras_RAF');
rrme_reac_1_k = addkineticlaw(rrme_reac_1,'MassAction');
set(rrme_reac_1_k, 'ParameterVariableNames', {'kon_ras_raf'});

rrme_reac_2 = addreaction(pfpc_model, 'ras_RAF -> ras_gtp + RAF');
rrme_reac_2_k = addkineticlaw(rrme_reac_2,'MassAction');
set(rrme_reac_2_k, 'ParameterVariableNames', {'koff_ras_raf'});

rrme_reac_3 = addreaction(pfpc_model, 'ras_RAF_act + MEK -> ras_RAF_MEK');
rrme_reac_3_k = addkineticlaw(rrme_reac_3,'MassAction');
set(rrme_reac_3_k, 'ParameterVariableNames', {'kon_raf_mek'});

rrme_reac_4 = addreaction(pfpc_model, 'ras_RAF_MEK -> ras_RAF_act + MEK');
rrme_reac_4_k = addkineticlaw(rrme_reac_4,'MassAction');
set(rrme_reac_4_k, 'ParameterVariableNames', {'koff_raf_mek'});

rrme_reac_5 = addreaction(pfpc_model, 'ras_RAF_MEK -> ras_RAF_act + pMEK');
rrme_reac_5_k = addkineticlaw(rrme_reac_5,'MassAction');
set(rrme_reac_5_k, 'ParameterVariableNames', {'kcat_raf'});

rrme_reac_6 = addreaction(pfpc_model, 'ras_RAF_act + pMEK -> ras_RAF_pMEK');
rrme_reac_6_k = addkineticlaw(rrme_reac_6,'MassAction');
set(rrme_reac_6_k, 'ParameterVariableNames', {'kon_raf_mek'});

rrme_reac_7 = addreaction(pfpc_model, 'ras_RAF_pMEK -> ras_RAF_act + pMEK');
rrme_reac_7_k = addkineticlaw(rrme_reac_7,'MassAction');
set(rrme_reac_7_k, 'ParameterVariableNames', {'koff_raf_mek'});

rrme_reac_8 = addreaction(pfpc_model, 'ras_RAF_pMEK -> ras_RAF_act + ppMEK');
rrme_reac_8_k = addkineticlaw(rrme_reac_8,'MassAction');
set(rrme_reac_8_k, 'ParameterVariableNames', {'kcat_raf'});

rrme_reac_9 = addreaction(pfpc_model, 'ppMEK + ERK -> ppMEK_ERK');
rrme_reac_9_k = addkineticlaw(rrme_reac_9,'MassAction');
set(rrme_reac_9_k, 'ParameterVariableNames', {'kon_mek_erk'});

rrme_reac_10 = addreaction(pfpc_model, 'ppMEK_ERK -> ppMEK + ERK');
rrme_reac_10_k = addkineticlaw(rrme_reac_10,'MassAction');
set(rrme_reac_10_k, 'ParameterVariableNames', {'koff_mek_erk'});

rrme_reac_11 = addreaction(pfpc_model, 'ppMEK_ERK -> ppMEK + pERK');
rrme_reac_11_k = addkineticlaw(rrme_reac_11,'MassAction');
set(rrme_reac_11_k, 'ParameterVariableNames', {'kcat_mek'});

rrme_reac_12 = addreaction(pfpc_model, 'ppMEK + pERK -> ppMEK_pERK');
rrme_reac_12_k = addkineticlaw(rrme_reac_12,'MassAction');
set(rrme_reac_12_k, 'ParameterVariableNames', {'kon_mek_erk'});

rrme_reac_13 = addreaction(pfpc_model, 'ppMEK_pERK -> ppMEK + pERK');
rrme_reac_13_k = addkineticlaw(rrme_reac_13,'MassAction');
set(rrme_reac_13_k, 'ParameterVariableNames', {'koff_mek_erk'});

rrme_reac_14 = addreaction(pfpc_model, 'ppMEK_pERK -> ppMEK + ppERK');
rrme_reac_14_k = addkineticlaw(rrme_reac_14,'MassAction');
set(rrme_reac_14_k, 'ParameterVariableNames', {'kcat_mek'});

rrme_reac_15 = addreaction(pfpc_model, 'ppMEK + PP2A -> ppMEK_PP2A');
rrme_reac_15_k = addkineticlaw(rrme_reac_15,'MassAction');
set(rrme_reac_15_k, 'ParameterVariableNames', {'kon_mek_pp2a'});

rrme_reac_16 = addreaction(pfpc_model, 'ppMEK_PP2A -> ppMEK + PP2A');
rrme_reac_16_k = addkineticlaw(rrme_reac_16,'MassAction');
set(rrme_reac_16_k, 'ParameterVariableNames', {'koff_mek_pp2a'});

rrme_reac_17 = addreaction(pfpc_model, 'ppMEK_PP2A -> pMEK + PP2A');
rrme_reac_17_k = addkineticlaw(rrme_reac_17,'MassAction');
set(rrme_reac_17_k, 'ParameterVariableNames', {'kcat_pp2a'});

rrme_reac_18 = addreaction(pfpc_model, 'pMEK + PP2A -> pMEK_PP2A');
rrme_reac_18_k = addkineticlaw(rrme_reac_18,'MassAction');
set(rrme_reac_18_k, 'ParameterVariableNames', {'kon_mek_pp2a'});

rrme_reac_19 = addreaction(pfpc_model, 'pMEK_PP2A -> pMEK + PP2A');
rrme_reac_19_k = addkineticlaw(rrme_reac_19,'MassAction');
set(rrme_reac_19_k, 'ParameterVariableNames', {'koff_mek_pp2a'});

rrme_reac_20 = addreaction(pfpc_model, 'pMEK_PP2A -> MEK + PP2A');
rrme_reac_20_k = addkineticlaw(rrme_reac_20,'MassAction');
set(rrme_reac_20_k, 'ParameterVariableNames', {'kcat_pp2a'});

rrme_reac_21 = addreaction(pfpc_model, 'ppERK + MKP -> ppERK_MKP');
rrme_reac_21_k = addkineticlaw(rrme_reac_21,'MassAction');
set(rrme_reac_21_k, 'ParameterVariableNames', {'kon_erk_mkp'});

rrme_reac_22 = addreaction(pfpc_model, 'ppERK_MKP -> ppERK + MKP');
rrme_reac_22_k = addkineticlaw(rrme_reac_22,'MassAction');
set(rrme_reac_22_k, 'ParameterVariableNames', {'koff_erk_mkp'});

rrme_reac_23 = addreaction(pfpc_model, 'ppERK_MKP -> pERK + MKP');
rrme_reac_23_k = addkineticlaw(rrme_reac_23,'MassAction');
set(rrme_reac_23_k, 'ParameterVariableNames', {'kcat_mkp'});

rrme_reac_24 = addreaction(pfpc_model, 'pERK + MKP -> pERK_MKP');
rrme_reac_24_k = addkineticlaw(rrme_reac_24,'MassAction');
set(rrme_reac_24_k, 'ParameterVariableNames', {'kon_erk_mkp'});

rrme_reac_25 = addreaction(pfpc_model, 'pERK_MKP -> pERK + MKP');
rrme_reac_25_k = addkineticlaw(rrme_reac_25,'MassAction');
set(rrme_reac_25_k, 'ParameterVariableNames', {'koff_erk_mkp'});

rrme_reac_26 = addreaction(pfpc_model, 'pERK_MKP -> ERK + MKP');
rrme_reac_26_k = addkineticlaw(rrme_reac_26,'MassAction');
set(rrme_reac_26_k, 'ParameterVariableNames', {'kcat_mkp'});

rrme_reac_27 = addreaction(pfpc_model, 'ras_gtp -> ras_gdp');
rrme_reac_27_k = addkineticlaw(rrme_reac_27,'MassAction');
set(rrme_reac_27_k, 'ParameterVariableNames', {'kras_deact'});

rrme_reac_28 = addreaction(pfpc_model, 'pGEF -> GEF');
rrme_reac_28_k = addkineticlaw(rrme_reac_28,'MassAction');
set(rrme_reac_28_k, 'ParameterVariableNames', {'k_gef_deact'});

rrme_reac_29 = addreaction(pfpc_model, 'ras_gdp + pGEF -> ras_gdp_pGEF');
rrme_reac_29_k = addkineticlaw(rrme_reac_29,'MassAction');
set(rrme_reac_29_k, 'ParameterVariableNames', {'kon_ras_gef'});

rrme_reac_30 = addreaction(pfpc_model, 'ras_gdp_pGEF -> ras_gdp + pGEF');
rrme_reac_30_k = addkineticlaw(rrme_reac_30,'MassAction');
set(rrme_reac_30_k, 'ParameterVariableNames', {'koff_ras_gef'});

rrme_reac_31 = addreaction(pfpc_model, 'ras_gdp_pGEF -> ras_gtp + pGEF');
rrme_reac_31_k = addkineticlaw(rrme_reac_31,'MassAction');
set(rrme_reac_31_k, 'ParameterVariableNames', {'kcat_gef'});

rrme_reac_32 = addreaction(pfpc_model, 'ras_gtp + pGAP -> ras_gtp_pGAP');
rrme_reac_32_k = addkineticlaw(rrme_reac_32,'MassAction');
set(rrme_reac_32_k, 'ParameterVariableNames', {'kon_ras_gap'});

rrme_reac_33 = addreaction(pfpc_model, 'ras_gtp_pGAP -> ras_gtp + pGAP');
rrme_reac_33_k = addkineticlaw(rrme_reac_33,'MassAction');
set(rrme_reac_33_k, 'ParameterVariableNames', {'koff_ras_gap'});

rrme_reac_34 = addreaction(pfpc_model, 'ras_gtp_pGAP -> ras_gdp + pGAP');
rrme_reac_34_k = addkineticlaw(rrme_reac_34,'MassAction');
set(rrme_reac_34_k, 'ParameterVariableNames', {'kcat_gap'});

rrme_reac_35 = addreaction(pfpc_model, 'ras_RAF -> ras_RAF_act');
rrme_reac_35_k = addkineticlaw(rrme_reac_35,'MassAction');
set(rrme_reac_35_k, 'ParameterVariableNames', {'k_raf_act'});

rrme_reac_36 = addreaction(pfpc_model, 'ras_RAF_act -> ras_RAF');
rrme_reac_36_k = addkineticlaw(rrme_reac_36,'MassAction');
set(rrme_reac_36_k, 'ParameterVariableNames', {'k_raf_deact'});

rrme_reac_37 = addreaction(pfpc_model, 'RAF + RKIP -> RAF_RKIP');
rrme_reac_37_k = addkineticlaw(rrme_reac_37,'MassAction');
set(rrme_reac_37_k, 'ParameterVariableNames', {'kon_raf_rkip'});

rrme_reac_38 = addreaction(pfpc_model, 'RAF_RKIP -> RAF + RKIP');
rrme_reac_38_k = addkineticlaw(rrme_reac_38,'MassAction');
set(rrme_reac_38_k, 'ParameterVariableNames', {'koff_raf_rkip'});

rrme_reac_39 = addreaction(pfpc_model, 'pRKIP + RP -> pRKIP_RP');
rrme_reac_39_k = addkineticlaw(rrme_reac_39,'MassAction');
set(rrme_reac_39_k, 'ParameterVariableNames', {'kon_rkip_rp'});

rrme_reac_40 = addreaction(pfpc_model, 'pRKIP_RP -> pRKIP + RP');
rrme_reac_40_k = addkineticlaw(rrme_reac_40,'MassAction');
set(rrme_reac_40_k, 'ParameterVariableNames', {'koff_rkip_rp'});

rrme_reac_41 = addreaction(pfpc_model, 'pRKIP_RP -> RKIP + RP');
rrme_reac_41_k = addkineticlaw(rrme_reac_41,'MassAction');
set(rrme_reac_41_k, 'ParameterVariableNames', {'kcat_rp'});

rrme_reac_42 = addreaction(pfpc_model, 'ras_RAF_act + RAFK -> ras_RAF_RAFK');
rrme_reac_42_k = addkineticlaw(rrme_reac_42,'MassAction');
set(rrme_reac_42_k, 'ParameterVariableNames', {'kon_raf_rafk'});

rrme_reac_43 = addreaction(pfpc_model, 'ras_RAF_RAFK -> ras_RAF_act + RAFK');
rrme_reac_43_k = addkineticlaw(rrme_reac_43,'MassAction');
set(rrme_reac_43_k, 'ParameterVariableNames', {'koff_raf_rafk'});

rrme_reac_44 = addreaction(pfpc_model, 'ras_RAF_RAFK -> ras_pRAF_act + RAFK');
rrme_reac_44_k = addkineticlaw(rrme_reac_44,'MassAction');
set(rrme_reac_44_k, 'ParameterVariableNames', {'kcat_rafk'});

rrme_reac_45 = addreaction(pfpc_model, 'ras_pRAF_act -> ras_gtp + pRAF_act');
rrme_reac_45_k = addkineticlaw(rrme_reac_45,'MassAction');
set(rrme_reac_45_k, 'ParameterVariableNames', {'koff_ras_raf'});

rrme_reac_46 = addreaction(pfpc_model, 'ras_gtp + pRAF_act -> ras_pRAF_act');
rrme_reac_46_k = addkineticlaw(rrme_reac_46,'MassAction');
set(rrme_reac_46_k, 'ParameterVariableNames', {'kon_ras_raf'});

rrme_reac_47 = addreaction(pfpc_model, 'ras_pRAF_act + MEK -> ras_pRAF_MEK');
rrme_reac_47_k = addkineticlaw(rrme_reac_47,'MassAction');
set(rrme_reac_47_k, 'ParameterVariableNames', {'kon_raf_mek'});

rrme_reac_48 = addreaction(pfpc_model, 'ras_pRAF_MEK -> ras_pRAF_act + MEK');
rrme_reac_48_k = addkineticlaw(rrme_reac_48,'MassAction');
set(rrme_reac_48_k, 'ParameterVariableNames', {'koff_raf_mek'});

rrme_reac_49 = addreaction(pfpc_model, 'ras_pRAF_MEK -> ras_pRAF_act + pMEK');
rrme_reac_49_k = addkineticlaw(rrme_reac_49,'MassAction');
set(rrme_reac_49_k, 'ParameterVariableNames', {'kcat_raf'});

rrme_reac_50 = addreaction(pfpc_model, 'ras_pRAF_act + pMEK -> ras_pRAF_pMEK');
rrme_reac_50_k = addkineticlaw(rrme_reac_50,'MassAction');
set(rrme_reac_50_k, 'ParameterVariableNames', {'kon_raf_mek'});

rrme_reac_51 = addreaction(pfpc_model, 'ras_pRAF_pMEK -> ras_pRAF_act + pMEK');
rrme_reac_51_k = addkineticlaw(rrme_reac_51,'MassAction');
set(rrme_reac_51_k, 'ParameterVariableNames', {'koff_raf_mek'});

rrme_reac_52 = addreaction(pfpc_model, 'ras_pRAF_pMEK -> ras_pRAF_act + ppMEK');
rrme_reac_52_k = addkineticlaw(rrme_reac_52,'MassAction');
set(rrme_reac_52_k, 'ParameterVariableNames', {'kcat_raf'});

rrme_reac_53 = addreaction(pfpc_model, 'pRAF_act + MEK -> pRAF_MEK');
rrme_reac_53_k = addkineticlaw(rrme_reac_53,'MassAction');
set(rrme_reac_53_k, 'ParameterVariableNames', {'kon_raf_mek'});

rrme_reac_54 = addreaction(pfpc_model, 'pRAF_MEK -> pRAF_act + MEK');
rrme_reac_54_k = addkineticlaw(rrme_reac_54,'MassAction');
set(rrme_reac_54_k, 'ParameterVariableNames', {'koff_raf_mek'});

rrme_reac_55 = addreaction(pfpc_model, 'pRAF_MEK -> pRAF_act + pMEK');
rrme_reac_55_k = addkineticlaw(rrme_reac_55,'MassAction');
set(rrme_reac_55_k, 'ParameterVariableNames', {'kcat_raf'});

rrme_reac_56 = addreaction(pfpc_model, 'pRAF_act + pMEK -> pRAF_pMEK');
rrme_reac_56_k = addkineticlaw(rrme_reac_56,'MassAction');
set(rrme_reac_56_k, 'ParameterVariableNames', {'kon_raf_mek'});

rrme_reac_57 = addreaction(pfpc_model, 'pRAF_pMEK -> pRAF_act + pMEK');
rrme_reac_57_k = addkineticlaw(rrme_reac_57,'MassAction');
set(rrme_reac_57_k, 'ParameterVariableNames', {'koff_raf_mek'});

rrme_reac_58 = addreaction(pfpc_model, 'pRAF_pMEK -> pRAF_act + ppMEK');
rrme_reac_58_k = addkineticlaw(rrme_reac_58,'MassAction');
set(rrme_reac_58_k, 'ParameterVariableNames', {'kcat_raf'});

rrme_reac_59 = addreaction(pfpc_model, 'pRAF_act + PP5 -> pRAF_PP5');
rrme_reac_59_k = addkineticlaw(rrme_reac_59,'MassAction');
set(rrme_reac_59_k, 'ParameterVariableNames', {'kon_raf_pp5'});

rrme_reac_60 = addreaction(pfpc_model, 'pRAF_PP5 -> pRAF_act + PP5');
rrme_reac_60_k = addkineticlaw(rrme_reac_60,'MassAction');
set(rrme_reac_60_k, 'ParameterVariableNames', {'koff_raf_pp5'});

rrme_reac_61 = addreaction(pfpc_model, 'pRAF_PP5 -> RAF + PP5');
rrme_reac_61_k = addkineticlaw(rrme_reac_61,'MassAction');
set(rrme_reac_61_k, 'ParameterVariableNames', {'kcat_pp5'});

%%%%%%PLA2-AA Reactions%%%%%%%%%%

rrme_reac_62 = addreaction(pfpc_model, 'PLA2c + ppERK -> PLA2c_ppERK');
rrme_reac_62_k = addkineticlaw(rrme_reac_62,'MassAction');
set(rrme_reac_62_k, 'ParameterVariableNames', {'kon_pla2_erk'});

rrme_reac_63 = addreaction(pfpc_model, 'PLA2c_ppERK -> PLA2c + ppERK');
rrme_reac_63_k = addkineticlaw(rrme_reac_63,'MassAction');
set(rrme_reac_63_k, 'ParameterVariableNames', {'koff_pla2_erk'});

rrme_reac_64 = addreaction(pfpc_model, 'PLA2c_ppERK -> pPLA2c + ppERK');
rrme_reac_64_k = addkineticlaw(rrme_reac_64,'MassAction');
set(rrme_reac_64_k, 'ParameterVariableNames', {'kcat_erk'});

rrme_reac_65 = addreaction(pfpc_model, 'CaPLA2c + ppERK -> CaPLA2c_ppERK');
rrme_reac_65_k = addkineticlaw(rrme_reac_65,'MassAction');
set(rrme_reac_65_k, 'ParameterVariableNames', {'kon_pla2_erk'});

rrme_reac_66 = addreaction(pfpc_model, 'CaPLA2c_ppERK -> CaPLA2c + ppERK');
rrme_reac_66_k = addkineticlaw(rrme_reac_66,'MassAction');
set(rrme_reac_66_k, 'ParameterVariableNames', {'koff_pla2_erk'});

rrme_reac_67 = addreaction(pfpc_model, 'CaPLA2c_ppERK -> CapPLA2c + ppERK');
rrme_reac_67_k = addkineticlaw(rrme_reac_67,'MassAction');
set(rrme_reac_67_k, 'ParameterVariableNames', {'kcat_erk'});

rrme_reac_68 = addreaction(pfpc_model, 'Ca2PLA2c + ppERK -> Ca2PLA2c_ppERK');
rrme_reac_68_k = addkineticlaw(rrme_reac_68,'MassAction');
set(rrme_reac_68_k, 'ParameterVariableNames', {'kon_pla2_erk'});

rrme_reac_69 = addreaction(pfpc_model, 'Ca2PLA2c_ppERK -> Ca2PLA2c + ppERK');
rrme_reac_69_k = addkineticlaw(rrme_reac_69,'MassAction');
set(rrme_reac_69_k, 'ParameterVariableNames', {'koff_pla2_erk'});

rrme_reac_70 = addreaction(pfpc_model, 'Ca2PLA2c_ppERK -> Ca2pPLA2c + ppERK');
rrme_reac_70_k = addkineticlaw(rrme_reac_70,'MassAction');
set(rrme_reac_70_k, 'ParameterVariableNames', {'kcat_erk'});

rrme_reac_71 = addreaction(pfpc_model, 'pPLA2c + PP2A -> pPLA2c_PP2A');
rrme_reac_71_k = addkineticlaw(rrme_reac_71,'MassAction');
set(rrme_reac_71_k, 'ParameterVariableNames', {'kon_pla2_pp2a'});

rrme_reac_72 = addreaction(pfpc_model, 'pPLA2c_PP2A -> pPLA2c + PP2A');
rrme_reac_72_k = addkineticlaw(rrme_reac_72,'MassAction');
set(rrme_reac_72_k, 'ParameterVariableNames', {'koff_pla2_pp2a'});

rrme_reac_73 = addreaction(pfpc_model, 'pPLA2c_PP2A -> PLA2c + PP2A');
rrme_reac_73_k = addkineticlaw(rrme_reac_73,'MassAction');
set(rrme_reac_73_k, 'ParameterVariableNames', {'kcat_pp2a'});

rrme_reac_74 = addreaction(pfpc_model, 'CapPLA2c + PP2A -> CapPLA2c_PP2A');
rrme_reac_74_k = addkineticlaw(rrme_reac_74,'MassAction');
set(rrme_reac_74_k, 'ParameterVariableNames', {'kon_pla2_pp2a'});

rrme_reac_75 = addreaction(pfpc_model, 'CapPLA2c_PP2A -> CapPLA2c + PP2A');
rrme_reac_75_k = addkineticlaw(rrme_reac_75,'MassAction');
set(rrme_reac_75_k, 'ParameterVariableNames', {'koff_pla2_pp2a'});

rrme_reac_76 = addreaction(pfpc_model, 'CapPLA2c_PP2A -> CaPLA2c + PP2A');
rrme_reac_76_k = addkineticlaw(rrme_reac_76,'MassAction');
set(rrme_reac_76_k, 'ParameterVariableNames', {'kcat_pp2a'});

rrme_reac_77 = addreaction(pfpc_model, 'Ca2pPLA2c + PP2A -> Ca2pPLA2c_PP2A');
rrme_reac_77_k = addkineticlaw(rrme_reac_77,'MassAction');
set(rrme_reac_77_k, 'ParameterVariableNames', {'kon_pla2_pp2a'});

rrme_reac_78 = addreaction(pfpc_model, 'Ca2pPLA2c_PP2A -> Ca2pPLA2c + PP2A');
rrme_reac_78_k = addkineticlaw(rrme_reac_78,'MassAction');
set(rrme_reac_78_k, 'ParameterVariableNames', {'koff_pla2_pp2a'});

rrme_reac_79 = addreaction(pfpc_model, 'Ca2pPLA2c_PP2A -> Ca2PLA2c + PP2A');
rrme_reac_79_k = addkineticlaw(rrme_reac_79,'MassAction');
set(rrme_reac_79_k, 'ParameterVariableNames', {'kcat_pp2a'});

rrme_reac_80 = addreaction(pfpc_model, 'pPLA2c + Ca -> CapPLA2c');
rrme_reac_80_k = addkineticlaw(rrme_reac_80,'MassAction');
set(rrme_reac_80_k, 'ParameterVariableNames', {'kon_pla2_ca1'});

rrme_reac_81 = addreaction(pfpc_model, 'CapPLA2c -> pPLA2c + Ca');
rrme_reac_81_k = addkineticlaw(rrme_reac_81,'MassAction');
set(rrme_reac_81_k, 'ParameterVariableNames', {'koff_pla2_ca1'});

rrme_reac_82 = addreaction(pfpc_model, 'CapPLA2c + Ca -> Ca2pPLA2c');
rrme_reac_82_k = addkineticlaw(rrme_reac_82,'MassAction');
set(rrme_reac_82_k, 'ParameterVariableNames', {'kon_pla2_ca2'});

rrme_reac_83 = addreaction(pfpc_model, 'Ca2pPLA2c -> CapPLA2c + Ca');
rrme_reac_83_k = addkineticlaw(rrme_reac_83,'MassAction');
set(rrme_reac_83_k, 'ParameterVariableNames', {'koff_pla2_ca2'});

rrme_reac_84 = addreaction(pfpc_model, 'Ca2pPLA2c -> Ca2pPLA2');
rrme_reac_84_k = addkineticlaw(rrme_reac_84,'MassAction');
set(rrme_reac_84_k, 'ParameterVariableNames', {'kon_pla2_pm'});

rrme_reac_85 = addreaction(pfpc_model, 'Ca2pPLA2 -> Ca2pPLA2_act');
rrme_reac_85_k = addkineticlaw(rrme_reac_85,'MassAction');
set(rrme_reac_85_k, 'ParameterVariableNames', {'kact_pla2_pS505'});

rrme_reac_86 = addreaction(pfpc_model, 'Ca2pPLA2_act -> Ca2pPLA2');
rrme_reac_86_k = addkineticlaw(rrme_reac_86,'MassAction');
set(rrme_reac_86_k, 'ParameterVariableNames', {'kdeact_pla2_pS505'});

rrme_reac_87 = addreaction(pfpc_model, 'Ca2pPLA2 -> Ca2pPLA2c');
rrme_reac_87_k = addkineticlaw(rrme_reac_87,'MassAction');
set(rrme_reac_87_k, 'ParameterVariableNames', {'koff1_pla2_pm'});

rrme_reac_88 = addreaction(pfpc_model, 'Ca2pPLA2_act + DAG -> DAG_Ca2pPLA2_act');
rrme_reac_88_k = addkineticlaw(rrme_reac_88,'MassAction');
set(rrme_reac_88_k, 'ParameterVariableNames', {'kon_pla2_dag'});

rrme_reac_89 = addreaction(pfpc_model, 'DAG_Ca2pPLA2_act -> Ca2pPLA2_act + DAG');
rrme_reac_89_k = addkineticlaw(rrme_reac_89,'MassAction');
set(rrme_reac_89_k, 'ParameterVariableNames', {'koff_pla2_dag'});

rrme_reac_90 = addreaction(pfpc_model, 'PLA2c + Ca -> CaPLA2c');
rrme_reac_90_k = addkineticlaw(rrme_reac_90,'MassAction');
set(rrme_reac_90_k, 'ParameterVariableNames', {'kon_pla2_ca1'});

rrme_reac_91 = addreaction(pfpc_model, 'CaPLA2c -> PLA2c + Ca');
rrme_reac_91_k = addkineticlaw(rrme_reac_91,'MassAction');
set(rrme_reac_91_k, 'ParameterVariableNames', {'koff_pla2_ca1'});

rrme_reac_92 = addreaction(pfpc_model, 'CaPLA2c + Ca -> Ca2PLA2c');
rrme_reac_92_k = addkineticlaw(rrme_reac_92,'MassAction');
set(rrme_reac_92_k, 'ParameterVariableNames', {'kon_pla2_ca2'});

rrme_reac_93 = addreaction(pfpc_model, 'Ca2PLA2c -> CaPLA2c + Ca');
rrme_reac_93_k = addkineticlaw(rrme_reac_93,'MassAction');
set(rrme_reac_93_k, 'ParameterVariableNames', {'koff_pla2_ca2'});

rrme_reac_94 = addreaction(pfpc_model, 'Ca2PLA2c -> Ca2PLA2');
rrme_reac_94_k = addkineticlaw(rrme_reac_94,'MassAction');
set(rrme_reac_94_k, 'ParameterVariableNames', {'kon_pla2_pm'});

rrme_reac_95 = addreaction(pfpc_model, 'Ca2PLA2 -> Ca2PLA2_act');
rrme_reac_95_k = addkineticlaw(rrme_reac_95,'MassAction');
set(rrme_reac_95_k, 'ParameterVariableNames', {'kact_pla2'});

rrme_reac_96 = addreaction(pfpc_model, 'Ca2PLA2_act -> Ca2PLA2');
rrme_reac_96_k = addkineticlaw(rrme_reac_96,'MassAction');
set(rrme_reac_96_k, 'ParameterVariableNames', {'kdeact_pla2'});

rrme_reac_97 = addreaction(pfpc_model, 'Ca2PLA2 -> Ca2PLA2c');
rrme_reac_97_k = addkineticlaw(rrme_reac_97,'MassAction');
set(rrme_reac_97_k, 'ParameterVariableNames', {'koff1_pla2_pm'});

rrme_reac_98 = addreaction(pfpc_model, 'Ca2PLA2_act + DAG -> DAG_Ca2PLA2_act');
rrme_reac_98_k = addkineticlaw(rrme_reac_98,'MassAction');
set(rrme_reac_98_k, 'ParameterVariableNames', {'kon_pla2_dag'});

rrme_reac_99 = addreaction(pfpc_model, 'DAG_Ca2PLA2_act -> Ca2PLA2_act + DAG');
rrme_reac_99_k = addkineticlaw(rrme_reac_99,'MassAction');
set(rrme_reac_99_k, 'ParameterVariableNames', {'koff_pla2_dag'});

rrme_reac_100 = addreaction(pfpc_model, 'Ca2PLA2_act -> Ca2PLA2_act_APC');
rrme_reac_100_k = addkineticlaw(rrme_reac_100,'MassAction');
set(rrme_reac_100_k, 'ParameterVariableNames', {'kon_pla2_apc'});

rrme_reac_101 = addreaction(pfpc_model, 'Ca2PLA2_act_APC -> Ca2PLA2_act');
rrme_reac_101_k = addkineticlaw(rrme_reac_101,'MassAction');
set(rrme_reac_101_k, 'ParameterVariableNames', {'koff_pla2_apc'});

rrme_reac_102 = addreaction(pfpc_model, 'Ca2PLA2_act_APC -> Ca2PLA2_act + AA');
rrme_reac_102_k = addkineticlaw(rrme_reac_102,'MassAction');
set(rrme_reac_102_k, 'ParameterVariableNames', {'kcat_pla2'});

rrme_reac_103 = addreaction(pfpc_model, 'DAG_Ca2PLA2_act -> DAG_Ca2PLA2_act_APC');
rrme_reac_103_k = addkineticlaw(rrme_reac_103,'MassAction');
set(rrme_reac_103_k, 'ParameterVariableNames', {'kon_pla2_apc'});

rrme_reac_104 = addreaction(pfpc_model, 'DAG_Ca2PLA2_act_APC -> DAG_Ca2PLA2_act');
rrme_reac_104_k = addkineticlaw(rrme_reac_104,'MassAction');
set(rrme_reac_104_k, 'ParameterVariableNames', {'koff_pla2_apc'});

rrme_reac_105 = addreaction(pfpc_model, 'DAG_Ca2PLA2_act_APC -> DAG_Ca2PLA2_act + AA');
rrme_reac_105_k = addkineticlaw(rrme_reac_105,'MassAction');
set(rrme_reac_105_k, 'ParameterVariableNames', {'kcat_pla2'});

rrme_reac_106 = addreaction(pfpc_model, 'Ca2pPLA2_act -> Ca2pPLA2_act_APC');
rrme_reac_106_k = addkineticlaw(rrme_reac_106,'MassAction');
set(rrme_reac_106_k, 'ParameterVariableNames', {'kon_pla2_apc'});

rrme_reac_107 = addreaction(pfpc_model, 'Ca2pPLA2_act_APC -> Ca2pPLA2_act');
rrme_reac_107_k = addkineticlaw(rrme_reac_107,'MassAction');
set(rrme_reac_107_k, 'ParameterVariableNames', {'koff_pla2_apc'});

rrme_reac_108 = addreaction(pfpc_model, 'Ca2pPLA2_act_APC -> Ca2pPLA2_act + AA');
rrme_reac_108_k = addkineticlaw(rrme_reac_108,'MassAction');
set(rrme_reac_108_k, 'ParameterVariableNames', {'kcat_pla2'});

rrme_reac_109 = addreaction(pfpc_model, 'DAG_Ca2pPLA2_act -> DAG_Ca2pPLA2_act_APC');
rrme_reac_109_k = addkineticlaw(rrme_reac_109,'MassAction');
set(rrme_reac_109_k, 'ParameterVariableNames', {'kon_pla2_apc'});

rrme_reac_110 = addreaction(pfpc_model, 'DAG_Ca2pPLA2_act_APC -> DAG_Ca2pPLA2_act');
rrme_reac_110_k = addkineticlaw(rrme_reac_110,'MassAction');
set(rrme_reac_110_k, 'ParameterVariableNames', {'koff_pla2_apc'});

rrme_reac_111 = addreaction(pfpc_model, 'DAG_Ca2pPLA2_act_APC -> DAG_Ca2pPLA2_act + AA');
rrme_reac_111_k = addkineticlaw(rrme_reac_111,'MassAction');
set(rrme_reac_111_k, 'ParameterVariableNames', {'kcat_pla2'});

rrme_reac_112 = addreaction(pfpc_model, 'AA -> null');
rrme_reac_112_k = addkineticlaw(rrme_reac_112,'MassAction');
set(rrme_reac_112_k, 'ParameterVariableNames', {'kdeg_aa'});

%%%%%%%PKC REACTIONS%%%%%%%%%
rrme_reac_144 = addreaction(pfpc_model, 'RAF_RKIP + CaPKC_act -> RAF_RKIP_CaPKC_act');
rrme_reac_144_k = addkineticlaw(rrme_reac_144,'MassAction');
set(rrme_reac_144_k, 'ParameterVariableNames', {'kon_rkip_pkc'});

rrme_reac_145 = addreaction(pfpc_model, 'RAF_RKIP_CaPKC_act -> RAF_RKIP + CaPKC_act');
rrme_reac_145_k = addkineticlaw(rrme_reac_145,'MassAction');
set(rrme_reac_145_k, 'ParameterVariableNames', {'koff_rkip_pkc'});

rrme_reac_146 = addreaction(pfpc_model, 'RAF_RKIP_CaPKC_act -> RAF + pRKIP + CaPKC_act');
rrme_reac_146_k = addkineticlaw(rrme_reac_146,'MassAction');
set(rrme_reac_146_k, 'ParameterVariableNames', {'kcat_pkc'});

rrme_reac_147 = addreaction(pfpc_model, 'RAF_RKIP + Ca2PKC_act -> RAF_RKIP_Ca2PKC_act');
rrme_reac_147_k = addkineticlaw(rrme_reac_147,'MassAction');
set(rrme_reac_147_k, 'ParameterVariableNames', {'kon_rkip_pkc'});

rrme_reac_148 = addreaction(pfpc_model, 'RAF_RKIP_Ca2PKC_act -> RAF_RKIP + Ca2PKC_act');
rrme_reac_148_k = addkineticlaw(rrme_reac_148, 'MassAction');
set(rrme_reac_148_k, 'ParameterVariableNames', {'koff_rkip_pkc'});

rrme_reac_149 = addreaction(pfpc_model, 'RAF_RKIP_Ca2PKC_act -> RAF + pRKIP + Ca2PKC_act');
rrme_reac_149_k = addkineticlaw(rrme_reac_149,'MassAction');
set(rrme_reac_149_k, 'ParameterVariableNames', {'kcat_pkc'});

rrme_reac_150 = addreaction(pfpc_model, 'RAF_RKIP + CaPKC_AA_act -> RAF_RKIP_CaPKC_AA_act');
rrme_reac_150_k = addkineticlaw(rrme_reac_150,'MassAction');
set(rrme_reac_150_k, 'ParameterVariableNames', {'kon_rkip_pkc'});

rrme_reac_151 = addreaction(pfpc_model, 'RAF_RKIP_CaPKC_AA_act -> RAF_RKIP + CaPKC_AA_act');
rrme_reac_151_k = addkineticlaw(rrme_reac_151,'MassAction');
set(rrme_reac_151_k, 'ParameterVariableNames', {'koff_rkip_pkc'});

rrme_reac_152 = addreaction(pfpc_model, 'RAF_RKIP_CaPKC_AA_act -> RAF + pRKIP + CaPKC_AA_act');
rrme_reac_152_k = addkineticlaw(rrme_reac_152,'MassAction');
set(rrme_reac_152_k, 'ParameterVariableNames', {'kcat_pkc'});

rrme_reac_153 = addreaction(pfpc_model, 'RAF_RKIP + Ca2PKC_AA_act -> RAF_RKIP_Ca2PKC_AA_act');
rrme_reac_153_k = addkineticlaw(rrme_reac_153,'MassAction');
set(rrme_reac_153_k, 'ParameterVariableNames', {'kon_rkip_pkc'});

rrme_reac_154 = addreaction(pfpc_model, 'RAF_RKIP_Ca2PKC_AA_act -> RAF_RKIP + Ca2PKC_AA_act');
rrme_reac_154_k = addkineticlaw(rrme_reac_154,'MassAction');
set(rrme_reac_154_k, 'ParameterVariableNames', {'koff_rkip_pkc'});

rrme_reac_155 = addreaction(pfpc_model, 'RAF_RKIP_Ca2PKC_AA_act -> RAF + pRKIP + Ca2PKC_AA_act');
rrme_reac_155_k = addkineticlaw(rrme_reac_155,'MassAction');
set(rrme_reac_155_k, 'ParameterVariableNames', {'kcat_pkc'});

rrme_reac_156 = addreaction(pfpc_model, 'RAF_RKIP + DAG_CaPKC_act -> RAF_RKIP_DAG_CaPKC_act');
rrme_reac_156_k = addkineticlaw(rrme_reac_156,'MassAction');
set(rrme_reac_156_k, 'ParameterVariableNames', {'kon_rkip_pkc'});

rrme_reac_157 = addreaction(pfpc_model, 'RAF_RKIP_DAG_CaPKC_act -> RAF_RKIP + DAG_CaPKC_act');
rrme_reac_157_k = addkineticlaw(rrme_reac_157,'MassAction');
set(rrme_reac_157_k, 'ParameterVariableNames', {'koff_rkip_pkc'});

rrme_reac_158 = addreaction(pfpc_model, 'RAF_RKIP_DAG_CaPKC_act -> RAF + pRKIP + DAG_CaPKC_act');
rrme_reac_158_k = addkineticlaw(rrme_reac_158,'MassAction');
set(rrme_reac_158_k, 'ParameterVariableNames', {'kcat_pkc'});

rrme_reac_159 = addreaction(pfpc_model, 'RAF_RKIP + DAG_Ca2PKC_act -> RAF_RKIP_DAG_Ca2PKC_act');
rrme_reac_159_k = addkineticlaw(rrme_reac_159,'MassAction');
set(rrme_reac_159_k, 'ParameterVariableNames', {'kon_rkip_pkc'});

rrme_reac_160 = addreaction(pfpc_model, 'RAF_RKIP_DAG_Ca2PKC_act -> RAF_RKIP + DAG_Ca2PKC_act');
rrme_reac_160_k = addkineticlaw(rrme_reac_160,'MassAction');
set(rrme_reac_160_k, 'ParameterVariableNames', {'koff_rkip_pkc'});

rrme_reac_161 = addreaction(pfpc_model, 'RAF_RKIP_DAG_Ca2PKC_act -> RAF + pRKIP + DAG_Ca2PKC_act');
rrme_reac_161_k = addkineticlaw(rrme_reac_161,'MassAction');
set(rrme_reac_161_k, 'ParameterVariableNames', {'kcat_pkc'});

rrme_reac_162 = addreaction(pfpc_model, 'RAF_RKIP + DAG_CaPKC_AA_act -> RAF_RKIP_DAG_CaPKC_AA_act');
rrme_reac_162_k = addkineticlaw(rrme_reac_162,'MassAction');
set(rrme_reac_162_k, 'ParameterVariableNames', {'kon_rkip_pkc'});

rrme_reac_163 = addreaction(pfpc_model, 'RAF_RKIP_DAG_CaPKC_AA_act -> RAF_RKIP + DAG_CaPKC_AA_act');
rrme_reac_163_k = addkineticlaw(rrme_reac_163,'MassAction');
set(rrme_reac_163_k, 'ParameterVariableNames', {'koff_rkip_pkc'});

rrme_reac_164 = addreaction(pfpc_model, 'RAF_RKIP_DAG_CaPKC_AA_act -> RAF + pRKIP + DAG_CaPKC_AA_act');
rrme_reac_164_k = addkineticlaw(rrme_reac_164,'MassAction');
set(rrme_reac_164_k, 'ParameterVariableNames', {'kcat_pkc'});

rrme_reac_165 = addreaction(pfpc_model, 'RAF_RKIP + DAG_Ca2PKC_AA_act -> RAF_RKIP_DAG_Ca2PKC_AA_act');
rrme_reac_165_k = addkineticlaw(rrme_reac_165,'MassAction');
set(rrme_reac_165_k, 'ParameterVariableNames', {'kon_rkip_pkc'});

rrme_reac_166 = addreaction(pfpc_model, 'RAF_RKIP_DAG_Ca2PKC_AA_act -> RAF_RKIP + DAG_Ca2PKC_AA_act');
rrme_reac_166_k = addkineticlaw(rrme_reac_166,'MassAction');
set(rrme_reac_166_k, 'ParameterVariableNames', {'koff_rkip_pkc'});

rrme_reac_167 = addreaction(pfpc_model, 'RAF_RKIP_DAG_Ca2PKC_AA_act -> RAF + pRKIP + DAG_Ca2PKC_AA_act');
rrme_reac_167_k = addkineticlaw(rrme_reac_167,'MassAction');
set(rrme_reac_167_k, 'ParameterVariableNames', {'kcat_pkc'});

rrme_reac_168 = addreaction(pfpc_model, 'GEF + CaPKC_act -> GEF_CaPKC_act');
rrme_reac_168_k = addkineticlaw(rrme_reac_168,'MassAction');
set(rrme_reac_168_k, 'ParameterVariableNames', {'kon_pkc_gef'});

rrme_reac_169 = addreaction(pfpc_model, 'GEF_CaPKC_act -> GEF + CaPKC_act');
rrme_reac_169_k = addkineticlaw(rrme_reac_169,'MassAction');
set(rrme_reac_169_k, 'ParameterVariableNames', {'koff_pkc_gef'});

rrme_reac_170 = addreaction(pfpc_model, 'GEF_CaPKC_act -> pGEF + CaPKC_act');
rrme_reac_170_k = addkineticlaw(rrme_reac_170,'MassAction');
set(rrme_reac_170_k, 'ParameterVariableNames', {'kcat_pkc_gef'});

rrme_reac_171 = addreaction(pfpc_model, 'GEF + Ca2PKC_act -> GEF_Ca2PKC_act');
rrme_reac_171_k = addkineticlaw(rrme_reac_171,'MassAction');
set(rrme_reac_171_k, 'ParameterVariableNames', {'kon_pkc_gef'});

rrme_reac_172 = addreaction(pfpc_model, 'GEF_Ca2PKC_act -> GEF + Ca2PKC_act');
rrme_reac_172_k = addkineticlaw(rrme_reac_172,'MassAction');
set(rrme_reac_172_k, 'ParameterVariableNames', {'koff_pkc_gef'});

rrme_reac_173 = addreaction(pfpc_model, 'GEF_Ca2PKC_act -> pGEF + Ca2PKC_act');
rrme_reac_173_k = addkineticlaw(rrme_reac_173,'MassAction');
set(rrme_reac_173_k, 'ParameterVariableNames', {'kcat_pkc_gef'});

rrme_reac_174 = addreaction(pfpc_model, 'GEF + CaPKC_AA_act -> GEF_CaPKC_AA_act');
rrme_reac_174_k = addkineticlaw(rrme_reac_174,'MassAction');
set(rrme_reac_174_k, 'ParameterVariableNames', {'kon_pkc_gef'});

rrme_reac_175 = addreaction(pfpc_model, 'GEF_CaPKC_AA_act -> GEF + CaPKC_AA_act');
rrme_reac_175_k = addkineticlaw(rrme_reac_175,'MassAction');
set(rrme_reac_175_k, 'ParameterVariableNames', {'koff_pkc_gef'});

rrme_reac_176 = addreaction(pfpc_model, 'GEF_CaPKC_AA_act -> pGEF + CaPKC_AA_act');
rrme_reac_176_k = addkineticlaw(rrme_reac_176,'MassAction');
set(rrme_reac_176_k, 'ParameterVariableNames', {'kcat_pkc_gef'});

rrme_reac_177 = addreaction(pfpc_model, 'GEF + Ca2PKC_AA_act -> GEF_Ca2PKC_AA_act');
rrme_reac_177_k = addkineticlaw(rrme_reac_177,'MassAction');
set(rrme_reac_177_k, 'ParameterVariableNames', {'kon_pkc_gef'});

rrme_reac_178 = addreaction(pfpc_model, 'GEF_Ca2PKC_AA_act -> GEF + Ca2PKC_AA_act');
rrme_reac_178_k = addkineticlaw(rrme_reac_178,'MassAction');
set(rrme_reac_178_k, 'ParameterVariableNames', {'koff_pkc_gef'});

rrme_reac_179 = addreaction(pfpc_model, 'GEF_Ca2PKC_AA_act -> pGEF + Ca2PKC_AA_act');
rrme_reac_179_k = addkineticlaw(rrme_reac_179,'MassAction');
set(rrme_reac_179_k, 'ParameterVariableNames', {'kcat_pkc_gef'});

rrme_reac_180 = addreaction(pfpc_model, 'GEF + DAG_CaPKC_act -> GEF_DAG_CaPKC_act');
rrme_reac_180_k = addkineticlaw(rrme_reac_180,'MassAction');
set(rrme_reac_180_k, 'ParameterVariableNames', {'kon_pkc_gef'});

rrme_reac_181 = addreaction(pfpc_model, 'GEF_DAG_CaPKC_act -> GEF + DAG_CaPKC_act');
rrme_reac_181_k = addkineticlaw(rrme_reac_181,'MassAction');
set(rrme_reac_181_k, 'ParameterVariableNames', {'koff_pkc_gef'});

rrme_reac_182 = addreaction(pfpc_model, 'GEF_DAG_CaPKC_act -> pGEF + DAG_CaPKC_act');
rrme_reac_182_k = addkineticlaw(rrme_reac_182,'MassAction');
set(rrme_reac_182_k, 'ParameterVariableNames', {'kcat_pkc_gef'});

rrme_reac_183 = addreaction(pfpc_model, 'GEF + DAG_Ca2PKC_act -> GEF_DAG_Ca2PKC_act');
rrme_reac_183_k = addkineticlaw(rrme_reac_183,'MassAction');
set(rrme_reac_183_k, 'ParameterVariableNames', {'kon_pkc_gef'});

rrme_reac_184 = addreaction(pfpc_model, 'GEF_DAG_Ca2PKC_act -> GEF + DAG_Ca2PKC_act');
rrme_reac_184_k = addkineticlaw(rrme_reac_184,'MassAction');
set(rrme_reac_184_k, 'ParameterVariableNames', {'koff_pkc_gef'});

rrme_reac_185 = addreaction(pfpc_model, 'GEF_DAG_Ca2PKC_act -> pGEF + DAG_Ca2PKC_act');
rrme_reac_185_k = addkineticlaw(rrme_reac_185,'MassAction');
set(rrme_reac_185_k, 'ParameterVariableNames', {'kcat_pkc_gef'});

rrme_reac_186 = addreaction(pfpc_model, 'GEF + DAG_CaPKC_AA_act -> GEF_DAG_CaPKC_AA_act');
rrme_reac_186_k = addkineticlaw(rrme_reac_186,'MassAction');
set(rrme_reac_186_k, 'ParameterVariableNames', {'kon_pkc_gef'});

rrme_reac_187 = addreaction(pfpc_model, 'GEF_DAG_CaPKC_AA_act -> GEF + DAG_CaPKC_AA_act');
rrme_reac_187_k = addkineticlaw(rrme_reac_187,'MassAction');
set(rrme_reac_187_k, 'ParameterVariableNames', {'koff_pkc_gef'});

rrme_reac_188 = addreaction(pfpc_model, 'GEF_DAG_CaPKC_AA_act -> pGEF + DAG_CaPKC_AA_act');
rrme_reac_188_k = addkineticlaw(rrme_reac_188,'MassAction');
set(rrme_reac_188_k, 'ParameterVariableNames', {'kcat_pkc_gef'});

rrme_reac_189 = addreaction(pfpc_model, 'GEF + DAG_Ca2PKC_AA_act -> GEF_DAG_Ca2PKC_AA_act');
rrme_reac_189_k = addkineticlaw(rrme_reac_189,'MassAction');
set(rrme_reac_189_k, 'ParameterVariableNames', {'kon_pkc_gef'});

rrme_reac_190 = addreaction(pfpc_model, 'GEF_DAG_Ca2PKC_AA_act -> GEF + DAG_Ca2PKC_AA_act');
rrme_reac_190_k = addkineticlaw(rrme_reac_190,'MassAction');
set(rrme_reac_190_k, 'ParameterVariableNames', {'koff_pkc_gef'});

rrme_reac_191 = addreaction(pfpc_model, 'GEF_DAG_Ca2PKC_AA_act -> pGEF + DAG_Ca2PKC_AA_act');
rrme_reac_191_k = addkineticlaw(rrme_reac_191,'MassAction');
set(rrme_reac_191_k, 'ParameterVariableNames', {'kcat_pkc_gef'});

rrme_reac_192 = addreaction(pfpc_model, 'pPLA2c + PP1 -> pPLA2c_PP1');
rrme_reac_192_k = addkineticlaw(rrme_reac_192,'MassAction');
set(rrme_reac_192_k, 'ParameterVariableNames', {'kon_pla2_pp1'});

rrme_reac_193 = addreaction(pfpc_model, 'pPLA2c_PP1 -> pPLA2c + PP1');
rrme_reac_193_k = addkineticlaw(rrme_reac_193,'MassAction');
set(rrme_reac_193_k, 'ParameterVariableNames', {'koff_pla2_pp1'});

rrme_reac_194 = addreaction(pfpc_model, 'pPLA2c_PP1 -> PLA2c + PP1');
rrme_reac_194_k = addkineticlaw(rrme_reac_194,'MassAction');
set(rrme_reac_194_k, 'ParameterVariableNames', {'kcat_pp1'});

rrme_reac_195 = addreaction(pfpc_model, 'CapPLA2c + PP1 -> CapPLA2c_PP1');
rrme_reac_195_k = addkineticlaw(rrme_reac_195,'MassAction');
set(rrme_reac_195_k, 'ParameterVariableNames', {'kon_pla2_pp1'});

rrme_reac_196 = addreaction(pfpc_model, 'CapPLA2c_PP1 -> CapPLA2c + PP1');
rrme_reac_196_k = addkineticlaw(rrme_reac_196,'MassAction');
set(rrme_reac_196_k, 'ParameterVariableNames', {'koff_pla2_pp1'});

rrme_reac_197 = addreaction(pfpc_model, 'CapPLA2c_PP1 -> CaPLA2c + PP1');
rrme_reac_197_k = addkineticlaw(rrme_reac_197,'MassAction');
set(rrme_reac_197_k, 'ParameterVariableNames', {'kcat_pp1'});

rrme_reac_198 = addreaction(pfpc_model, 'Ca2pPLA2c + PP1 -> Ca2pPLA2c_PP1');
rrme_reac_198_k = addkineticlaw(rrme_reac_198,'MassAction');
set(rrme_reac_198_k, 'ParameterVariableNames', {'kon_pla2_pp1'});

rrme_reac_199 = addreaction(pfpc_model, 'Ca2pPLA2c_PP1 -> Ca2pPLA2c + PP1');
rrme_reac_199_k = addkineticlaw(rrme_reac_199,'MassAction');
set(rrme_reac_199_k, 'ParameterVariableNames', {'koff_pla2_pp1'});

rrme_reac_200 = addreaction(pfpc_model, 'Ca2pPLA2c_PP1 -> Ca2PLA2c + PP1');
rrme_reac_200_k = addkineticlaw(rrme_reac_200,'MassAction');
set(rrme_reac_200_k, 'ParameterVariableNames', {'kcat_pp1'});

rrme_reac_201 = addreaction(pfpc_model, 'PKCc + Ca -> CaPKCc');
rrme_reac_201_k = addkineticlaw(rrme_reac_201,'MassAction');
set(rrme_reac_201_k, 'ParameterVariableNames', {'kon_pkc_ca1'});

rrme_reac_202 = addreaction(pfpc_model, 'CaPKCc -> PKCc + Ca');
rrme_reac_202_k = addkineticlaw(rrme_reac_202,'MassAction');
set(rrme_reac_202_k, 'ParameterVariableNames', {'koff_pkc_ca1'});

rrme_reac_203 = addreaction(pfpc_model, 'CaPKCc -> CaPKC');
rrme_reac_203_k = addkineticlaw(rrme_reac_203,'MassAction');
set(rrme_reac_203_k, 'ParameterVariableNames', {'kon_pkc_pm'});

rrme_reac_204 = addreaction(pfpc_model, 'CaPKC -> CaPKCc');
rrme_reac_204_k = addkineticlaw(rrme_reac_204,'MassAction');
set(rrme_reac_204_k, 'ParameterVariableNames', {'koff_pkc_pm'});

rrme_reac_205 = addreaction(pfpc_model, 'CaPKC + Ca -> Ca2PKC');
rrme_reac_205_k = addkineticlaw(rrme_reac_205,'MassAction');
set(rrme_reac_205_k, 'ParameterVariableNames', {'kon_pkc_ca2'});

rrme_reac_206 = addreaction(pfpc_model, 'Ca2PKC -> CaPKC + Ca');
rrme_reac_206_k = addkineticlaw(rrme_reac_206,'MassAction');
set(rrme_reac_206_k, 'ParameterVariableNames', {'koff_pkc_ca2'});

rrme_reac_207 = addreaction(pfpc_model, 'Ca2PKC -> Ca2PKC_act');
rrme_reac_207_k = addkineticlaw(rrme_reac_207,'MassAction');
set(rrme_reac_207_k, 'ParameterVariableNames', {'kact_pkc'});

rrme_reac_208 = addreaction(pfpc_model, 'Ca2PKC_act -> Ca2PKC');
rrme_reac_208_k = addkineticlaw(rrme_reac_208,'MassAction');
set(rrme_reac_208_k, 'ParameterVariableNames', {'kdeact_pkc'});

rrme_reac_209 = addreaction(pfpc_model, 'Ca2PKC_act + DAG -> DAG_Ca2PKC_act');
rrme_reac_209_k = addkineticlaw(rrme_reac_209,'MassAction');
set(rrme_reac_209_k, 'ParameterVariableNames', {'kon_pkc_dag'});

rrme_reac_210 = addreaction(pfpc_model, 'DAG_Ca2PKC_act -> Ca2PKC_act + DAG');
rrme_reac_210_k = addkineticlaw(rrme_reac_210,'MassAction');
set(rrme_reac_210_k, 'ParameterVariableNames', {'koff_pkc_dag'});

rrme_reac_211 = addreaction(pfpc_model, 'DAG_CaPKC_act + Ca -> DAG_Ca2PKC_act');
rrme_reac_211_k = addkineticlaw(rrme_reac_211,'MassAction');
set(rrme_reac_211_k, 'ParameterVariableNames', {'kon_pkc_ca2'});

rrme_reac_212 = addreaction(pfpc_model, 'DAG_Ca2PKC_act -> DAG_CaPKC_act + Ca');
rrme_reac_212_k = addkineticlaw(rrme_reac_212,'MassAction');
set(rrme_reac_212_k, 'ParameterVariableNames', {'koff_pkc_ca2'});

rrme_reac_213 = addreaction(pfpc_model, 'CaPKC_act + DAG -> DAG_CaPKC_act');
rrme_reac_213_k = addkineticlaw(rrme_reac_213,'MassAction');
set(rrme_reac_213_k, 'ParameterVariableNames', {'kon_pkc_dag'});

rrme_reac_214 = addreaction(pfpc_model, 'DAG_CaPKC_act -> CaPKC_act + DAG');
rrme_reac_214_k = addkineticlaw(rrme_reac_214,'MassAction');
set(rrme_reac_214_k, 'ParameterVariableNames', {'koff_pkc_dag'});

rrme_reac_215 = addreaction(pfpc_model, 'CaPKC_act + Ca -> Ca2PKC_act');
rrme_reac_215_k = addkineticlaw(rrme_reac_215,'MassAction');
set(rrme_reac_215_k, 'ParameterVariableNames', {'kon_pkc_ca2'});

rrme_reac_216 = addreaction(pfpc_model, 'Ca2PKC_act -> CaPKC_act + Ca');
rrme_reac_216_k = addkineticlaw(rrme_reac_216,'MassAction');
set(rrme_reac_216_k, 'ParameterVariableNames', {'koff_pkc_ca2'});

rrme_reac_217 = addreaction(pfpc_model, 'CaPKC_act -> CaPKC');
rrme_reac_217_k = addkineticlaw(rrme_reac_217,'MassAction');
set(rrme_reac_217_k, 'ParameterVariableNames', {'kdeact_pkc'});

rrme_reac_218 = addreaction(pfpc_model, 'CaPKC + AA -> CaPKC_AA');
rrme_reac_218_k = addkineticlaw(rrme_reac_218,'MassAction');
set(rrme_reac_218_k, 'ParameterVariableNames', {'kon_pkc_aa'});

rrme_reac_219 = addreaction(pfpc_model, 'CaPKC_AA -> CaPKC + AA');
rrme_reac_219_k = addkineticlaw(rrme_reac_219,'MassAction');
set(rrme_reac_219_k, 'ParameterVariableNames', {'koff_pkc_aa'});

rrme_reac_220 = addreaction(pfpc_model, 'Ca2PKC + AA -> Ca2PKC_AA');
rrme_reac_220_k = addkineticlaw(rrme_reac_220,'MassAction');
set(rrme_reac_220_k, 'ParameterVariableNames', {'kon_pkc_aa'});

rrme_reac_221 = addreaction(pfpc_model, 'Ca2PKC_AA -> Ca2PKC + AA');
rrme_reac_221_k = addkineticlaw(rrme_reac_221,'MassAction');
set(rrme_reac_221_k, 'ParameterVariableNames', {'koff_pkc_aa'});

rrme_reac_222 = addreaction(pfpc_model, 'CaPKC_AA + Ca -> Ca2PKC_AA');
rrme_reac_222_k = addkineticlaw(rrme_reac_222,'MassAction');
set(rrme_reac_222_k, 'ParameterVariableNames', {'kon_pkc_ca2_aa'});

rrme_reac_223 = addreaction(pfpc_model, 'Ca2PKC_AA -> CaPKC_AA + Ca');
rrme_reac_223_k = addkineticlaw(rrme_reac_223,'MassAction');
set(rrme_reac_223_k, 'ParameterVariableNames', {'koff_pkc_ca2'});

rrme_reac_224 = addreaction(pfpc_model, 'Ca2PKC_AA -> Ca2PKC_AA_act');
rrme_reac_224_k = addkineticlaw(rrme_reac_224,'MassAction');
set(rrme_reac_224_k, 'ParameterVariableNames', {'kact_pkc'});

rrme_reac_225 = addreaction(pfpc_model, 'Ca2PKC_AA_act -> Ca2PKC_AA');
rrme_reac_225_k = addkineticlaw(rrme_reac_225,'MassAction');
set(rrme_reac_225_k, 'ParameterVariableNames', {'kdeact_pkc_aa'});

rrme_reac_226 = addreaction(pfpc_model, 'Ca2PKC_AA_act + DAG -> DAG_Ca2PKC_AA_act');
rrme_reac_226_k = addkineticlaw(rrme_reac_226,'MassAction');
set(rrme_reac_226_k, 'ParameterVariableNames', {'kon_pkc_dag'});

rrme_reac_227 = addreaction(pfpc_model, 'DAG_Ca2PKC_AA_act -> Ca2PKC_AA_act + DAG');
rrme_reac_227_k = addkineticlaw(rrme_reac_227,'MassAction');
set(rrme_reac_227_k, 'ParameterVariableNames', {'koff_pkc_dag'});

rrme_reac_228 = addreaction(pfpc_model, 'DAG_CaPKC_AA_act + Ca -> DAG_Ca2PKC_AA_act');
rrme_reac_228_k = addkineticlaw(rrme_reac_228,'MassAction');
set(rrme_reac_228_k, 'ParameterVariableNames', {'kon_pkc_ca2_aa'});

rrme_reac_229 = addreaction(pfpc_model, 'DAG_Ca2PKC_AA_act -> DAG_CaPKC_AA_act + Ca');
rrme_reac_229_k = addkineticlaw(rrme_reac_229,'MassAction');
set(rrme_reac_229_k, 'ParameterVariableNames', {'koff_pkc_ca2'});

rrme_reac_230 = addreaction(pfpc_model, 'CaPKC_AA_act + DAG -> DAG_CaPKC_AA_act');
rrme_reac_230_k = addkineticlaw(rrme_reac_230,'MassAction');
set(rrme_reac_230_k, 'ParameterVariableNames', {'kon_pkc_dag'});

rrme_reac_231 = addreaction(pfpc_model, 'DAG_CaPKC_AA_act -> CaPKC_AA_act + DAG');
rrme_reac_231_k = addkineticlaw(rrme_reac_231,'MassAction');
set(rrme_reac_231_k, 'ParameterVariableNames', {'koff_pkc_dag'});

rrme_reac_232 = addreaction(pfpc_model, 'CaPKC_AA_act + DAG -> DAG_CaPKC_AA_act');
rrme_reac_232_k = addkineticlaw(rrme_reac_232,'MassAction');
set(rrme_reac_232_k, 'ParameterVariableNames', {'kon_pkc_dag'});

rrme_reac_233 = addreaction(pfpc_model, 'DAG_CaPKC_AA_act -> CaPKC_AA_act + DAG');
rrme_reac_233_k = addkineticlaw(rrme_reac_233,'MassAction');
set(rrme_reac_233_k, 'ParameterVariableNames', {'koff_pkc_dag'});

rrme_reac_234 = addreaction(pfpc_model, 'CaPKC_AA_act + Ca -> Ca2PKC_AA_act');
rrme_reac_234_k = addkineticlaw(rrme_reac_234,'MassAction');
set(rrme_reac_234_k, 'ParameterVariableNames', {'kon_pkc_ca2_aa'});

rrme_reac_235 = addreaction(pfpc_model, 'Ca2PKC_AA_act -> CaPKC_AA_act + Ca');
rrme_reac_235_k = addkineticlaw(rrme_reac_235,'MassAction');
set(rrme_reac_235_k, 'ParameterVariableNames', {'koff_pkc_ca2'});

rrme_reac_236 = addreaction(pfpc_model, 'CaPKC_AA_act -> CaPKC_AA');
rrme_reac_236_k = addkineticlaw(rrme_reac_236,'MassAction');
set(rrme_reac_236_k, 'ParameterVariableNames', {'kdeact_pkc_aa'});

rrme_reac_237 = addreaction(pfpc_model, 'CaPKC_AA -> CaPKC_AA_act');
rrme_reac_237_k = addkineticlaw(rrme_reac_237,'MassAction');
set(rrme_reac_237_k, 'ParameterVariableNames', {'kact_pkc'});

rrme_reac_238 = addreaction(pfpc_model, 'DAG_Ca2PKC_act + AA -> DAG_Ca2PKC_AA_act');
rrme_reac_238_k = addkineticlaw(rrme_reac_238,'MassAction');
set(rrme_reac_238_k, 'ParameterVariableNames', {'kon_pkc_aa'});

rrme_reac_239 = addreaction(pfpc_model, 'DAG_Ca2PKC_AA_act -> DAG_Ca2PKC_act + AA');
rrme_reac_239_k = addkineticlaw(rrme_reac_239,'MassAction');
set(rrme_reac_239_k, 'ParameterVariableNames', {'koff_pkc_aa'});

rrme_reac_240 = addreaction(pfpc_model, 'DAG_CaPKC_act + AA -> DAG_CaPKC_AA_act');
rrme_reac_240_k = addkineticlaw(rrme_reac_240,'MassAction');
set(rrme_reac_240_k, 'ParameterVariableNames', {'kon_pkc_aa'});

rrme_reac_241 = addreaction(pfpc_model, 'DAG_CaPKC_AA_act -> DAG_CaPKC_act + AA');
rrme_reac_241_k = addkineticlaw(rrme_reac_241,'MassAction');
set(rrme_reac_241_k, 'ParameterVariableNames', {'koff_pkc_aa'});

rrme_reac_242 = addreaction(pfpc_model, 'CaPKC_act + AA -> CaPKC_AA_act');
rrme_reac_242_k = addkineticlaw(rrme_reac_242,'MassAction');
set(rrme_reac_242_k, 'ParameterVariableNames', {'kon_pkc_aa'});

rrme_reac_243 = addreaction(pfpc_model, 'CaPKC_AA_act -> CaPKC_act + AA');
rrme_reac_243_k = addkineticlaw(rrme_reac_243,'MassAction');
set(rrme_reac_243_k, 'ParameterVariableNames', {'koff_pkc_aa'});

%DGK reactions
%Phosphorylation of DGK by PKC
rrme_reac_244 = addreaction(pfpc_model, 'CaPKC_act + DGK -> CaPKC_act_DGK');
rrme_reac_244_k = addkineticlaw(rrme_reac_244,'MassAction');
set(rrme_reac_244_k, 'ParameterVariableNames', {'kon_pkc_dgk'});

rrme_reac_245 = addreaction(pfpc_model, 'CaPKC_act_DGK -> CaPKC_act + DGK');
rrme_reac_245_k = addkineticlaw(rrme_reac_245,'MassAction');
set(rrme_reac_245_k, 'ParameterVariableNames', {'koff_pkc_dgk'});

rrme_reac_246 = addreaction(pfpc_model, 'CaPKC_act_DGK -> CaPKC_act + pDGK');
rrme_reac_246_k = addkineticlaw(rrme_reac_246,'MassAction');
set(rrme_reac_246_k, 'ParameterVariableNames', {'kcat_pkc_dgk'});


rrme_reac_247 = addreaction(pfpc_model, 'Ca2PKC_act + DGK -> Ca2PKC_act_DGK');
rrme_reac_247_k = addkineticlaw(rrme_reac_247,'MassAction');
set(rrme_reac_247_k, 'ParameterVariableNames', {'kon_pkc_dgk'});

rrme_reac_248 = addreaction(pfpc_model, 'Ca2PKC_act_DGK -> Ca2PKC_act + DGK');
rrme_reac_248_k = addkineticlaw(rrme_reac_248,'MassAction');
set(rrme_reac_248_k, 'ParameterVariableNames', {'koff_pkc_dgk'});

rrme_reac_249 = addreaction(pfpc_model, 'Ca2PKC_act_DGK -> Ca2PKC_act + pDGK');
rrme_reac_249_k = addkineticlaw(rrme_reac_249,'MassAction');
set(rrme_reac_249_k, 'ParameterVariableNames', {'kcat_pkc_dgk'});


rrme_reac_250 = addreaction(pfpc_model, 'CaPKC_AA_act + DGK -> CaPKC_AA_act_DGK');
rrme_reac_250_k = addkineticlaw(rrme_reac_250,'MassAction');
set(rrme_reac_250_k, 'ParameterVariableNames', {'kon_pkc_dgk'});

rrme_reac_251 = addreaction(pfpc_model, 'CaPKC_AA_act_DGK -> CaPKC_AA_act + DGK');
rrme_reac_251_k = addkineticlaw(rrme_reac_251,'MassAction');
set(rrme_reac_251_k, 'ParameterVariableNames', {'koff_pkc_dgk'});

rrme_reac_252 = addreaction(pfpc_model, 'CaPKC_AA_act_DGK -> CaPKC_AA_act + pDGK');
rrme_reac_252_k = addkineticlaw(rrme_reac_252,'MassAction');
set(rrme_reac_252_k, 'ParameterVariableNames', {'kcat_pkc_dgk'});


rrme_reac_253 = addreaction(pfpc_model, 'Ca2PKC_AA_act + DGK -> Ca2PKC_AA_act_DGK');
rrme_reac_253_k = addkineticlaw(rrme_reac_253,'MassAction');
set(rrme_reac_253_k, 'ParameterVariableNames', {'kon_pkc_dgk'});

rrme_reac_254 = addreaction(pfpc_model, 'Ca2PKC_AA_act_DGK -> Ca2PKC_AA_act + DGK');
rrme_reac_254_k = addkineticlaw(rrme_reac_254,'MassAction');
set(rrme_reac_254_k, 'ParameterVariableNames', {'koff_pkc_dgk'});

rrme_reac_255 = addreaction(pfpc_model, 'Ca2PKC_AA_act_DGK -> Ca2PKC_AA_act + pDGK');
rrme_reac_255_k = addkineticlaw(rrme_reac_255,'MassAction');
set(rrme_reac_255_k, 'ParameterVariableNames', {'kcat_pkc_dgk'});


rrme_reac_256 = addreaction(pfpc_model, 'DAG_CaPKC_act + DGK -> DAG_CaPKC_act_DGK');
rrme_reac_256_k = addkineticlaw(rrme_reac_256,'MassAction');
set(rrme_reac_256_k, 'ParameterVariableNames', {'kon_pkc_dgk'});

rrme_reac_257 = addreaction(pfpc_model, 'DAG_CaPKC_act_DGK -> DAG_CaPKC_act + DGK');
rrme_reac_257_k = addkineticlaw(rrme_reac_257,'MassAction');
set(rrme_reac_257_k, 'ParameterVariableNames', {'koff_pkc_dgk'});

rrme_reac_258 = addreaction(pfpc_model, 'DAG_CaPKC_act_DGK -> DAG_CaPKC_act + pDGK');
rrme_reac_258_k = addkineticlaw(rrme_reac_258,'MassAction');
set(rrme_reac_258_k, 'ParameterVariableNames', {'kcat_pkc_dgk'});


rrme_reac_259 = addreaction(pfpc_model, 'DAG_Ca2PKC_act + DGK -> DAG_Ca2PKC_act_DGK');
rrme_reac_259_k = addkineticlaw(rrme_reac_259,'MassAction');
set(rrme_reac_259_k, 'ParameterVariableNames', {'kon_pkc_dgk'});

rrme_reac_260 = addreaction(pfpc_model, 'DAG_Ca2PKC_act_DGK -> DAG_Ca2PKC_act + DGK');
rrme_reac_260_k = addkineticlaw(rrme_reac_260,'MassAction');
set(rrme_reac_260_k, 'ParameterVariableNames', {'koff_pkc_dgk'});

rrme_reac_261 = addreaction(pfpc_model, 'DAG_Ca2PKC_act_DGK -> DAG_Ca2PKC_act + pDGK');
rrme_reac_261_k = addkineticlaw(rrme_reac_261,'MassAction');
set(rrme_reac_261_k, 'ParameterVariableNames', {'kcat_pkc_dgk'});


rrme_reac_262 = addreaction(pfpc_model, 'DAG_CaPKC_AA_act + DGK -> DAG_CaPKC_AA_act_DGK');
rrme_reac_262_k = addkineticlaw(rrme_reac_262,'MassAction');
set(rrme_reac_262_k, 'ParameterVariableNames', {'kon_pkc_dgk'});

rrme_reac_263 = addreaction(pfpc_model, 'DAG_CaPKC_AA_act_DGK -> DAG_CaPKC_AA_act + DGK');
rrme_reac_263_k = addkineticlaw(rrme_reac_263,'MassAction');
set(rrme_reac_263_k, 'ParameterVariableNames', {'koff_pkc_dgk'});

rrme_reac_264 = addreaction(pfpc_model, 'DAG_CaPKC_AA_act_DGK -> DAG_CaPKC_AA_act + pDGK');
rrme_reac_264_k = addkineticlaw(rrme_reac_264,'MassAction');
set(rrme_reac_264_k, 'ParameterVariableNames', {'kcat_pkc_dgk'});


rrme_reac_265 = addreaction(pfpc_model, 'DAG_Ca2PKC_AA_act + DGK -> DAG_Ca2PKC_AA_act_DGK');
rrme_reac_265_k = addkineticlaw(rrme_reac_265,'MassAction');
set(rrme_reac_265_k, 'ParameterVariableNames', {'kon_pkc_dgk'});

rrme_reac_266 = addreaction(pfpc_model, 'DAG_Ca2PKC_AA_act_DGK -> DAG_Ca2PKC_AA_act + DGK');
rrme_reac_266_k = addkineticlaw(rrme_reac_266,'MassAction');
set(rrme_reac_266_k, 'ParameterVariableNames', {'koff_pkc_dgk'});

rrme_reac_267 = addreaction(pfpc_model, 'DAG_Ca2PKC_AA_act_DGK -> DAG_Ca2PKC_AA_act + pDGK');
rrme_reac_267_k = addkineticlaw(rrme_reac_267,'MassAction');
set(rrme_reac_267_k, 'ParameterVariableNames', {'kcat_pkc_dgk'});

%Breakdown of DAG by DGK
rrme_reac_268 = addreaction(pfpc_model, 'DGK + DAG -> DGK_DAG');
rrme_reac_268_k = addkineticlaw(rrme_reac_268,'MassAction');
set(rrme_reac_268_k, 'ParameterVariableNames', {'kon_dgk_dag'});

rrme_reac_269 = addreaction(pfpc_model, 'DGK_DAG -> DGK + DAG');
rrme_reac_269_k = addkineticlaw(rrme_reac_269,'MassAction');
set(rrme_reac_269_k, 'ParameterVariableNames', {'koff_dgk_dag'});

rrme_reac_270 = addreaction(pfpc_model, 'DGK_DAG -> DGK');
rrme_reac_270_k = addkineticlaw(rrme_reac_270,'MassAction');
set(rrme_reac_270_k, 'ParameterVariableNames', {'kcat_dgk'});


%Dephosp. of pDGK by PP2A
rrme_reac_271 = addreaction(pfpc_model, 'PP2A + pDGK -> PP2A_pDGK');
rrme_reac_271_k = addkineticlaw(rrme_reac_271,'MassAction');
set(rrme_reac_271_k, 'ParameterVariableNames', {'kon_pp2a_dgk'});

rrme_reac_272 = addreaction(pfpc_model, 'PP2A_pDGK -> PP2A + pDGK');
rrme_reac_272_k = addkineticlaw(rrme_reac_272,'MassAction');
set(rrme_reac_272_k, 'ParameterVariableNames', {'koff_pp2a_dgk'});

rrme_reac_273 = addreaction(pfpc_model, 'PP2A_pDGK -> PP2A + DGK');
rrme_reac_273_k = addkineticlaw(rrme_reac_273,'MassAction');
set(rrme_reac_273_k, 'ParameterVariableNames', {'kcat_pp2a_dgk'});




%%%%NEW AC REACTIONS%%%%%

AC = addspecies(cytosol, 'AC', 'InitialAmount', 0.1, 'InitialAmountUnits', 'micromolarity');
ACact = addspecies(cytosol, 'ACact', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
PDE4 = addspecies(cytosol, 'PDE4', 'InitialAmount', 0.2, 'InitialAmountUnits', 'micromolarity');
cAMP = addspecies(cytosol, 'cAMP', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
PDE4_cAMP = addspecies(cytosol, 'PDE4_cAMP', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');

kon_ac_ga = addparameter(pfpc_model,'kon_ac_ga', 'Value', 1, 'ValueUnits', '1/(micromolarity*second)');
koff_ac_ga = addparameter(pfpc_model,'koff_ac_ga', 'Value', 0.02, 'ValueUnits', '1/second');
kcat_ac = addparameter(pfpc_model,'kcat_ac', 'Value', 12, 'ValueUnits', '1/second');

kon_pde_camp = addparameter(pfpc_model,'kon_pde_camp', 'Value', 2.5, 'ValueUnits', '1/(micromolarity*second)');
koff_pde_camp = addparameter(pfpc_model,'koff_pde_camp', 'Value', 0.4, 'ValueUnits', '1/second');
kcat_pde4 = addparameter(pfpc_model,'kcat_pde4', 'Value', 0.5, 'ValueUnits', '1/second');

mglur_reac_30 = addreaction(pfpc_model, 'AC + ga_gtp -> ACact');
mglur_reac_30_k = addkineticlaw(mglur_reac_30,'MassAction');
set(mglur_reac_30_k, 'ParameterVariableNames', {'kon_ac_ga'});

mglur_reac_31 = addreaction(pfpc_model, 'ACact -> AC + ga_gtp');
mglur_reac_31_k = addkineticlaw(mglur_reac_31,'MassAction');
set(mglur_reac_31_k, 'ParameterVariableNames', {'koff_ac_ga'});

mglur_reac_32 = addreaction(pfpc_model, 'ACact -> ACact + cAMP');
mglur_reac_32_k = addkineticlaw(mglur_reac_32,'MassAction');
set(mglur_reac_32_k, 'ParameterVariableNames', {'kcat_ac'});

mglur_reac_33a = addreaction(pfpc_model, 'PDE4 + cAMP -> PDE4_cAMP');
mglur_reac_33a_k = addkineticlaw(mglur_reac_33a,'MassAction');
set(mglur_reac_33a_k, 'ParameterVariableNames', {'kon_pde_camp'});

mglur_reac_33b = addreaction(pfpc_model, 'PDE4_cAMP -> PDE4 + cAMP');
mglur_reac_33b_k = addkineticlaw(mglur_reac_33b,'MassAction');
set(mglur_reac_33b_k, 'ParameterVariableNames', {'koff_pde_camp'});

mglur_reac_33c = addreaction(pfpc_model, 'PDE4_cAMP -> PDE4');
mglur_reac_33c_k = addkineticlaw(mglur_reac_33c,'MassAction');
set(mglur_reac_33c_k, 'ParameterVariableNames', {'kcat_pde4'});

%%%%%AMPAR Trafficking Model%%%%%%%
%Species
ampar = addspecies(cytosol, 'ampar', 'InitialAmount', 0.5, 'InitialAmountUnits', 'micromolarity'); 
ampar_x = addspecies(cytosol, 'ampar_x', 'InitialAmount', 0.5, 'InitialAmountUnits', 'micromolarity');
amparp = addspecies(cytosol, 'amparp', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
amparp_x = addspecies(cytosol, 'amparp_x', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
pick = addspecies(cytosol, 'pick', 'InitialAmount', 0.7, 'InitialAmountUnits', 'micromolarity');
ampar_pick = addspecies(cytosol, 'ampar_pick', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
amparp_pick = addspecies(cytosol, 'amparp_pick', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
ampar_s = addspecies(cytosol, 'ampar_s', 'InitialAmount', 0.15, 'InitialAmountUnits', 'micromolarity');
ampar_s_pick = addspecies(cytosol, 'ampar_s_pick', 'InitialAmount', 0.05, 'InitialAmountUnits', 'micromolarity');
amparp_s = addspecies(cytosol, 'amparp_s', 'InitialAmount', 0.15, 'InitialAmountUnits', 'micromolarity');

amparp_s_pick = addspecies(cytosol, 'amparp_s_pick', 'InitialAmount', 0.05, 'InitialAmountUnits', 'micromolarity');
ampar_ez = addspecies(cytosol, 'ampar_ez', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
ampar_ez_pick = addspecies(cytosol, 'ampar_ez_pick', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
amparp_ez = addspecies(cytosol, 'amparp_ez', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
amparp_ez_pick = addspecies(cytosol, 'amparp_ez_pick', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
ampar_c = addspecies(cytosol, 'ampar_c', 'InitialAmount', 0.25, 'InitialAmountUnits', 'micromolarity');
amparp_c = addspecies(cytosol, 'amparp_c', 'InitialAmount', 0.25, 'InitialAmountUnits', 'micromolarity');
ampar_CaPKC_act = addspecies(cytosol, 'ampar_CaPKC_act', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
ampar_x_CaPKC_act = addspecies(cytosol, 'ampar_x_CaPKC_act', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
ampar_pick_CaPKC_act = addspecies(cytosol, 'ampar_pick_CaPKC_act', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');

ampar_s_CaPKC_act = addspecies(cytosol, 'ampar_s_CaPKC_act', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
ampar_s_pick_CaPKC_act = addspecies(cytosol, 'ampar_s_pick_CaPKC_act', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
ampar_ez_CaPKC_act = addspecies(cytosol, 'ampar_ez_CaPKC_act', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
ampar_ez_pick_CaPKC_act = addspecies(cytosol, 'ampar_ez_pick_CaPKC_act', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
ampar_Ca2PKC_act = addspecies(cytosol, 'ampar_Ca2PKC_act', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
ampar_x_Ca2PKC_act = addspecies(cytosol, 'ampar_x_Ca2PKC_act', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
ampar_pick_Ca2PKC_act = addspecies(cytosol, 'ampar_pick_Ca2PKC_act', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
ampar_s_Ca2PKC_act = addspecies(cytosol, 'ampar_s_Ca2PKC_act', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
ampar_s_pick_Ca2PKC_act = addspecies(cytosol, 'ampar_s_pick_Ca2PKC_act', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
ampar_ez_Ca2PKC_act = addspecies(cytosol, 'ampar_ez_Ca2PKC_act', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');

ampar_ez_pick_Ca2PKC_act = addspecies(cytosol, 'ampar_ez_pick_Ca2PKC_act', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
ampar_CaPKC_AA_act = addspecies(cytosol, 'ampar_CaPKC_AA_act', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
ampar_x_CaPKC_AA_act = addspecies(cytosol, 'ampar_x_CaPKC_AA_act', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
ampar_pick_CaPKC_AA_act = addspecies(cytosol, 'ampar_pick_CaPKC_AA_act', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
ampar_s_CaPKC_AA_act = addspecies(cytosol, 'ampar_s_CaPKC_AA_act', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
ampar_s_pick_CaPKC_AA_act = addspecies(cytosol, 'ampar_s_pick_CaPKC_AA_act', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
ampar_ez_CaPKC_AA_act = addspecies(cytosol, 'ampar_ez_CaPKC_AA_act', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
ampar_ez_pick_CaPKC_AA_act = addspecies(cytosol, 'ampar_ez_pick_CaPKC_AA_act', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
ampar_Ca2PKC_AA_act = addspecies(cytosol, 'ampar_Ca2PKC_AA_act', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
ampar_x_Ca2PKC_AA_act = addspecies(cytosol, 'ampar_x_Ca2PKC_AA_act', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');

ampar_pick_Ca2PKC_AA_act  = addspecies(cytosol, 'ampar_pick_Ca2PKC_AA_act ', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
ampar_s_Ca2PKC_AA_act = addspecies(cytosol, 'ampar_s_Ca2PKC_AA_act', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
ampar_s_pick_Ca2PKC_AA_act = addspecies(cytosol, 'ampar_s_pick_Ca2PKC_AA_act', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
ampar_ez_Ca2PKC_AA_act = addspecies(cytosol, 'ampar_ez_Ca2PKC_AA_act', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
ampar_ez_pick_Ca2PKC_AA_act = addspecies(cytosol, 'ampar_ez_pick_Ca2PKC_AA_act', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
ampar_DAG_CaPKC_act = addspecies(cytosol, 'ampar_DAG_CaPKC_act', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
ampar_x_DAG_CaPKC_act = addspecies(cytosol, 'ampar_x_DAG_CaPKC_act', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
ampar_pick_DAG_CaPKC_act = addspecies(cytosol, 'ampar_pick_DAG_CaPKC_act', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
ampar_s_DAG_CaPKC_act = addspecies(cytosol, 'ampar_s_DAG_CaPKC_act', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
ampar_s_pick_DAG_CaPKC_act = addspecies(cytosol, 'ampar_s_pick_DAG_CaPKC_act', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');

ampar_ez_DAG_CaPKC_act = addspecies(cytosol, 'ampar_ez_DAG_CaPKC_act', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
ampar_ez_pick_DAG_CaPKC_act = addspecies(cytosol, 'ampar_ez_pick_DAG_CaPKC_act', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
ampar_DAG_Ca2PKC_act = addspecies(cytosol, 'ampar_DAG_Ca2PKC_act', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
ampar_x_DAG_Ca2PKC_act = addspecies(cytosol, 'ampar_x_DAG_Ca2PKC_act', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
ampar_pick_DAG_Ca2PKC_act = addspecies(cytosol, 'ampar_pick_DAG_Ca2PKC_act', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
ampar_s_DAG_Ca2PKC_act = addspecies(cytosol, 'ampar_s_DAG_Ca2PKC_act', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
ampar_s_pick_DAG_Ca2PKC_act = addspecies(cytosol, 'ampar_s_pick_DAG_Ca2PKC_act', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
ampar_ez_DAG_Ca2PKC_act = addspecies(cytosol, 'ampar_ez_DAG_Ca2PKC_act', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
ampar_ez_pick_DAG_Ca2PKC_act = addspecies(cytosol, 'ampar_ez_pick_DAG_Ca2PKC_act', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
ampar_DAG_CaPKC_AA_act = addspecies(cytosol, 'ampar_DAG_CaPKC_AA_act', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');

ampar_x_DAG_CaPKC_AA_act = addspecies(cytosol, 'ampar_x_DAG_CaPKC_AA_act', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
ampar_pick_DAG_CaPKC_AA_act = addspecies(cytosol, 'ampar_pick_DAG_CaPKC_AA_act', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
ampar_s_DAG_CaPKC_AA_act = addspecies(cytosol, 'ampar_s_DAG_CaPKC_AA_act', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
ampar_s_pick_DAG_CaPKC_AA_act = addspecies(cytosol, 'ampar_s_pick_DAG_CaPKC_AA_act', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
ampar_ez_DAG_CaPKC_AA_act = addspecies(cytosol, 'ampar_ez_DAG_CaPKC_AA_act', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
ampar_ez_pick_DAG_CaPKC_AA_act = addspecies(cytosol, 'ampar_ez_pick_DAG_CaPKC_AA_act', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
ampar_DAG_Ca2PKC_AA_act = addspecies(cytosol, 'ampar_DAG_Ca2PKC_AA_act', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
ampar_x_DAG_Ca2PKC_AA_act = addspecies(cytosol, 'ampar_x_DAG_Ca2PKC_AA_act', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
ampar_pick_DAG_Ca2PKC_AA_act = addspecies(cytosol, 'ampar_pick_DAG_Ca2PKC_AA_act', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
ampar_s_DAG_Ca2PKC_AA_act = addspecies(cytosol, 'ampar_s_DAG_Ca2PKC_AA_act', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');

ampar_s_pick_DAG_Ca2PKC_AA_act = addspecies(cytosol, 'ampar_s_pick_DAG_Ca2PKC_AA_act', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
ampar_ez_DAG_Ca2PKC_AA_act = addspecies(cytosol, 'ampar_ez_DAG_Ca2PKC_AA_act', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
ampar_ez_pick_DAG_Ca2PKC_AA_act = addspecies(cytosol, 'ampar_ez_pick_DAG_Ca2PKC_AA_act', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
amparp_pp2a = addspecies(cytosol, 'amparp_pp2a', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
amparp_pick_pp2a = addspecies(cytosol, 'amparp_pick_pp2a', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
amparp_s_pp2a = addspecies(cytosol, 'amparp_s_pp2a', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
amparp_s_pick_pp2a = addspecies(cytosol, 'amparp_s_pick_pp2a', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
amparp_ez_pp2a = addspecies(cytosol, 'amparp_ez_pp2a', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
amparp_ez_pick_pp2a = addspecies(cytosol, 'amparp_ez_pick_pp2a', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
amparp_pp1 = addspecies(cytosol, 'amparp_pp1', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
amparp_pick_pp1 = addspecies(cytosol, 'amparp_pick_pp1', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
amparp_s_pp1 = addspecies(cytosol, 'amparp_s_pp1', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
amparp_s_pick_pp1 = addspecies(cytosol, 'amparp_s_pick_pp1', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
amparp_ez_pp1 = addspecies(cytosol, 'amparp_ez_pp1', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
amparp_ez_pick_pp1 = addspecies(cytosol, 'amparp_ez_pick_pp1', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
nsf = addspecies(cytosol, 'nsf', 'InitialAmount', 4, 'InitialAmountUnits', 'micromolarity');
nsfp = addspecies(cytosol, 'nsfp', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');

ampar_pick_nsf = addspecies(cytosol, 'ampar_pick_nsf', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
amparp_pick_nsf = addspecies(cytosol, 'amparp_pick_nsf', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
ampar_s_pick_nsf = addspecies(cytosol, 'ampar_s_pick_nsf', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
amparp_s_pick_nsf = addspecies(cytosol, 'amparp_s_pick_nsf', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
ampar_ez_pick_nsf = addspecies(cytosol, 'ampar_ez_pick_nsf', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
amparp_ez_pick_nsf = addspecies(cytosol, 'amparp_ez_pick_nsf', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
ampar_c_pick_nsf = addspecies(cytosol, 'ampar_c_pick_nsf', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
amparp_c_pick_nsf = addspecies(cytosol, 'amparp_c_pick_nsf', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
ampar_c_pick = addspecies(cytosol, 'ampar_c_pick', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
amparp_c_pick = addspecies(cytosol, 'amparp_c_pick', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');

ampar_psd = addspecies(cytosol, 'ampar_psd', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
ampar_cytosol = addspecies(cytosol, 'ampar_cytosol', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
ampar_endozone = addspecies(cytosol, 'ampar_endozone', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
ampar_synapse = addspecies(cytosol, 'ampar_synapse', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%NSF Nitrosylation%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%nsf = addspecies(cytosol, 'nsf', 'InitialAmount', 2, 'InitialAmountUnits', 'micromolarity');
NSF_NO_tot = addspecies(cytosol, 'NSF_NO_tot', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
NSF_NO3_tot = addspecies(cytosol, 'NSF_NO3_tot', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
NSFP_NO3_tot = addspecies(cytosol, 'NSFP_NO3_tot', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
nsf_no = addspecies(cytosol, 'nsf_0no', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
nsf_2no = addspecies(cytosol, 'nsf_2no', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
nsf_3no = addspecies(cytosol, 'nsf_3no', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');

nsfp_no = addspecies(cytosol, 'nsfp_no', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
nsfp_2no = addspecies(cytosol, 'nsfp_2no', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
nsfp_3no = addspecies(cytosol, 'nsfp_3no', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');

kon_nsf_no = addparameter(pfpc_model,'kon_nsf_no', 'Value', 4, 'ValueUnits', '1/(micromolarity*second)', 'ConstantValue', false);
koff_nsf_no = addparameter(pfpc_model,'koff_nsf_no', 'Value', 17, 'ValueUnits', '1/second');
kon_nsf_no2 = addparameter(pfpc_model,'kon_nsf_no2', 'Value', 3.2, 'ValueUnits', '1/(micromolarity*second)', 'ConstantValue', false);
koff_nsf_no2 = addparameter(pfpc_model,'koff_nsf_no2', 'Value', 7, 'ValueUnits', '1/second','ConstantValue', false);%7
kon_nsf_no3 = addparameter(pfpc_model,'kon_nsf_no3', 'Value', 2, 'ValueUnits', '1/(micromolarity*second)','ConstantValue', false);
koff_nsf_no3 = addparameter(pfpc_model,'koff_nsf_no3', 'Value', 0.00015, 'ValueUnits', '1/second', 'ConstantValue', false);

%D3
nsf_reac_5 = addreaction(pfpc_model, 'nsf + NO -> nsf_no')
nsf_reac_5_k = addkineticlaw(nsf_reac_5,'MassAction');
set(nsf_reac_5_k, 'ParameterVariableNames', {'kon_nsf_no'});
nsf_reac_6 = addreaction(pfpc_model, 'nsf_no -> nsf + NO')
nsf_reac_6_k = addkineticlaw(nsf_reac_6,'MassAction');
set(nsf_reac_6_k, 'ParameterVariableNames', {'koff_nsf_no'});
%D4
nsf_reac_7 = addreaction(pfpc_model, 'nsf_no + NO -> nsf_2no')
nsf_reac_7_k = addkineticlaw(nsf_reac_7,'MassAction');
set(nsf_reac_7_k, 'ParameterVariableNames', {'kon_nsf_no2'});
nsf_reac_8 = addreaction(pfpc_model, 'nsf_2no -> nsf_no + NO')
nsf_reac_8_k = addkineticlaw(nsf_reac_8,'MassAction');
set(nsf_reac_8_k, 'ParameterVariableNames', {'koff_nsf_no2'});
%D5
nsf_reac_9 = addreaction(pfpc_model, 'nsf_2no + NO -> nsf_3no')
nsf_reac_9_k = addkineticlaw(nsf_reac_9,'MassAction');
set(nsf_reac_9_k, 'ParameterVariableNames', {'kon_nsf_no3'});
nsf_reac_10 = addreaction(pfpc_model, 'nsf_3no -> nsf_2no + NO')
nsf_reac_10_k = addkineticlaw(nsf_reac_10,'MassAction');
set(nsf_reac_10_k, 'ParameterVariableNames', {'koff_nsf_no3'});


ampar_pick_nsf_no = addspecies(cytosol, 'ampar_pick_nsf_no', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
amparp_pick_nsf_no = addspecies(cytosol, 'amparp_pick_nsf_no', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
ampar_s_pick_nsf_no = addspecies(cytosol, 'ampar_s_pick_nsf_no', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
amparp_s_pick_nsf_no = addspecies(cytosol, 'amparp_s_pick_nsf_no', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
ampar_ez_pick_nsf_no = addspecies(cytosol, 'ampar_ez_pick_nsf_no', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
amparp_ez_pick_nsf_no = addspecies(cytosol, 'amparp_ez_pick_nsf_no', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
ampar_c_pick_nsf_no = addspecies(cytosol, 'ampar_c_pick_nsf_no', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
amparp_c_pick_nsf_no = addspecies(cytosol, 'amparp_c_pick_nsf_no', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');

ampar_pick_nsf_2no = addspecies(cytosol, 'ampar_pick_nsf_2no', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
amparp_pick_nsf_2no = addspecies(cytosol, 'amparp_pick_nsf_2no', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
ampar_s_pick_nsf_2no = addspecies(cytosol, 'ampar_s_pick_nsf_2no', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
amparp_s_pick_nsf_2no = addspecies(cytosol, 'amparp_s_pick_nsf_2no', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
ampar_ez_pick_nsf_2no = addspecies(cytosol, 'ampar_ez_pick_nsf_2no', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
amparp_ez_pick_nsf_2no = addspecies(cytosol, 'amparp_ez_pick_nsf_2no', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
ampar_c_pick_nsf_2no = addspecies(cytosol, 'ampar_c_pick_nsf_2no', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
amparp_c_pick_nsf_2no = addspecies(cytosol, 'amparp_c_pick_nsf_2no', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');

ampar_pick_nsf_3no = addspecies(cytosol, 'ampar_pick_nsf_3no', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
amparp_pick_nsf_3no = addspecies(cytosol, 'amparp_pick_nsf_3no', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
ampar_s_pick_nsf_3no = addspecies(cytosol, 'ampar_s_pick_nsf_3no', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
amparp_s_pick_nsf_3no = addspecies(cytosol, 'amparp_s_pick_nsf_3no', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
ampar_ez_pick_nsf_3no = addspecies(cytosol, 'ampar_ez_pick_nsf_3no', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
amparp_ez_pick_nsf_3no = addspecies(cytosol, 'amparp_ez_pick_nsf_3no', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
ampar_c_pick_nsf_3no = addspecies(cytosol, 'ampar_c_pick_nsf_3no', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
amparp_c_pick_nsf_3no = addspecies(cytosol, 'amparp_c_pick_nsf_3no', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');

%Parameters
kon_psd = addparameter(pfpc_model,'kon_psd', 'Value', 5.5, 'ValueUnits', '1/second');
koff_psd = addparameter(pfpc_model,'koff_psd', 'Value', 0.3, 'ValueUnits', '1/second');
kon_ampar_pick = addparameter(pfpc_model,'kon_ampar_pick', 'Value', 25, 'ValueUnits', '1/(micromolarity*second)','ConstantValue', false);
koff_ampar_pick = addparameter(pfpc_model,'koff_ampar_pick', 'Value',5, 'ValueUnits', '1/second');
koff_psd_s880 = addparameter(pfpc_model,'koff_psd_s880', 'Value', 70, 'ValueUnits', '1/second');
kendo = addparameter(pfpc_model,'kendo', 'Value', 0.1, 'ValueUnits', '1/second');
kexo = addparameter(pfpc_model,'kexo', 'Value', 0.0002, 'ValueUnits', '1/second');
kon_ampar_pkc = addparameter(pfpc_model,'kon_ampar_pkc', 'Value', 0.1, 'ValueUnits', '1/(micromolarity*second)');
koff_ampar_pkc = addparameter(pfpc_model,'koff_ampar_pkc', 'Value', 1, 'ValueUnits', '1/second');
kcat_pkc_ampar = addparameter(pfpc_model,'kcat_pkc_ampar', 'Value', 4.7, 'ValueUnits', '1/second');

kon_ampar_pp1 = addparameter(pfpc_model,'kon_ampar_pp1', 'Value', 1.5, 'ValueUnits', '1/(micromolarity*second)', 'ConstantValue', false);
koff_ampar_pp1 = addparameter(pfpc_model,'koff_ampar_pp1', 'Value', 10, 'ValueUnits', '1/second');
kcat_pp1_ampar = addparameter(pfpc_model,'kcat_pp1_ampar', 'Value', 0.1, 'ValueUnits', '1/second');

kon_ampar_pp2a = addparameter(pfpc_model,'kon_ampar_pp2a', 'Value', 1.5, 'ValueUnits', '1/(micromolarity*second)', 'ConstantValue', false);
koff_ampar_pp2a = addparameter(pfpc_model,'koff_ampar_pp2a', 'Value', 5, 'ValueUnits', '1/second');
kcat_pp2a_ampar = addparameter(pfpc_model,'kcat_pp2a_ampar', 'Value', 0.15, 'ValueUnits', '1/second');
kdiff_psd_syn = addparameter(pfpc_model,'kdiff_psd_syn', 'Value', 0.01, 'ValueUnits', '1/second');
kdiff_syn_psd = addparameter(pfpc_model,'kdiff_syn_psd', 'Value', 0.02, 'ValueUnits', '1/second');
kdiff_syn_ez = addparameter(pfpc_model,'kdiff_syn_ez', 'Value', 0.02, 'ValueUnits', '1/second');
kdiff_ez_syn = addparameter(pfpc_model,'kdiff_ez_syn', 'Value', 0.01, 'ValueUnits', '1/second');
kon_ampar_nsf = addparameter(pfpc_model,'kon_ampar_nsf', 'Value', 106, 'ValueUnits', '1/(micromolarity*second)','ConstantValue', false);
kon_ampar_nsf_no = addparameter(pfpc_model,'kon_ampar_nsf_no', 'Value', 2000, 'ValueUnits', '1/(micromolarity*second)','ConstantValue', false);
koff_ampar_nsf = addparameter(pfpc_model,'koff_ampar_nsf', 'Value', 1, 'ValueUnits', '1/second');
kcat_nsf = addparameter(pfpc_model,'kcat_nsf', 'Value', 2, 'ValueUnits', '1/second');
kcat_nsf_no = addparameter(pfpc_model,'kcat_nsf_no', 'Value', 100, 'ValueUnits', '1/second');



%Reactions

ampar_reac_1 = addreaction(pfpc_model, 'ampar -> ampar_x');
ampar_reac_1_k = addkineticlaw(ampar_reac_1,'MassAction');
set(ampar_reac_1_k, 'ParameterVariableNames', {'kon_psd'});

ampar_reac_2 = addreaction(pfpc_model, 'ampar_x -> ampar');
ampar_reac_2_k = addkineticlaw(ampar_reac_2,'MassAction');
set(ampar_reac_2_k, 'ParameterVariableNames', {'koff_psd'});

ampar_reac_3 = addreaction(pfpc_model, 'amparp -> amparp_x');
ampar_reac_3_k = addkineticlaw(ampar_reac_3,'MassAction');
set(ampar_reac_3_k, 'ParameterVariableNames', {'kon_psd'});

ampar_reac_4 = addreaction(pfpc_model, 'amparp_x -> amparp');
ampar_reac_4_k = addkineticlaw(ampar_reac_4,'MassAction');
set(ampar_reac_4_k, 'ParameterVariableNames', {'koff_psd_s880'});

ampar_reac_5 = addreaction(pfpc_model, 'ampar + pick -> ampar_pick');
ampar_reac_5_k = addkineticlaw(ampar_reac_5,'MassAction');
set(ampar_reac_5_k, 'ParameterVariableNames', {'kon_ampar_pick'});

ampar_reac_6 = addreaction(pfpc_model, 'ampar_pick -> ampar + pick');
ampar_reac_6_k = addkineticlaw(ampar_reac_6,'MassAction');
set(ampar_reac_6_k, 'ParameterVariableNames', {'koff_ampar_pick'});

ampar_reac_7 = addreaction(pfpc_model, 'amparp + pick -> amparp_pick');
ampar_reac_7_k = addkineticlaw(ampar_reac_7,'MassAction');
set(ampar_reac_7_k, 'ParameterVariableNames', {'kon_ampar_pick'});

ampar_reac_8 = addreaction(pfpc_model, 'amparp_pick -> amparp + pick');
ampar_reac_8_k = addkineticlaw(ampar_reac_8,'MassAction');
set(ampar_reac_8_k, 'ParameterVariableNames', {'koff_ampar_pick'});

ampar_reac_9 = addreaction(pfpc_model, 'ampar_s + pick -> ampar_s_pick');
ampar_reac_9_k = addkineticlaw(ampar_reac_9,'MassAction');
set(ampar_reac_9_k, 'ParameterVariableNames', {'kon_ampar_pick'});

ampar_reac_10 = addreaction(pfpc_model, 'ampar_s_pick -> ampar_s + pick');
ampar_reac_10_k = addkineticlaw(ampar_reac_10,'MassAction');
set(ampar_reac_10_k, 'ParameterVariableNames', {'koff_ampar_pick'});

ampar_reac_11 = addreaction(pfpc_model, 'amparp_s + pick -> amparp_s_pick');
ampar_reac_11_k = addkineticlaw(ampar_reac_11,'MassAction');
set(ampar_reac_11_k, 'ParameterVariableNames', {'kon_ampar_pick'});

ampar_reac_12 = addreaction(pfpc_model, 'amparp_s_pick -> amparp_s + pick');
ampar_reac_12_k = addkineticlaw(ampar_reac_12,'MassAction');
set(ampar_reac_12_k, 'ParameterVariableNames', {'koff_ampar_pick'});

ampar_reac_13 = addreaction(pfpc_model, 'ampar_ez + pick -> ampar_ez_pick');
ampar_reac_13_k = addkineticlaw(ampar_reac_13,'MassAction');
set(ampar_reac_13_k, 'ParameterVariableNames', {'kon_ampar_pick'});

ampar_reac_14 = addreaction(pfpc_model, 'ampar_ez_pick -> ampar_ez + pick');
ampar_reac_14_k = addkineticlaw(ampar_reac_14,'MassAction');
set(ampar_reac_14_k, 'ParameterVariableNames', {'koff_ampar_pick'});

ampar_reac_15 = addreaction(pfpc_model, 'amparp_ez + pick -> amparp_ez_pick');
ampar_reac_15_k = addkineticlaw(ampar_reac_15,'MassAction');
set(ampar_reac_15_k, 'ParameterVariableNames', {'kon_ampar_pick'});

ampar_reac_16 = addreaction(pfpc_model, 'amparp_ez_pick -> amparp_ez + pick');
ampar_reac_16_k = addkineticlaw(ampar_reac_16,'MassAction');
set(ampar_reac_16_k, 'ParameterVariableNames', {'koff_ampar_pick'});

ampar_reac_17 = addreaction(pfpc_model, 'ampar_ez_pick -> ampar_c_pick');
ampar_reac_17_k = addkineticlaw(ampar_reac_17,'MassAction');
set(ampar_reac_17_k, 'ParameterVariableNames', {'kendo'});

ampar_reac_18 = addreaction(pfpc_model, 'amparp_ez_pick -> amparp_c_pick');
ampar_reac_18_k = addkineticlaw(ampar_reac_18,'MassAction');
set(ampar_reac_18_k, 'ParameterVariableNames', {'kendo'});

%%%%%%NEW NSF MODEL%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ampar_sx = addspecies(cytosol, 'ampar_sx', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
amparp_sx = addspecies(cytosol, 'amparp_sx', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');

ampar_sx_nsf = addspecies(cytosol, 'ampar_sx_nsf', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
amparp_sx_nsf = addspecies(cytosol, 'amparp_sx_nsf', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');

ampar_sx_nsf_no = addspecies(cytosol, 'ampar_sx_nsf_no', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
amparp_sx_nsf_no = addspecies(cytosol, 'amparp_sx_nsf_no', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');

ampar_sx_nsf_2no = addspecies(cytosol, 'ampar_sx_nsf_2no', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
amparp_sx_nsf_2no = addspecies(cytosol, 'amparp_sx_nsf_2no', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');

ampar_sx_nsf_3no = addspecies(cytosol, 'ampar_sx_nsf_3no', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
amparp_sx_nsf_3no = addspecies(cytosol, 'amparp_sx_nsf_3no', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');

koff_ampar_scaff = addparameter(pfpc_model,'koff_ampar_scaff', 'Value', 0.001, 'ValueUnits', '1/second');

kon_ampar_nsf_no_scaff = addparameter(pfpc_model,'kon_ampar_nsf_no_scaff', 'Value', 1000, 'ValueUnits', '1/(micromolarity*second)');



ampar_reac_19 = addreaction(pfpc_model, 'ampar_c -> ampar_sx'); %%%%%%% 'ampar_c -> ampar'
ampar_reac_19_k = addkineticlaw(ampar_reac_19,'MassAction');
set(ampar_reac_19_k, 'ParameterVariableNames', {'kexo'});

ampar_reac_20 = addreaction(pfpc_model, 'amparp_c -> amparp_sx');%%%%%%% 'amparp_c -> amparp'
ampar_reac_20_k = addkineticlaw(ampar_reac_20,'MassAction');
set(ampar_reac_20_k, 'ParameterVariableNames', {'kexo'});

%slow spontaneous dissociation of AMPARs from scaffold

ampar_reac_20a = addreaction(pfpc_model, 'ampar_sx -> ampar_s');%%%%%%% 'amparp_c -> amparp'
ampar_reac_20a_k = addkineticlaw(ampar_reac_20a,'MassAction');
set(ampar_reac_20a_k, 'ParameterVariableNames', {'koff_ampar_scaff'});

ampar_reac_20b = addreaction(pfpc_model, 'amparp_sx -> amparp_s');%%%%%%% 'amparp_c -> amparp'
ampar_reac_20b_k = addkineticlaw(ampar_reac_20b,'MassAction');
set(ampar_reac_20b_k, 'ParameterVariableNames', {'koff_ampar_scaff'});




%fast NSF-mediated dissociation of AMPARs from scaffold

%NSF
ampar_reac_20c = addreaction(pfpc_model, 'ampar_sx + nsf -> ampar_sx_nsf');%%%%%%% 'amparp_c -> amparp'
ampar_reac_20c_k = addkineticlaw(ampar_reac_20c,'MassAction');
set(ampar_reac_20c_k, 'ParameterVariableNames', {'kon_ampar_nsf'});

ampar_reac_20d = addreaction(pfpc_model, 'amparp_sx + nsf -> amparp_sx_nsf');%%%%%%% 'amparp_c -> amparp'
ampar_reac_20d_k = addkineticlaw(ampar_reac_20d,'MassAction');
set(ampar_reac_20d_k, 'ParameterVariableNames', {'kon_ampar_nsf'});

ampar_reac_20e = addreaction(pfpc_model, 'ampar_sx_nsf -> ampar_sx + nsf');%%%%%%% 'amparp_c -> amparp'
ampar_reac_20e_k = addkineticlaw(ampar_reac_20e,'MassAction');
set(ampar_reac_20e_k, 'ParameterVariableNames', {'koff_ampar_nsf'});

ampar_reac_20f = addreaction(pfpc_model, 'amparp_sx_nsf -> amparp_sx + nsf');%%%%%%% 'amparp_c -> amparp'
ampar_reac_20f_k = addkineticlaw(ampar_reac_20f,'MassAction');
set(ampar_reac_20f_k, 'ParameterVariableNames', {'koff_ampar_nsf'});

ampar_reac_20g = addreaction(pfpc_model, 'ampar_sx_nsf -> ampar_s + nsf');%%%%%%% 'amparp_c -> amparp'
ampar_reac_20g_k = addkineticlaw(ampar_reac_20g,'MassAction');
set(ampar_reac_20g_k, 'ParameterVariableNames', {'kcat_nsf'});

ampar_reac_20h = addreaction(pfpc_model, 'amparp_sx_nsf -> amparp_s + nsf');%%%%%%% 'amparp_c -> amparp'
ampar_reac_20h_k = addkineticlaw(ampar_reac_20h,'MassAction');
set(ampar_reac_20h_k, 'ParameterVariableNames', {'kcat_nsf'});

%NSF_NO
ampar_reac_20nsfc = addreaction(pfpc_model, 'ampar_sx + nsf_no -> ampar_sx_nsf_no');%%%%%%% 'amparp_c -> amparp'
ampar_reac_20nsfc_k = addkineticlaw(ampar_reac_20nsfc,'MassAction');
set(ampar_reac_20nsfc_k, 'ParameterVariableNames', {'kon_ampar_nsf'});

ampar_reac_20nsfd = addreaction(pfpc_model, 'amparp_sx + nsf_no -> amparp_sx_nsf_no');%%%%%%% 'amparp_c -> amparp'
ampar_reac_20nsfd_k = addkineticlaw(ampar_reac_20nsfd,'MassAction');
set(ampar_reac_20nsfd_k, 'ParameterVariableNames', {'kon_ampar_nsf'});

ampar_reac_20nsfe = addreaction(pfpc_model, 'ampar_sx_nsf_no -> ampar_sx + nsf_no');%%%%%%% 'amparp_c -> amparp'
ampar_reac_20nsfe_k = addkineticlaw(ampar_reac_20nsfe,'MassAction');
set(ampar_reac_20nsfe_k, 'ParameterVariableNames', {'koff_ampar_nsf'});

ampar_reac_20nsff = addreaction(pfpc_model, 'amparp_sx_nsf_no -> amparp_sx + nsf_no');%%%%%%% 'amparp_c -> amparp'
ampar_reac_20nsff_k = addkineticlaw(ampar_reac_20nsff,'MassAction');
set(ampar_reac_20nsff_k, 'ParameterVariableNames', {'koff_ampar_nsf'});

ampar_reac_20nsfg = addreaction(pfpc_model, 'ampar_sx_nsf_no -> ampar_s + nsf_no');%%%%%%% 'amparp_c -> amparp'
ampar_reac_20nsfg_k = addkineticlaw(ampar_reac_20nsfg,'MassAction');
set(ampar_reac_20nsfg_k, 'ParameterVariableNames', {'kcat_nsf'});

ampar_reac_20nsfh = addreaction(pfpc_model, 'amparp_sx_nsf_no -> amparp_s + nsf_no');%%%%%%% 'amparp_c -> amparp'
ampar_reac_20nsfh_k = addkineticlaw(ampar_reac_20nsfh,'MassAction');
set(ampar_reac_20nsfh_k, 'ParameterVariableNames', {'kcat_nsf'});

%NSF_2NO
ampar_reac_20nsf2c = addreaction(pfpc_model, 'ampar_sx + nsf_2no -> ampar_sx_nsf_2no');%%%%%%% 'amparp_c -> amparp'
ampar_reac_20nsf2c_k = addkineticlaw(ampar_reac_20nsf2c,'MassAction');
set(ampar_reac_20nsf2c_k, 'ParameterVariableNames', {'kon_ampar_nsf'});

ampar_reac_20nsf2d = addreaction(pfpc_model, 'amparp_sx + nsf_2no -> amparp_sx_nsf_2no');%%%%%%% 'amparp_c -> amparp'
ampar_reac_20nsf2d_k = addkineticlaw(ampar_reac_20nsf2d,'MassAction');
set(ampar_reac_20nsf2d_k, 'ParameterVariableNames', {'kon_ampar_nsf'});

ampar_reac_20nsf2e = addreaction(pfpc_model, 'ampar_sx_nsf_2no -> ampar_sx + nsf_2no');%%%%%%% 'amparp_c -> amparp'
ampar_reac_20nsf2e_k = addkineticlaw(ampar_reac_20nsf2e,'MassAction');
set(ampar_reac_20nsf2e_k, 'ParameterVariableNames', {'koff_ampar_nsf'});

ampar_reac_20nsf2f = addreaction(pfpc_model, 'amparp_sx_nsf_2no -> amparp_sx + nsf_2no');%%%%%%% 'amparp_c -> amparp'
ampar_reac_20nsf2f_k = addkineticlaw(ampar_reac_20nsf2f,'MassAction');
set(ampar_reac_20nsf2f_k, 'ParameterVariableNames', {'koff_ampar_nsf'});

ampar_reac_20nsf2g = addreaction(pfpc_model, 'ampar_sx_nsf_2no -> ampar_s + nsf_2no');%%%%%%% 'amparp_c -> amparp'
ampar_reac_20nsf2g_k = addkineticlaw(ampar_reac_20nsf2g,'MassAction');
set(ampar_reac_20nsf2g_k, 'ParameterVariableNames', {'kcat_nsf'});

ampar_reac_20nsf2h = addreaction(pfpc_model, 'amparp_sx_nsf_2no -> amparp_s + nsf_2no');%%%%%%% 'amparp_c -> amparp'
ampar_reac_20nsf2h_k = addkineticlaw(ampar_reac_20nsf2h,'MassAction');
set(ampar_reac_20nsf2h_k, 'ParameterVariableNames', {'kcat_nsf'});



%faster NSF_3NO-mediated dissociation of AMPARs from scaffold
%NSF_3NO
ampar_reac_20i = addreaction(pfpc_model, 'ampar_sx + nsf_3no -> ampar_sx_nsf_3no');%%%%%%% 'amparp_c -> amparp'
ampar_reac_20i_k = addkineticlaw(ampar_reac_20i,'MassAction');
set(ampar_reac_20i_k, 'ParameterVariableNames', {'kon_ampar_nsf_no_scaff'});

ampar_reac_20j = addreaction(pfpc_model, 'amparp_sx + nsf_3no -> amparp_sx_nsf_3no');%%%%%%% 'amparp_c -> amparp'
ampar_reac_20j_k = addkineticlaw(ampar_reac_20j,'MassAction');
set(ampar_reac_20j_k, 'ParameterVariableNames', {'kon_ampar_nsf_no_scaff'});

ampar_reac_20k = addreaction(pfpc_model, 'ampar_sx_nsf_3no -> ampar_sx + nsf_3no');%%%%%%% 'amparp_c -> amparp'
ampar_reac_20k_k = addkineticlaw(ampar_reac_20k,'MassAction');
set(ampar_reac_20k_k, 'ParameterVariableNames', {'koff_ampar_nsf'});

ampar_reac_20l = addreaction(pfpc_model, 'amparp_sx_nsf_3no -> amparp_sx + nsf_3no');%%%%%%% 'amparp_c -> amparp'
ampar_reac_20l_k = addkineticlaw(ampar_reac_20l,'MassAction');
set(ampar_reac_20l_k, 'ParameterVariableNames', {'koff_ampar_nsf'});

ampar_reac_20m = addreaction(pfpc_model, 'ampar_sx_nsf_3no -> ampar_s + nsf_3no');%%%%%%% 'amparp_c -> amparp'
ampar_reac_20m_k = addkineticlaw(ampar_reac_20m,'MassAction');
set(ampar_reac_20m_k, 'ParameterVariableNames', {'kcat_nsf_no'});

ampar_reac_20n = addreaction(pfpc_model, 'amparp_sx_nsf_3no -> amparp_s + nsf_3no');%%%%%%% 'amparp_c -> amparp'
ampar_reac_20n_k = addkineticlaw(ampar_reac_20n,'MassAction');
set(ampar_reac_20n_k, 'ParameterVariableNames', {'kcat_nsf_no'});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ampar_reac_19 = addreaction(pfpc_model, 'ampar_c -> ampar');
ampar_reac_19_k = addkineticlaw(ampar_reac_19,'MassAction');
set(ampar_reac_19_k, 'ParameterVariableNames', {'kexo'});

ampar_reac_20 = addreaction(pfpc_model, 'amparp_c -> amparp');
ampar_reac_20_k = addkineticlaw(ampar_reac_20,'MassAction');
set(ampar_reac_20_k, 'ParameterVariableNames', {'kexo'});

ampar_reac_21a = addreaction(pfpc_model, 'ampar + CaPKC_act -> ampar_CaPKC_act');
ampar_reac_21a_k = addkineticlaw(ampar_reac_21a,'MassAction');
set(ampar_reac_21a_k, 'ParameterVariableNames', {'kon_ampar_pkc'});

ampar_reac_21b = addreaction(pfpc_model, 'ampar + Ca2PKC_act -> ampar_Ca2PKC_act');
ampar_reac_21b_k = addkineticlaw(ampar_reac_21b,'MassAction');
set(ampar_reac_21b_k, 'ParameterVariableNames', {'kon_ampar_pkc'});

ampar_reac_21c = addreaction(pfpc_model, 'ampar + CaPKC_AA_act -> ampar_CaPKC_AA_act');
ampar_reac_21c_k = addkineticlaw(ampar_reac_21c,'MassAction');
set(ampar_reac_21c_k, 'ParameterVariableNames', {'kon_ampar_pkc'});

ampar_reac_21d = addreaction(pfpc_model, 'ampar + Ca2PKC_AA_act -> ampar_Ca2PKC_AA_act');
ampar_reac_21d_k = addkineticlaw(ampar_reac_21d,'MassAction');
set(ampar_reac_21d_k, 'ParameterVariableNames', {'kon_ampar_pkc'});

ampar_reac_21e = addreaction(pfpc_model, 'ampar + DAG_CaPKC_act -> ampar_DAG_CaPKC_act');
ampar_reac_21e_k = addkineticlaw(ampar_reac_21e,'MassAction');
set(ampar_reac_21e_k, 'ParameterVariableNames', {'kon_ampar_pkc'});

ampar_reac_21f = addreaction(pfpc_model, 'ampar + DAG_Ca2PKC_act -> ampar_DAG_Ca2PKC_act');
ampar_reac_21f_k = addkineticlaw(ampar_reac_21f,'MassAction');
set(ampar_reac_21f_k, 'ParameterVariableNames', {'kon_ampar_pkc'});

ampar_reac_21g = addreaction(pfpc_model, 'ampar + DAG_CaPKC_AA_act -> ampar_DAG_CaPKC_AA_act');
ampar_reac_21g_k = addkineticlaw(ampar_reac_21g,'MassAction');
set(ampar_reac_21g_k, 'ParameterVariableNames', {'kon_ampar_pkc'});

ampar_reac_21h = addreaction(pfpc_model, 'ampar + DAG_Ca2PKC_AA_act -> ampar_DAG_Ca2PKC_AA_act');
ampar_reac_21h_k = addkineticlaw(ampar_reac_21h,'MassAction');
set(ampar_reac_21h_k, 'ParameterVariableNames', {'kon_ampar_pkc'});

ampar_reac_22a = addreaction(pfpc_model, 'ampar_CaPKC_act -> ampar + CaPKC_act');
ampar_reac_22a_k = addkineticlaw(ampar_reac_22a,'MassAction');
set(ampar_reac_22a_k, 'ParameterVariableNames', {'koff_ampar_pkc'});

ampar_reac_22b = addreaction(pfpc_model, 'ampar_Ca2PKC_act -> ampar + Ca2PKC_act');
ampar_reac_22b_k = addkineticlaw(ampar_reac_22b,'MassAction');
set(ampar_reac_22b_k, 'ParameterVariableNames', {'koff_ampar_pkc'});

ampar_reac_22c = addreaction(pfpc_model, 'ampar_CaPKC_AA_act -> ampar + CaPKC_AA_act');
ampar_reac_22c_k = addkineticlaw(ampar_reac_22c,'MassAction');
set(ampar_reac_22c_k, 'ParameterVariableNames', {'koff_ampar_pkc'});

ampar_reac_22d = addreaction(pfpc_model, 'ampar_Ca2PKC_AA_act -> ampar + Ca2PKC_AA_act');
ampar_reac_22d_k = addkineticlaw(ampar_reac_22d,'MassAction');
set(ampar_reac_22d_k, 'ParameterVariableNames', {'koff_ampar_pkc'});

ampar_reac_22e = addreaction(pfpc_model, 'ampar_DAG_CaPKC_act -> ampar + DAG_CaPKC_act');
ampar_reac_22e_k = addkineticlaw(ampar_reac_22e,'MassAction');
set(ampar_reac_22e_k, 'ParameterVariableNames', {'koff_ampar_pkc'});

ampar_reac_22f = addreaction(pfpc_model, 'ampar_DAG_Ca2PKC_act -> ampar + DAG_Ca2PKC_act');
ampar_reac_22f_k = addkineticlaw(ampar_reac_22f,'MassAction');
set(ampar_reac_22f_k, 'ParameterVariableNames', {'koff_ampar_pkc'});

ampar_reac_22g = addreaction(pfpc_model, 'ampar_DAG_CaPKC_AA_act -> ampar + DAG_CaPKC_AA_act');
ampar_reac_22g_k = addkineticlaw(ampar_reac_22g,'MassAction');
set(ampar_reac_22g_k, 'ParameterVariableNames', {'koff_ampar_pkc'});

ampar_reac_22h = addreaction(pfpc_model, 'ampar_DAG_Ca2PKC_AA_act -> ampar + DAG_Ca2PKC_AA_act');
ampar_reac_22h_k = addkineticlaw(ampar_reac_22h,'MassAction');
set(ampar_reac_22h_k, 'ParameterVariableNames', {'koff_ampar_pkc'});

ampar_reac_23a = addreaction(pfpc_model, 'ampar_CaPKC_act -> amparp + CaPKC_act');
ampar_reac_23a_k = addkineticlaw(ampar_reac_23a,'MassAction');
set(ampar_reac_23a_k, 'ParameterVariableNames', {'kcat_pkc_ampar'});

ampar_reac_23b = addreaction(pfpc_model, 'ampar_Ca2PKC_act -> amparp + Ca2PKC_act');
ampar_reac_23b_k = addkineticlaw(ampar_reac_23b,'MassAction');
set(ampar_reac_23b_k, 'ParameterVariableNames', {'kcat_pkc_ampar'});

ampar_reac_23c = addreaction(pfpc_model, 'ampar_CaPKC_AA_act -> amparp + CaPKC_AA_act');
ampar_reac_23c_k = addkineticlaw(ampar_reac_23c,'MassAction');
set(ampar_reac_23c_k, 'ParameterVariableNames', {'kcat_pkc_ampar'});

ampar_reac_23d = addreaction(pfpc_model, 'ampar_Ca2PKC_AA_act -> amparp + Ca2PKC_AA_act');
ampar_reac_23d_k = addkineticlaw(ampar_reac_23d,'MassAction');
set(ampar_reac_23d_k, 'ParameterVariableNames', {'kcat_pkc_ampar'});

ampar_reac_23e = addreaction(pfpc_model, 'ampar_DAG_CaPKC_act -> amparp + DAG_CaPKC_act');
ampar_reac_23e_k = addkineticlaw(ampar_reac_23e,'MassAction');
set(ampar_reac_23e_k, 'ParameterVariableNames', {'kcat_pkc_ampar'});

ampar_reac_23f = addreaction(pfpc_model, 'ampar_DAG_Ca2PKC_act -> amparp + DAG_Ca2PKC_act');
ampar_reac_23f_k = addkineticlaw(ampar_reac_23f,'MassAction');
set(ampar_reac_23f_k, 'ParameterVariableNames', {'kcat_pkc_ampar'});

ampar_reac_23g = addreaction(pfpc_model, 'ampar_DAG_CaPKC_AA_act -> amparp + DAG_CaPKC_AA_act');
ampar_reac_23g_k = addkineticlaw(ampar_reac_23g,'MassAction');
set(ampar_reac_23g_k, 'ParameterVariableNames', {'kcat_pkc_ampar'});

ampar_reac_23h = addreaction(pfpc_model, 'ampar_DAG_Ca2PKC_AA_act -> amparp + DAG_Ca2PKC_AA_act');
ampar_reac_23h_k = addkineticlaw(ampar_reac_23h,'MassAction');
set(ampar_reac_23h_k, 'ParameterVariableNames', {'kcat_pkc_ampar'});

ampar_reac_24a = addreaction(pfpc_model, 'ampar_x + CaPKC_act -> ampar_x_CaPKC_act');
ampar_reac_24a_k = addkineticlaw(ampar_reac_24a,'MassAction');
set(ampar_reac_24a_k, 'ParameterVariableNames', {'kon_ampar_pkc'});

ampar_reac_24b = addreaction(pfpc_model, 'ampar_x + Ca2PKC_act -> ampar_x_Ca2PKC_act');
ampar_reac_24b_k = addkineticlaw(ampar_reac_24b,'MassAction');
set(ampar_reac_24b_k, 'ParameterVariableNames', {'kon_ampar_pkc'});

ampar_reac_24c = addreaction(pfpc_model, 'ampar_x + CaPKC_AA_act -> ampar_x_CaPKC_AA_act');
ampar_reac_24c_k = addkineticlaw(ampar_reac_24c,'MassAction');
set(ampar_reac_24c_k, 'ParameterVariableNames', {'kon_ampar_pkc'});

ampar_reac_24d = addreaction(pfpc_model, 'ampar_x + Ca2PKC_AA_act -> ampar_x_Ca2PKC_AA_act');
ampar_reac_24d_k = addkineticlaw(ampar_reac_24d,'MassAction');
set(ampar_reac_24d_k, 'ParameterVariableNames', {'kon_ampar_pkc'});

ampar_reac_24e = addreaction(pfpc_model, 'ampar_x + DAG_CaPKC_act -> ampar_x_DAG_CaPKC_act');
ampar_reac_24e_k = addkineticlaw(ampar_reac_24e,'MassAction');
set(ampar_reac_24e_k, 'ParameterVariableNames', {'kon_ampar_pkc'});

ampar_reac_24f = addreaction(pfpc_model, 'ampar_x + DAG_Ca2PKC_act -> ampar_x_DAG_Ca2PKC_act');
ampar_reac_24f_k = addkineticlaw(ampar_reac_24f,'MassAction');
set(ampar_reac_24f_k, 'ParameterVariableNames', {'kon_ampar_pkc'});

ampar_reac_24g = addreaction(pfpc_model, 'ampar_x + DAG_CaPKC_AA_act -> ampar_x_DAG_CaPKC_AA_act');
ampar_reac_24g_k = addkineticlaw(ampar_reac_24g,'MassAction');
set(ampar_reac_24g_k, 'ParameterVariableNames', {'kon_ampar_pkc'});

ampar_reac_24h = addreaction(pfpc_model, 'ampar_x + DAG_Ca2PKC_AA_act -> ampar_x_DAG_Ca2PKC_AA_act');
ampar_reac_24h_k = addkineticlaw(ampar_reac_24h,'MassAction');
set(ampar_reac_24h_k, 'ParameterVariableNames', {'kon_ampar_pkc'});

ampar_reac_25a = addreaction(pfpc_model, 'ampar_x_CaPKC_act -> ampar_x + CaPKC_act');
ampar_reac_25a_k = addkineticlaw(ampar_reac_25a,'MassAction');
set(ampar_reac_25a_k, 'ParameterVariableNames', {'koff_ampar_pkc'});

ampar_reac_25b = addreaction(pfpc_model, 'ampar_x_Ca2PKC_act -> ampar_x + Ca2PKC_act');
ampar_reac_25b_k = addkineticlaw(ampar_reac_25b,'MassAction');
set(ampar_reac_25b_k, 'ParameterVariableNames', {'koff_ampar_pkc'});

ampar_reac_25c = addreaction(pfpc_model, 'ampar_x_CaPKC_AA_act -> ampar_x + CaPKC_AA_act');
ampar_reac_25c_k = addkineticlaw(ampar_reac_25c,'MassAction');
set(ampar_reac_25c_k, 'ParameterVariableNames', {'koff_ampar_pkc'});

ampar_reac_25d = addreaction(pfpc_model, 'ampar_x_Ca2PKC_AA_act -> ampar_x + Ca2PKC_AA_act');
ampar_reac_25d_k = addkineticlaw(ampar_reac_25d,'MassAction');
set(ampar_reac_25d_k, 'ParameterVariableNames', {'koff_ampar_pkc'});

ampar_reac_25e = addreaction(pfpc_model, 'ampar_x_DAG_CaPKC_act -> ampar_x + DAG_CaPKC_act');
ampar_reac_25e_k = addkineticlaw(ampar_reac_25e,'MassAction');
set(ampar_reac_25e_k, 'ParameterVariableNames', {'koff_ampar_pkc'});

ampar_reac_25f = addreaction(pfpc_model, 'ampar_x_DAG_Ca2PKC_act -> ampar_x + DAG_Ca2PKC_act');
ampar_reac_25f_k = addkineticlaw(ampar_reac_25f,'MassAction');
set(ampar_reac_25f_k, 'ParameterVariableNames', {'koff_ampar_pkc'});

ampar_reac_25g = addreaction(pfpc_model, 'ampar_x_DAG_CaPKC_AA_act -> ampar_x + DAG_CaPKC_AA_act');
ampar_reac_25g_k = addkineticlaw(ampar_reac_25g,'MassAction');
set(ampar_reac_25g_k, 'ParameterVariableNames', {'koff_ampar_pkc'});

ampar_reac_25h = addreaction(pfpc_model, 'ampar_x_DAG_Ca2PKC_AA_act -> ampar_x + DAG_Ca2PKC_AA_act');
ampar_reac_25h_k = addkineticlaw(ampar_reac_25h,'MassAction');
set(ampar_reac_25h_k, 'ParameterVariableNames', {'koff_ampar_pkc'});

ampar_reac_26a = addreaction(pfpc_model, 'ampar_x_CaPKC_act -> amparp_x + CaPKC_act');
ampar_reac_26a_k = addkineticlaw(ampar_reac_26a,'MassAction');
set(ampar_reac_26a_k, 'ParameterVariableNames', {'kcat_pkc_ampar'});

ampar_reac_26b = addreaction(pfpc_model, 'ampar_x_Ca2PKC_act -> amparp_x + Ca2PKC_act');
ampar_reac_26b_k = addkineticlaw(ampar_reac_26b,'MassAction');
set(ampar_reac_26b_k, 'ParameterVariableNames', {'kcat_pkc_ampar'});

ampar_reac_26c = addreaction(pfpc_model, 'ampar_x_CaPKC_AA_act -> amparp_x + CaPKC_AA_act');
ampar_reac_26c_k = addkineticlaw(ampar_reac_26c,'MassAction');
set(ampar_reac_26c_k, 'ParameterVariableNames', {'kcat_pkc_ampar'});

ampar_reac_26d = addreaction(pfpc_model, 'ampar_x_Ca2PKC_AA_act -> amparp_x + Ca2PKC_AA_act');
ampar_reac_26d_k = addkineticlaw(ampar_reac_26d,'MassAction');
set(ampar_reac_26d_k, 'ParameterVariableNames', {'kcat_pkc_ampar'});

ampar_reac_26e = addreaction(pfpc_model, 'ampar_x_DAG_CaPKC_act -> amparp_x + DAG_CaPKC_act');
ampar_reac_26e_k = addkineticlaw(ampar_reac_26e,'MassAction');
set(ampar_reac_26e_k, 'ParameterVariableNames', {'kcat_pkc_ampar'});

ampar_reac_26f = addreaction(pfpc_model, 'ampar_x_DAG_Ca2PKC_act -> amparp_x + DAG_Ca2PKC_act');
ampar_reac_26f_k = addkineticlaw(ampar_reac_26f,'MassAction');
set(ampar_reac_26f_k, 'ParameterVariableNames', {'kcat_pkc_ampar'});

ampar_reac_26g = addreaction(pfpc_model, 'ampar_x_DAG_CaPKC_AA_act -> amparp_x + DAG_CaPKC_AA_act');
ampar_reac_26g_k = addkineticlaw(ampar_reac_26g,'MassAction');
set(ampar_reac_26g_k, 'ParameterVariableNames', {'kcat_pkc_ampar'});

ampar_reac_26h = addreaction(pfpc_model, 'ampar_x_DAG_Ca2PKC_AA_act -> amparp_x + DAG_Ca2PKC_AA_act');
ampar_reac_26h_k = addkineticlaw(ampar_reac_26h,'MassAction');
set(ampar_reac_26h_k, 'ParameterVariableNames', {'kcat_pkc_ampar'});

ampar_reac_27a = addreaction(pfpc_model, 'ampar_pick + CaPKC_act -> ampar_pick_CaPKC_act');
ampar_reac_27a_k = addkineticlaw(ampar_reac_27a,'MassAction');
set(ampar_reac_27a_k, 'ParameterVariableNames', {'kon_ampar_pkc'});

ampar_reac_27b = addreaction(pfpc_model, 'ampar_pick + Ca2PKC_act -> ampar_pick_Ca2PKC_act');
ampar_reac_27b_k = addkineticlaw(ampar_reac_27b,'MassAction');
set(ampar_reac_27b_k, 'ParameterVariableNames', {'kon_ampar_pkc'});

ampar_reac_27c = addreaction(pfpc_model, 'ampar_pick + CaPKC_AA_act -> ampar_pick_CaPKC_AA_act');
ampar_reac_27c_k = addkineticlaw(ampar_reac_27c,'MassAction');
set(ampar_reac_27c_k, 'ParameterVariableNames', {'kon_ampar_pkc'});

ampar_reac_27d = addreaction(pfpc_model, 'ampar_pick + Ca2PKC_AA_act -> ampar_pick_Ca2PKC_AA_act');
ampar_reac_27d_k = addkineticlaw(ampar_reac_27d,'MassAction');
set(ampar_reac_27d_k, 'ParameterVariableNames', {'kon_ampar_pkc'});

ampar_reac_27e = addreaction(pfpc_model, 'ampar_pick + DAG_CaPKC_act -> ampar_pick_DAG_CaPKC_act');
ampar_reac_27e_k = addkineticlaw(ampar_reac_27e,'MassAction');
set(ampar_reac_27e_k, 'ParameterVariableNames', {'kon_ampar_pkc'});

ampar_reac_27f = addreaction(pfpc_model, 'ampar_pick + DAG_Ca2PKC_act -> ampar_pick_DAG_Ca2PKC_act');
ampar_reac_27f_k = addkineticlaw(ampar_reac_27f,'MassAction');
set(ampar_reac_27f_k, 'ParameterVariableNames', {'kon_ampar_pkc'});

ampar_reac_27g = addreaction(pfpc_model, 'ampar_pick + DAG_CaPKC_AA_act -> ampar_pick_DAG_CaPKC_AA_act');
ampar_reac_27g_k = addkineticlaw(ampar_reac_27g,'MassAction');
set(ampar_reac_27g_k, 'ParameterVariableNames', {'kon_ampar_pkc'});

ampar_reac_27h = addreaction(pfpc_model, 'ampar_pick + DAG_Ca2PKC_AA_act -> ampar_pick_DAG_Ca2PKC_AA_act');
ampar_reac_27h_k = addkineticlaw(ampar_reac_27h,'MassAction');
set(ampar_reac_27h_k, 'ParameterVariableNames', {'kon_ampar_pkc'});

ampar_reac_28a = addreaction(pfpc_model, 'ampar_pick_CaPKC_act -> ampar_pick + CaPKC_act');
ampar_reac_28a_k = addkineticlaw(ampar_reac_28a,'MassAction');
set(ampar_reac_28a_k, 'ParameterVariableNames', {'koff_ampar_pkc'});

ampar_reac_28b = addreaction(pfpc_model, 'ampar_pick_Ca2PKC_act -> ampar_pick + Ca2PKC_act');
ampar_reac_28b_k = addkineticlaw(ampar_reac_28b,'MassAction');
set(ampar_reac_28b_k, 'ParameterVariableNames', {'koff_ampar_pkc'});

ampar_reac_28c = addreaction(pfpc_model, 'ampar_pick_CaPKC_AA_act -> ampar_pick + CaPKC_AA_act');
ampar_reac_28c_k = addkineticlaw(ampar_reac_28c,'MassAction');
set(ampar_reac_28c_k, 'ParameterVariableNames', {'koff_ampar_pkc'});

ampar_reac_28d = addreaction(pfpc_model, 'ampar_pick_Ca2PKC_AA_act -> ampar_pick + Ca2PKC_AA_act');
ampar_reac_28d_k = addkineticlaw(ampar_reac_28d,'MassAction');
set(ampar_reac_28d_k, 'ParameterVariableNames', {'koff_ampar_pkc'});

ampar_reac_28e = addreaction(pfpc_model, 'ampar_pick_DAG_CaPKC_act -> ampar_pick + DAG_CaPKC_act');
ampar_reac_28e_k = addkineticlaw(ampar_reac_28e,'MassAction');
set(ampar_reac_28e_k, 'ParameterVariableNames', {'koff_ampar_pkc'});

ampar_reac_28f = addreaction(pfpc_model, 'ampar_pick_DAG_Ca2PKC_act -> ampar_pick + DAG_Ca2PKC_act');
ampar_reac_28f_k = addkineticlaw(ampar_reac_28f,'MassAction');
set(ampar_reac_28f_k, 'ParameterVariableNames', {'koff_ampar_pkc'});

ampar_reac_28g = addreaction(pfpc_model, 'ampar_pick_DAG_CaPKC_AA_act -> ampar_pick + DAG_CaPKC_AA_act');
ampar_reac_28g_k = addkineticlaw(ampar_reac_28g,'MassAction');
set(ampar_reac_28g_k, 'ParameterVariableNames', {'koff_ampar_pkc'});

ampar_reac_28h = addreaction(pfpc_model, 'ampar_pick_DAG_Ca2PKC_AA_act -> ampar_pick + DAG_Ca2PKC_AA_act');
ampar_reac_28h_k = addkineticlaw(ampar_reac_28h,'MassAction');
set(ampar_reac_28h_k, 'ParameterVariableNames', {'koff_ampar_pkc'});

ampar_reac_29a = addreaction(pfpc_model, 'ampar_pick_CaPKC_act -> amparp_pick + CaPKC_act');
ampar_reac_29a_k = addkineticlaw(ampar_reac_29a,'MassAction');
set(ampar_reac_29a_k, 'ParameterVariableNames', {'kcat_pkc_ampar'});

ampar_reac_29b = addreaction(pfpc_model, 'ampar_pick_Ca2PKC_act -> amparp_pick + Ca2PKC_act');
ampar_reac_29b_k = addkineticlaw(ampar_reac_29b,'MassAction');
set(ampar_reac_29b_k, 'ParameterVariableNames', {'kcat_pkc_ampar'});

ampar_reac_29c = addreaction(pfpc_model, 'ampar_pick_CaPKC_AA_act -> amparp_pick + CaPKC_AA_act');
ampar_reac_29c_k = addkineticlaw(ampar_reac_29c,'MassAction');
set(ampar_reac_29c_k, 'ParameterVariableNames', {'kcat_pkc_ampar'});

ampar_reac_29d = addreaction(pfpc_model, 'ampar_pick_Ca2PKC_AA_act -> amparp_pick + Ca2PKC_AA_act');
ampar_reac_29d_k = addkineticlaw(ampar_reac_29d,'MassAction');
set(ampar_reac_29d_k, 'ParameterVariableNames', {'kcat_pkc_ampar'});

ampar_reac_29e = addreaction(pfpc_model, 'ampar_pick_DAG_CaPKC_act -> amparp_pick + DAG_CaPKC_act');
ampar_reac_29e_k = addkineticlaw(ampar_reac_29e,'MassAction');
set(ampar_reac_29e_k, 'ParameterVariableNames', {'kcat_pkc_ampar'});

ampar_reac_29f = addreaction(pfpc_model, 'ampar_pick_DAG_Ca2PKC_act -> amparp_pick + DAG_Ca2PKC_act');
ampar_reac_29f_k = addkineticlaw(ampar_reac_29f,'MassAction');
set(ampar_reac_29f_k, 'ParameterVariableNames', {'kcat_pkc_ampar'});

ampar_reac_29g = addreaction(pfpc_model, 'ampar_pick_DAG_CaPKC_AA_act -> amparp_pick + DAG_CaPKC_AA_act');
ampar_reac_29g_k = addkineticlaw(ampar_reac_29g,'MassAction');
set(ampar_reac_29g_k, 'ParameterVariableNames', {'kcat_pkc_ampar'});

ampar_reac_29h = addreaction(pfpc_model, 'ampar_pick_DAG_Ca2PKC_AA_act -> amparp_pick + DAG_Ca2PKC_AA_act');
ampar_reac_29h_k = addkineticlaw(ampar_reac_29h,'MassAction');
set(ampar_reac_29h_k, 'ParameterVariableNames', {'kcat_pkc_ampar'});

ampar_reac_30a = addreaction(pfpc_model, 'ampar_s + CaPKC_act -> ampar_s_CaPKC_act');
ampar_reac_30a_k = addkineticlaw(ampar_reac_30a,'MassAction');
set(ampar_reac_30a_k, 'ParameterVariableNames', {'kon_ampar_pkc'});

ampar_reac_30b = addreaction(pfpc_model, 'ampar_s + Ca2PKC_act -> ampar_s_Ca2PKC_act');
ampar_reac_30b_k = addkineticlaw(ampar_reac_30b,'MassAction');
set(ampar_reac_30b_k, 'ParameterVariableNames', {'kon_ampar_pkc'});

ampar_reac_30c = addreaction(pfpc_model, 'ampar_s + CaPKC_AA_act -> ampar_s_CaPKC_AA_act');
ampar_reac_30c_k = addkineticlaw(ampar_reac_30c,'MassAction');
set(ampar_reac_30c_k, 'ParameterVariableNames', {'kon_ampar_pkc'});

ampar_reac_30d = addreaction(pfpc_model, 'ampar_s + Ca2PKC_AA_act -> ampar_s_Ca2PKC_AA_act');
ampar_reac_30d_k = addkineticlaw(ampar_reac_30d,'MassAction');
set(ampar_reac_30d_k, 'ParameterVariableNames', {'kon_ampar_pkc'});

ampar_reac_30e = addreaction(pfpc_model, 'ampar_s + DAG_CaPKC_act -> ampar_s_DAG_CaPKC_act');
ampar_reac_30e_k = addkineticlaw(ampar_reac_30e,'MassAction');
set(ampar_reac_30e_k, 'ParameterVariableNames', {'kon_ampar_pkc'});

ampar_reac_30f = addreaction(pfpc_model, 'ampar_s + DAG_Ca2PKC_act -> ampar_s_DAG_Ca2PKC_act');
ampar_reac_30f_k = addkineticlaw(ampar_reac_30f,'MassAction');
set(ampar_reac_30f_k, 'ParameterVariableNames', {'kon_ampar_pkc'});

ampar_reac_30g = addreaction(pfpc_model, 'ampar_s + DAG_CaPKC_AA_act -> ampar_s_DAG_CaPKC_AA_act');
ampar_reac_30g_k = addkineticlaw(ampar_reac_30g,'MassAction');
set(ampar_reac_30g_k, 'ParameterVariableNames', {'kon_ampar_pkc'});

ampar_reac_30h = addreaction(pfpc_model, 'ampar_s + DAG_Ca2PKC_AA_act -> ampar_s_DAG_Ca2PKC_AA_act');
ampar_reac_30h_k = addkineticlaw(ampar_reac_30h,'MassAction');
set(ampar_reac_30h_k, 'ParameterVariableNames', {'kon_ampar_pkc'});

ampar_reac_31a = addreaction(pfpc_model, 'ampar_s_CaPKC_act -> ampar_s + CaPKC_act');
ampar_reac_31a_k = addkineticlaw(ampar_reac_31a,'MassAction');
set(ampar_reac_31a_k, 'ParameterVariableNames', {'koff_ampar_pkc'});

ampar_reac_31b = addreaction(pfpc_model, 'ampar_s_Ca2PKC_act -> ampar_s + Ca2PKC_act');
ampar_reac_31b_k = addkineticlaw(ampar_reac_31b,'MassAction');
set(ampar_reac_31b_k, 'ParameterVariableNames', {'koff_ampar_pkc'});

ampar_reac_31c = addreaction(pfpc_model, 'ampar_s_CaPKC_AA_act -> ampar_s + CaPKC_AA_act');
ampar_reac_31c_k = addkineticlaw(ampar_reac_31c,'MassAction');
set(ampar_reac_31c_k, 'ParameterVariableNames', {'koff_ampar_pkc'});

ampar_reac_31d = addreaction(pfpc_model, 'ampar_s_Ca2PKC_AA_act -> ampar_s + Ca2PKC_AA_act');
ampar_reac_31d_k = addkineticlaw(ampar_reac_31d,'MassAction');
set(ampar_reac_31d_k, 'ParameterVariableNames', {'koff_ampar_pkc'});

ampar_reac_31e = addreaction(pfpc_model, 'ampar_s_DAG_CaPKC_act -> ampar_s + DAG_CaPKC_act');
ampar_reac_31e_k = addkineticlaw(ampar_reac_31e,'MassAction');
set(ampar_reac_31e_k, 'ParameterVariableNames', {'koff_ampar_pkc'});

ampar_reac_31f = addreaction(pfpc_model, 'ampar_s_DAG_Ca2PKC_act -> ampar_s + DAG_Ca2PKC_act');
ampar_reac_31f_k = addkineticlaw(ampar_reac_31f,'MassAction');
set(ampar_reac_31f_k, 'ParameterVariableNames', {'koff_ampar_pkc'});

ampar_reac_31g = addreaction(pfpc_model, 'ampar_s_DAG_CaPKC_AA_act -> ampar_s + DAG_CaPKC_AA_act');
ampar_reac_31g_k = addkineticlaw(ampar_reac_31g,'MassAction');
set(ampar_reac_31g_k, 'ParameterVariableNames', {'koff_ampar_pkc'});

ampar_reac_31h = addreaction(pfpc_model, 'ampar_s_DAG_Ca2PKC_AA_act -> ampar_s + DAG_Ca2PKC_AA_act');
ampar_reac_31h_k = addkineticlaw(ampar_reac_31h,'MassAction');
set(ampar_reac_31h_k, 'ParameterVariableNames', {'koff_ampar_pkc'});

ampar_reac_32a = addreaction(pfpc_model, 'ampar_s_pCaPKC_actkc -> amparp_s + CaPKC_act');
ampar_reac_32a_k = addkineticlaw(ampar_reac_32a,'MassAction');
set(ampar_reac_32a_k, 'ParameterVariableNames', {'kcat_pkc_ampar'});

ampar_reac_32b = addreaction(pfpc_model, 'ampar_s_Ca2PKC_act -> amparp_s + Ca2PKC_act');
ampar_reac_32b_k = addkineticlaw(ampar_reac_32b,'MassAction');
set(ampar_reac_32b_k, 'ParameterVariableNames', {'kcat_pkc_ampar'});

ampar_reac_32c = addreaction(pfpc_model, 'ampar_s_CaPKC_AA_act -> amparp_s + CaPKC_AA_act');
ampar_reac_32c_k = addkineticlaw(ampar_reac_32c,'MassAction');
set(ampar_reac_32c_k, 'ParameterVariableNames', {'kcat_pkc_ampar'});

ampar_reac_32d = addreaction(pfpc_model, 'ampar_s_Ca2PKC_AA_act -> amparp_s + Ca2PKC_AA_act');
ampar_reac_32d_k = addkineticlaw(ampar_reac_32d,'MassAction');
set(ampar_reac_32d_k, 'ParameterVariableNames', {'kcat_pkc_ampar'});

ampar_reac_32e = addreaction(pfpc_model, 'ampar_s_DAG_CaPKC_act -> amparp_s + DAG_CaPKC_act');
ampar_reac_32e_k = addkineticlaw(ampar_reac_32e,'MassAction');
set(ampar_reac_32e_k, 'ParameterVariableNames', {'kcat_pkc_ampar'});

ampar_reac_32f = addreaction(pfpc_model, 'ampar_s_DAG_Ca2PKC_act -> amparp_s + DAG_Ca2PKC_act');
ampar_reac_32f_k = addkineticlaw(ampar_reac_32f,'MassAction');
set(ampar_reac_32f_k, 'ParameterVariableNames', {'kcat_pkc_ampar'});

ampar_reac_32g = addreaction(pfpc_model, 'ampar_s_DAG_CaPKC_AA_act -> amparp_s + DAG_CaPKC_AA_act');
ampar_reac_32g_k = addkineticlaw(ampar_reac_32g,'MassAction');
set(ampar_reac_32g_k, 'ParameterVariableNames', {'kcat_pkc_ampar'});

ampar_reac_32h = addreaction(pfpc_model, 'ampar_s_DAG_Ca2PKC_AA_act -> amparp_s + DAG_Ca2PKC_AA_act');
ampar_reac_32h_k = addkineticlaw(ampar_reac_32h,'MassAction');
set(ampar_reac_32h_k, 'ParameterVariableNames', {'kcat_pkc_ampar'});

ampar_reac_33a = addreaction(pfpc_model, 'ampar_s_pick + CaPKC_act -> ampar_s_pick_CaPKC_act');
ampar_reac_33a_k = addkineticlaw(ampar_reac_33a,'MassAction');
set(ampar_reac_33a_k, 'ParameterVariableNames', {'kon_ampar_pkc'});

ampar_reac_33b = addreaction(pfpc_model, 'ampar_s_pick + Ca2PKC_act -> ampar_s_pick_Ca2PKC_act');
ampar_reac_33b_k = addkineticlaw(ampar_reac_33b,'MassAction');
set(ampar_reac_33b_k, 'ParameterVariableNames', {'kon_ampar_pkc'});

ampar_reac_33c = addreaction(pfpc_model, 'ampar_s_pick + CaPKC_AA_act -> ampar_s_pick_CaPKC_AA_act');
ampar_reac_33c_k = addkineticlaw(ampar_reac_33c,'MassAction');
set(ampar_reac_33c_k, 'ParameterVariableNames', {'kon_ampar_pkc'});

ampar_reac_33d = addreaction(pfpc_model, 'ampar_s_pick + Ca2PKC_AA_act -> ampar_s_pick_Ca2PKC_AA_act');
ampar_reac_33d_k = addkineticlaw(ampar_reac_33d,'MassAction');
set(ampar_reac_33d_k, 'ParameterVariableNames', {'kon_ampar_pkc'});

ampar_reac_33e = addreaction(pfpc_model, 'ampar_s_pick + DAG_CaPKC_act -> ampar_s_pick_DAG_CaPKC_act');
ampar_reac_33e_k = addkineticlaw(ampar_reac_33e,'MassAction');
set(ampar_reac_33e_k, 'ParameterVariableNames', {'kon_ampar_pkc'});

ampar_reac_33f = addreaction(pfpc_model, 'ampar_s_pick + DAG_Ca2PKC_act -> ampar_s_pick_DAG_Ca2PKC_act');
ampar_reac_33f_k = addkineticlaw(ampar_reac_33f,'MassAction');
set(ampar_reac_33f_k, 'ParameterVariableNames', {'kon_ampar_pkc'});

ampar_reac_33g = addreaction(pfpc_model, 'ampar_s_pick + DAG_CaPKC_AA_act -> ampar_s_pick_DAG_CaPKC_AA_act');
ampar_reac_33g_k = addkineticlaw(ampar_reac_33g,'MassAction');
set(ampar_reac_33g_k, 'ParameterVariableNames', {'kon_ampar_pkc'});

ampar_reac_33h = addreaction(pfpc_model, 'ampar_s_pick + DAG_Ca2PKC_AA_act -> ampar_s_pick_DAG_Ca2PKC_AA_act');
ampar_reac_33h_k = addkineticlaw(ampar_reac_33h,'MassAction');
set(ampar_reac_33h_k, 'ParameterVariableNames', {'kon_ampar_pkc'});

ampar_reac_34a = addreaction(pfpc_model, 'ampar_s_pick_CaPKC_act -> ampar_s_pick + CaPKC_act');
ampar_reac_34a_k = addkineticlaw(ampar_reac_34a,'MassAction');
set(ampar_reac_34a_k, 'ParameterVariableNames', {'koff_ampar_pkc'});

ampar_reac_34b = addreaction(pfpc_model, 'ampar_s_pick_Ca2PKC_act -> ampar_s_pick + Ca2PKC_act');
ampar_reac_34b_k = addkineticlaw(ampar_reac_34b,'MassAction');
set(ampar_reac_34b_k, 'ParameterVariableNames', {'koff_ampar_pkc'});

ampar_reac_34c = addreaction(pfpc_model, 'ampar_s_pick_CaPKC_AA_act -> ampar_s_pick + CaPKC_AA_act');
ampar_reac_34c_k = addkineticlaw(ampar_reac_34c,'MassAction');
set(ampar_reac_34c_k, 'ParameterVariableNames', {'koff_ampar_pkc'});

ampar_reac_34d = addreaction(pfpc_model, 'ampar_s_pick_Ca2PKC_AA_act -> ampar_s_pick + Ca2PKC_AA_act');
ampar_reac_34d_k = addkineticlaw(ampar_reac_34d,'MassAction');
set(ampar_reac_34d_k, 'ParameterVariableNames', {'koff_ampar_pkc'});

ampar_reac_34e = addreaction(pfpc_model, 'ampar_s_pick_DAG_CaPKC_act -> ampar_s_pick + DAG_CaPKC_act');
ampar_reac_34e_k = addkineticlaw(ampar_reac_34e,'MassAction');
set(ampar_reac_34e_k, 'ParameterVariableNames', {'koff_ampar_pkc'});

ampar_reac_34f = addreaction(pfpc_model, 'ampar_s_pick_DAG_Ca2PKC_act -> ampar_s_pick + DAG_Ca2PKC_act');
ampar_reac_34f_k = addkineticlaw(ampar_reac_34f,'MassAction');
set(ampar_reac_34f_k, 'ParameterVariableNames', {'koff_ampar_pkc'});

ampar_reac_34g = addreaction(pfpc_model, 'ampar_s_pick_DAG_CaPKC_AA_act -> ampar_s_pick + DAG_CaPKC_AA_act');
ampar_reac_34g_k = addkineticlaw(ampar_reac_34g,'MassAction');
set(ampar_reac_34g_k, 'ParameterVariableNames', {'koff_ampar_pkc'});

ampar_reac_34h = addreaction(pfpc_model, 'ampar_s_pick_DAG_Ca2PKC_AA_act -> ampar_s_pick + DAG_Ca2PKC_AA_act');
ampar_reac_34h_k = addkineticlaw(ampar_reac_34h,'MassAction');
set(ampar_reac_34h_k, 'ParameterVariableNames', {'koff_ampar_pkc'});

ampar_reac_35a = addreaction(pfpc_model, 'ampar_s_pick_CaPKC_act -> amparp_s_pick + CaPKC_act');
ampar_reac_35a_k = addkineticlaw(ampar_reac_35a,'MassAction');
set(ampar_reac_35a_k, 'ParameterVariableNames', {'kcat_pkc_ampar'});

ampar_reac_35b = addreaction(pfpc_model, 'ampar_s_pick_Ca2PKC_act -> amparp_s_pick + Ca2PKC_act');
ampar_reac_35b_k = addkineticlaw(ampar_reac_35b,'MassAction');
set(ampar_reac_35b_k, 'ParameterVariableNames', {'kcat_pkc_ampar'});

ampar_reac_35c = addreaction(pfpc_model, 'ampar_s_pick_CaPKC_AA_act -> amparp_s_pick + CaPKC_AA_act');
ampar_reac_35c_k = addkineticlaw(ampar_reac_35c,'MassAction');
set(ampar_reac_35c_k, 'ParameterVariableNames', {'kcat_pkc_ampar'});

ampar_reac_35d = addreaction(pfpc_model, 'ampar_s_pick_Ca2PKC_AA_act -> amparp_s_pick + Ca2PKC_AA_act');
ampar_reac_35d_k = addkineticlaw(ampar_reac_35d,'MassAction');
set(ampar_reac_35d_k, 'ParameterVariableNames', {'kcat_pkc_ampar'});

ampar_reac_35e = addreaction(pfpc_model, 'ampar_s_pick_DAG_CaPKC_act -> amparp_s_pick + DAG_CaPKC_act');
ampar_reac_35e_k = addkineticlaw(ampar_reac_35e,'MassAction');
set(ampar_reac_35e_k, 'ParameterVariableNames', {'kcat_pkc_ampar'});

ampar_reac_35f = addreaction(pfpc_model, 'ampar_s_pick_DAG_Ca2PKC_act -> amparp_s_pick + DAG_Ca2PKC_act');
ampar_reac_35f_k = addkineticlaw(ampar_reac_35f,'MassAction');
set(ampar_reac_35f_k, 'ParameterVariableNames', {'kcat_pkc_ampar'});

ampar_reac_35g = addreaction(pfpc_model, 'ampar_s_pick_DAG_CaPKC_AA_act -> amparp_s_pick + DAG_CaPKC_AA_act');
ampar_reac_35g_k = addkineticlaw(ampar_reac_35g,'MassAction');
set(ampar_reac_35g_k, 'ParameterVariableNames', {'kcat_pkc_ampar'});

ampar_reac_35h = addreaction(pfpc_model, 'ampar_s_pick_DAG_Ca2PKC_AA_act -> amparp_s_pick + DAG_Ca2PKC_AA_act');
ampar_reac_35h_k = addkineticlaw(ampar_reac_35h,'MassAction');
set(ampar_reac_35h_k, 'ParameterVariableNames', {'kcat_pkc_ampar'});

ampar_reac_36a = addreaction(pfpc_model, 'ampar_ez + CaPKC_act -> ampar_ez_CaPKC_act');
ampar_reac_36a_k = addkineticlaw(ampar_reac_36a,'MassAction');
set(ampar_reac_36a_k, 'ParameterVariableNames', {'kon_ampar_pkc'});

ampar_reac_36b = addreaction(pfpc_model, 'ampar_ez + Ca2PKC_act -> ampar_ez_Ca2PKC_act');
ampar_reac_36b_k = addkineticlaw(ampar_reac_36b,'MassAction');
set(ampar_reac_36b_k, 'ParameterVariableNames', {'kon_ampar_pkc'});

ampar_reac_36c = addreaction(pfpc_model, 'ampar_ez + CaPKC_AA_act -> ampar_ez_CaPKC_AA_act');
ampar_reac_36c_k = addkineticlaw(ampar_reac_36c,'MassAction');
set(ampar_reac_36c_k, 'ParameterVariableNames', {'kon_ampar_pkc'});

ampar_reac_36d = addreaction(pfpc_model, 'ampar_ez + Ca2PKC_AA_act -> ampar_ez_Ca2PKC_AA_act');
ampar_reac_36d_k = addkineticlaw(ampar_reac_36d,'MassAction');
set(ampar_reac_36d_k, 'ParameterVariableNames', {'kon_ampar_pkc'});

ampar_reac_36e = addreaction(pfpc_model, 'ampar_ez + DAG_CaPKC_act -> ampar_ez_DAG_CaPKC_act');
ampar_reac_36e_k = addkineticlaw(ampar_reac_36e,'MassAction');
set(ampar_reac_36e_k, 'ParameterVariableNames', {'kon_ampar_pkc'});

ampar_reac_36f = addreaction(pfpc_model, 'ampar_ez + DAG_Ca2PKC_act -> ampar_ez_DAG_Ca2PKC_act');
ampar_reac_36f_k = addkineticlaw(ampar_reac_36f,'MassAction');
set(ampar_reac_36f_k, 'ParameterVariableNames', {'kon_ampar_pkc'});

ampar_reac_36g = addreaction(pfpc_model, 'ampar_ez + DAG_CaPKC_AA_act -> ampar_ez_DAG_CaPKC_AA_act');
ampar_reac_36g_k = addkineticlaw(ampar_reac_36g,'MassAction');
set(ampar_reac_36g_k, 'ParameterVariableNames', {'kon_ampar_pkc'});

ampar_reac_36h = addreaction(pfpc_model, 'ampar_ez + DAG_Ca2PKC_AA_act -> ampar_ez_DAG_Ca2PKC_AA_act');
ampar_reac_36h_k = addkineticlaw(ampar_reac_36h,'MassAction');
set(ampar_reac_36h_k, 'ParameterVariableNames', {'kon_ampar_pkc'});

ampar_reac_37a = addreaction(pfpc_model, 'ampar_ez_CaPKC_act -> ampar_ez + CaPKC_act');
ampar_reac_37a_k = addkineticlaw(ampar_reac_37a,'MassAction');
set(ampar_reac_37a_k, 'ParameterVariableNames', {'koff_ampar_pkc'});

ampar_reac_37b = addreaction(pfpc_model, 'ampar_ez_Ca2PKC_act -> ampar_ez + Ca2PKC_act');
ampar_reac_37b_k = addkineticlaw(ampar_reac_37b,'MassAction');
set(ampar_reac_37b_k, 'ParameterVariableNames', {'koff_ampar_pkc'});

ampar_reac_37c = addreaction(pfpc_model, 'ampar_ez_CaPKC_AA_act -> ampar_ez + CaPKC_AA_act');
ampar_reac_37c_k = addkineticlaw(ampar_reac_37c,'MassAction');
set(ampar_reac_37c_k, 'ParameterVariableNames', {'koff_ampar_pkc'});

ampar_reac_37d = addreaction(pfpc_model, 'ampar_ez_Ca2PKC_AA_act -> ampar_ez + Ca2PKC_AA_act');
ampar_reac_37d_k = addkineticlaw(ampar_reac_37d,'MassAction');
set(ampar_reac_37d_k, 'ParameterVariableNames', {'koff_ampar_pkc'});

ampar_reac_37e = addreaction(pfpc_model, 'ampar_ez_DAG_CaPKC_act -> ampar_ez + DAG_CaPKC_act');
ampar_reac_37e_k = addkineticlaw(ampar_reac_37e,'MassAction');
set(ampar_reac_37e_k, 'ParameterVariableNames', {'koff_ampar_pkc'});

ampar_reac_37f = addreaction(pfpc_model, 'ampar_ez_DAG_Ca2PKC_act -> ampar_ez + DAG_Ca2PKC_act');
ampar_reac_37f_k = addkineticlaw(ampar_reac_37f,'MassAction');
set(ampar_reac_37f_k, 'ParameterVariableNames', {'koff_ampar_pkc'});

ampar_reac_37g = addreaction(pfpc_model, 'ampar_ez_DAG_CaPKC_AA_act -> ampar_ez + DAG_CaPKC_AA_act');
ampar_reac_37g_k = addkineticlaw(ampar_reac_37g,'MassAction');
set(ampar_reac_37g_k, 'ParameterVariableNames', {'koff_ampar_pkc'});

ampar_reac_37h = addreaction(pfpc_model, 'ampar_ez_DAG_Ca2PKC_AA_act -> ampar_ez + DAG_Ca2PKC_AA_act');
ampar_reac_37h_k = addkineticlaw(ampar_reac_37h,'MassAction');
set(ampar_reac_37h_k, 'ParameterVariableNames', {'koff_ampar_pkc'});

ampar_reac_38a = addreaction(pfpc_model, 'ampar_ez_CaPKC_act -> amparp_ez + CaPKC_act');
ampar_reac_38a_k = addkineticlaw(ampar_reac_38a,'MassAction');
set(ampar_reac_38a_k, 'ParameterVariableNames', {'kcat_pkc_ampar'});

ampar_reac_38b = addreaction(pfpc_model, 'ampar_ez_Ca2PKC_act -> amparp_ez + Ca2PKC_act');
ampar_reac_38b_k = addkineticlaw(ampar_reac_38b,'MassAction');
set(ampar_reac_38b_k, 'ParameterVariableNames', {'kcat_pkc_ampar'});

ampar_reac_38c = addreaction(pfpc_model, 'ampar_ez_CaPKC_AA_act -> amparp_ez + CaPKC_AA_act');
ampar_reac_38c_k = addkineticlaw(ampar_reac_38c,'MassAction');
set(ampar_reac_38c_k, 'ParameterVariableNames', {'kcat_pkc_ampar'});

ampar_reac_38d = addreaction(pfpc_model, 'ampar_ez_Ca2PKC_AA_act -> amparp_ez + Ca2PKC_AA_act');
ampar_reac_38d_k = addkineticlaw(ampar_reac_38d,'MassAction');
set(ampar_reac_38d_k, 'ParameterVariableNames', {'kcat_pkc_ampar'});

ampar_reac_38e = addreaction(pfpc_model, 'ampar_ez_DAG_CaPKC_act -> amparp_ez + DAG_CaPKC_act');
ampar_reac_38e_k = addkineticlaw(ampar_reac_38e,'MassAction');
set(ampar_reac_38e_k, 'ParameterVariableNames', {'kcat_pkc_ampar'});

ampar_reac_38f = addreaction(pfpc_model, 'ampar_ez_DAG_Ca2PKC_act -> amparp_ez + DAG_Ca2PKC_act');
ampar_reac_38f_k = addkineticlaw(ampar_reac_38f,'MassAction');
set(ampar_reac_38f_k, 'ParameterVariableNames', {'kcat_pkc_ampar'});

ampar_reac_38g = addreaction(pfpc_model, 'ampar_ez_DAG_CaPKC_AA_act -> amparp_ez + DAG_CaPKC_AA_act');
ampar_reac_38g_k = addkineticlaw(ampar_reac_38g,'MassAction');
set(ampar_reac_38g_k, 'ParameterVariableNames', {'kcat_pkc_ampar'});

ampar_reac_38h = addreaction(pfpc_model, 'ampar_ez_DAG_Ca2PKC_AA_act -> amparp_ez + DAG_Ca2PKC_AA_act');
ampar_reac_38h_k = addkineticlaw(ampar_reac_38h,'MassAction');
set(ampar_reac_38h_k, 'ParameterVariableNames', {'kcat_pkc_ampar'});

ampar_reac_39a = addreaction(pfpc_model, 'ampar_ez_pick + CaPKC_act -> ampar_ez_pick_CaPKC_act');
ampar_reac_39a_k = addkineticlaw(ampar_reac_39a,'MassAction');
set(ampar_reac_39a_k, 'ParameterVariableNames', {'kon_ampar_pkc'});

ampar_reac_39b = addreaction(pfpc_model, 'ampar_ez_pick + Ca2PKC_act -> ampar_ez_pick_Ca2PKC_act');
ampar_reac_39b_k = addkineticlaw(ampar_reac_39b,'MassAction');
set(ampar_reac_39b_k, 'ParameterVariableNames', {'kon_ampar_pkc'});

ampar_reac_39c = addreaction(pfpc_model, 'ampar_ez_pick + CaPKC_AA_act -> ampar_ez_pick_CaPKC_AA_act');
ampar_reac_39c_k = addkineticlaw(ampar_reac_39c,'MassAction');
set(ampar_reac_39c_k, 'ParameterVariableNames', {'kon_ampar_pkc'});

ampar_reac_39d = addreaction(pfpc_model, 'ampar_ez_pick + Ca2PKC_AA_act -> ampar_ez_pick_Ca2PKC_AA_act');
ampar_reac_39d_k = addkineticlaw(ampar_reac_39d,'MassAction');
set(ampar_reac_39d_k, 'ParameterVariableNames', {'kon_ampar_pkc'});

ampar_reac_39e = addreaction(pfpc_model, 'ampar_ez_pick + DAG_CaPKC_act -> ampar_ez_pick_DAG_CaPKC_act');
ampar_reac_39e_k = addkineticlaw(ampar_reac_39e,'MassAction');
set(ampar_reac_39e_k, 'ParameterVariableNames', {'kon_ampar_pkc'});

ampar_reac_39f = addreaction(pfpc_model, 'ampar_ez_pick + DAG_Ca2PKC_act -> ampar_ez_pick_DAG_Ca2PKC_act');
ampar_reac_39f_k = addkineticlaw(ampar_reac_39f,'MassAction');
set(ampar_reac_39f_k, 'ParameterVariableNames', {'kon_ampar_pkc'});

ampar_reac_39g = addreaction(pfpc_model, 'ampar_ez_pick + DAG_CaPKC_AA_act -> ampar_ez_pick_DAG_CaPKC_AA_act');
ampar_reac_39g_k = addkineticlaw(ampar_reac_39g,'MassAction');
set(ampar_reac_39g_k, 'ParameterVariableNames', {'kon_ampar_pkc'});

ampar_reac_39h = addreaction(pfpc_model, 'ampar_ez_pick + DAG_Ca2PKC_AA_act -> ampar_ez_pick_DAG_Ca2PKC_AA_act');
ampar_reac_39h_k = addkineticlaw(ampar_reac_39h,'MassAction');
set(ampar_reac_39h_k, 'ParameterVariableNames', {'kon_ampar_pkc'});

ampar_reac_40a = addreaction(pfpc_model, 'ampar_ez_pick_CaPKC_act -> ampar_ez_pick + CaPKC_act');
ampar_reac_40a_k = addkineticlaw(ampar_reac_40a,'MassAction');
set(ampar_reac_40a_k, 'ParameterVariableNames', {'koff_ampar_pkc'});

ampar_reac_40b = addreaction(pfpc_model, 'ampar_ez_pick_Ca2PKC_act -> ampar_ez_pick + Ca2PKC_act');
ampar_reac_40b_k = addkineticlaw(ampar_reac_40b,'MassAction');
set(ampar_reac_40b_k, 'ParameterVariableNames', {'koff_ampar_pkc'});

ampar_reac_40c = addreaction(pfpc_model, 'ampar_ez_pick_CaPKC_AA_act -> ampar_ez_pick + CaPKC_AA_act');
ampar_reac_40c_k = addkineticlaw(ampar_reac_40c,'MassAction');
set(ampar_reac_40c_k, 'ParameterVariableNames', {'koff_ampar_pkc'});

ampar_reac_40d = addreaction(pfpc_model, 'ampar_ez_pick_Ca2PKC_AA_act -> ampar_ez_pick + Ca2PKC_AA_act');
ampar_reac_40d_k = addkineticlaw(ampar_reac_40d,'MassAction');
set(ampar_reac_40d_k, 'ParameterVariableNames', {'koff_ampar_pkc'});

ampar_reac_40e = addreaction(pfpc_model, 'ampar_ez_pick_DAG_CaPKC_act -> ampar_ez_pick + DAG_CaPKC_act');
ampar_reac_40e_k = addkineticlaw(ampar_reac_40e,'MassAction');
set(ampar_reac_40e_k, 'ParameterVariableNames', {'koff_ampar_pkc'});

ampar_reac_40f = addreaction(pfpc_model, 'ampar_ez_pick_DAG_Ca2PKC_act -> ampar_ez_pick + DAG_Ca2PKC_act');
ampar_reac_40f_k = addkineticlaw(ampar_reac_40f,'MassAction');
set(ampar_reac_40f_k, 'ParameterVariableNames', {'koff_ampar_pkc'});

ampar_reac_40g = addreaction(pfpc_model, 'ampar_ez_pick_DAG_CaPKC_AA_act -> ampar_ez_pick + DAG_CaPKC_AA_act');
ampar_reac_40g_k = addkineticlaw(ampar_reac_40g,'MassAction');
set(ampar_reac_40g_k, 'ParameterVariableNames', {'koff_ampar_pkc'});

ampar_reac_40h = addreaction(pfpc_model, 'ampar_ez_pick_DAG_Ca2PKC_AA_act -> ampar_ez_pick + DAG_Ca2PKC_AA_act');
ampar_reac_40h_k = addkineticlaw(ampar_reac_40h,'MassAction');
set(ampar_reac_40h_k, 'ParameterVariableNames', {'koff_ampar_pkc'});

ampar_reac_41a = addreaction(pfpc_model, 'ampar_ez_pick_CaPKC_act -> amparp_ez_pick + CaPKC_act');
ampar_reac_41a_k = addkineticlaw(ampar_reac_41a,'MassAction');
set(ampar_reac_41a_k, 'ParameterVariableNames', {'kcat_pkc_ampar'});

ampar_reac_41b = addreaction(pfpc_model, 'ampar_ez_pick_Ca2PKC_act -> amparp_ez_pick + Ca2PKC_act');
ampar_reac_41b_k = addkineticlaw(ampar_reac_41b,'MassAction');
set(ampar_reac_41b_k, 'ParameterVariableNames', {'kcat_pkc_ampar'});

ampar_reac_41c = addreaction(pfpc_model, 'ampar_ez_pick_CaPKC_AA_act -> amparp_ez_pick + CaPKC_AA_act');
ampar_reac_41c_k = addkineticlaw(ampar_reac_41c,'MassAction');
set(ampar_reac_41c_k, 'ParameterVariableNames', {'kcat_pkc_ampar'});

ampar_reac_41d = addreaction(pfpc_model, 'ampar_ez_pick_Ca2PKC_AA_act -> amparp_ez_pick + Ca2PKC_AA_act');
ampar_reac_41d_k = addkineticlaw(ampar_reac_41d,'MassAction');
set(ampar_reac_41d_k, 'ParameterVariableNames', {'kcat_pkc_ampar'});

ampar_reac_41e = addreaction(pfpc_model, 'ampar_ez_pick_DAG_CaPKC_act -> amparp_ez_pick + DAG_CaPKC_act');
ampar_reac_41e_k = addkineticlaw(ampar_reac_41e,'MassAction');
set(ampar_reac_41e_k, 'ParameterVariableNames', {'kcat_pkc_ampar'});

ampar_reac_41f = addreaction(pfpc_model, 'ampar_ez_pick_DAG_Ca2PKC_act -> amparp_ez_pick + DAG_Ca2PKC_act');
ampar_reac_41f_k = addkineticlaw(ampar_reac_41f,'MassAction');
set(ampar_reac_41f_k, 'ParameterVariableNames', {'kcat_pkc_ampar'});

ampar_reac_41g = addreaction(pfpc_model, 'ampar_ez_pick_DAG_CaPKC_AA_act -> amparp_ez_pick + DAG_CaPKC_AA_act');
ampar_reac_41g_k = addkineticlaw(ampar_reac_41g,'MassAction');
set(ampar_reac_41g_k, 'ParameterVariableNames', {'kcat_pkc_ampar'});

ampar_reac_41h = addreaction(pfpc_model, 'ampar_ez_pick_DAG_Ca2PKC_AA_act -> amparp_ez_pick + DAG_Ca2PKC_AA_act');
ampar_reac_41h_k = addkineticlaw(ampar_reac_41h,'MassAction');
set(ampar_reac_41h_k, 'ParameterVariableNames', {'kcat_pkc_ampar'});

ampar_reac_42 = addreaction(pfpc_model, 'amparp + PP2A -> amparp_pp2a');
ampar_reac_42_k = addkineticlaw(ampar_reac_42,'MassAction');
set(ampar_reac_42_k, 'ParameterVariableNames', {'kon_ampar_pp2a'});

ampar_reac_43 = addreaction(pfpc_model, 'amparp_pp2a -> amparp + PP2A');
ampar_reac_43_k = addkineticlaw(ampar_reac_43,'MassAction');
set(ampar_reac_43_k, 'ParameterVariableNames', {'koff_ampar_pp2a'});

ampar_reac_44 = addreaction(pfpc_model, 'amparp_pp2a -> ampar + PP2A');
ampar_reac_44_k = addkineticlaw(ampar_reac_44,'MassAction');
set(ampar_reac_44_k, 'ParameterVariableNames', {'kcat_pp2a_ampar'});

ampar_reac_45 = addreaction(pfpc_model, 'amparp_pick + PP2A -> amparp_pick_pp2a');
ampar_reac_45_k = addkineticlaw(ampar_reac_45,'MassAction');
set(ampar_reac_45_k, 'ParameterVariableNames', {'kon_ampar_pp2a'});

ampar_reac_46 = addreaction(pfpc_model, 'amparp_pick_pp2a -> amparp_pick + PP2A');
ampar_reac_46_k = addkineticlaw(ampar_reac_46,'MassAction');
set(ampar_reac_46_k, 'ParameterVariableNames', {'koff_ampar_pp2a'});

ampar_reac_47 = addreaction(pfpc_model, 'amparp_pick_pp2a -> ampar_pick + PP2A');
ampar_reac_47_k = addkineticlaw(ampar_reac_47,'MassAction');
set(ampar_reac_47_k, 'ParameterVariableNames', {'kcat_pp2a_ampar'});

ampar_reac_48 = addreaction(pfpc_model, 'amparp_s + PP2A -> amparp_s_pp2a');
ampar_reac_48_k = addkineticlaw(ampar_reac_48,'MassAction');
set(ampar_reac_48_k, 'ParameterVariableNames', {'kon_ampar_pp2a'});

ampar_reac_49 = addreaction(pfpc_model, 'amparp_s_pp2a -> amparp_s + PP2A');
ampar_reac_49_k = addkineticlaw(ampar_reac_49,'MassAction');
set(ampar_reac_49_k, 'ParameterVariableNames', {'koff_ampar_pp2a'});

ampar_reac_50 = addreaction(pfpc_model, 'amparp_s_pp2a -> ampar_s + PP2A');
ampar_reac_50_k = addkineticlaw(ampar_reac_50,'MassAction');
set(ampar_reac_50_k, 'ParameterVariableNames', {'kcat_pp2a_ampar'});

ampar_reac_51 = addreaction(pfpc_model, 'amparp_s_pick + PP2A -> amparp_s_pick_pp2a');
ampar_reac_51_k = addkineticlaw(ampar_reac_51,'MassAction');
set(ampar_reac_51_k, 'ParameterVariableNames', {'kon_ampar_pp2a'});

ampar_reac_52 = addreaction(pfpc_model, 'amparp_s_pick_pp2a -> amparp_s_pick + PP2A');
ampar_reac_52_k = addkineticlaw(ampar_reac_52,'MassAction');
set(ampar_reac_52_k, 'ParameterVariableNames', {'koff_ampar_pp2a'});

ampar_reac_53 = addreaction(pfpc_model, 'amparp_s_pick_pp2a -> ampar_pick + PP2A');
ampar_reac_53_k = addkineticlaw(ampar_reac_53,'MassAction');
set(ampar_reac_53_k, 'ParameterVariableNames', {'kcat_pp2a_ampar'});

ampar_reac_54 = addreaction(pfpc_model, 'amparp_ez + PP2A -> amparp_ez_pp2a');
ampar_reac_54_k = addkineticlaw(ampar_reac_54,'MassAction');
set(ampar_reac_54_k, 'ParameterVariableNames', {'kon_ampar_pp2a'});

ampar_reac_55 = addreaction(pfpc_model, 'amparp_ez_pp2a -> amparp_ez + PP2A');
ampar_reac_55_k = addkineticlaw(ampar_reac_55,'MassAction');
set(ampar_reac_55_k, 'ParameterVariableNames', {'koff_ampar_pp2a'});

ampar_reac_56 = addreaction(pfpc_model, 'amparp_ez_pp2a -> ampar_ez + PP2A');
ampar_reac_56_k = addkineticlaw(ampar_reac_56,'MassAction');
set(ampar_reac_56_k, 'ParameterVariableNames', {'kcat_pp2a_ampar'});

ampar_reac_57 = addreaction(pfpc_model, 'amparp_ez_pick + PP2A -> amparp_ez_pick_pp2a');
ampar_reac_57_k = addkineticlaw(ampar_reac_57,'MassAction');
set(ampar_reac_57_k, 'ParameterVariableNames', {'kon_ampar_pp2a'});

ampar_reac_58 = addreaction(pfpc_model, 'amparp_ez_pick_pp2a -> amparp_ez_pick + PP2A');
ampar_reac_58_k = addkineticlaw(ampar_reac_58,'MassAction');
set(ampar_reac_58_k, 'ParameterVariableNames', {'koff_ampar_pp2a'});

ampar_reac_59 = addreaction(pfpc_model, 'amparp_ez_pick_pp2a -> ampar_ez_pick + PP2A');
ampar_reac_59_k = addkineticlaw(ampar_reac_59,'MassAction');
set(ampar_reac_59_k, 'ParameterVariableNames', {'kcat_pp2a_ampar'});


ampar_reac_42a = addreaction(pfpc_model, 'amparp + PP1 -> amparp_pp1');
ampar_reac_42a_k = addkineticlaw(ampar_reac_42a,'MassAction');
set(ampar_reac_42a_k, 'ParameterVariableNames', {'kon_ampar_pp1'});

ampar_reac_43a = addreaction(pfpc_model, 'amparp_pp1 -> amparp + PP1');
ampar_reac_43a_k = addkineticlaw(ampar_reac_43a,'MassAction');
set(ampar_reac_43a_k, 'ParameterVariableNames', {'koff_ampar_pp1'});

ampar_reac_44a = addreaction(pfpc_model, 'amparp_pp1 -> ampar + PP1');
ampar_reac_44a_k = addkineticlaw(ampar_reac_44a,'MassAction');
set(ampar_reac_44a_k, 'ParameterVariableNames', {'kcat_pp1_ampar'});

ampar_reac_45a = addreaction(pfpc_model, 'amparp_pick + PP1 -> amparp_pick_pp1');
ampar_reac_45a_k = addkineticlaw(ampar_reac_45a,'MassAction');
set(ampar_reac_45a_k, 'ParameterVariableNames', {'kon_ampar_pp1'});

ampar_reac_46a = addreaction(pfpc_model, 'amparp_pick_pp1 -> amparp_pick + PP1');
ampar_reac_46a_k = addkineticlaw(ampar_reac_46a,'MassAction');
set(ampar_reac_46a_k, 'ParameterVariableNames', {'koff_ampar_pp1'});

ampar_reac_47a = addreaction(pfpc_model, 'amparp_pick_pp1 -> ampar_pick + PP1');
ampar_reac_47a_k = addkineticlaw(ampar_reac_47a,'MassAction');
set(ampar_reac_47a_k, 'ParameterVariableNames', {'kcat_pp1_ampar'});

ampar_reac_48a = addreaction(pfpc_model, 'amparp_s + PP1 -> amparp_s_pp1');
ampar_reac_48a_k = addkineticlaw(ampar_reac_48a,'MassAction');
set(ampar_reac_48a_k, 'ParameterVariableNames', {'kon_ampar_pp1'});

ampar_reac_49a = addreaction(pfpc_model, 'amparp_s_pp1 -> amparp_s + PP1');
ampar_reac_49a_k = addkineticlaw(ampar_reac_49a,'MassAction');
set(ampar_reac_49a_k, 'ParameterVariableNames', {'koff_ampar_pp1'});

ampar_reac_50a = addreaction(pfpc_model, 'amparp_s_pp1 -> ampar_s + PP1');
ampar_reac_50a_k = addkineticlaw(ampar_reac_50a,'MassAction');
set(ampar_reac_50a_k, 'ParameterVariableNames', {'kcat_pp1_ampar'});

ampar_reac_51a = addreaction(pfpc_model, 'amparp_s_pick + PP1 -> amparp_s_pick_pp1');
ampar_reac_51a_k = addkineticlaw(ampar_reac_51a,'MassAction');
set(ampar_reac_51a_k, 'ParameterVariableNames', {'kon_ampar_pp1'});

ampar_reac_52a = addreaction(pfpc_model, 'amparp_s_pick_pp1 -> amparp_s_pick + PP1');
ampar_reac_52a_k = addkineticlaw(ampar_reac_52a,'MassAction');
set(ampar_reac_52a_k, 'ParameterVariableNames', {'koff_ampar_pp1'});

ampar_reac_53a = addreaction(pfpc_model, 'amparp_s_pick_pp1 -> ampar_pick + PP1');
ampar_reac_53a_k = addkineticlaw(ampar_reac_53a,'MassAction');
set(ampar_reac_53a_k, 'ParameterVariableNames', {'kcat_pp1_ampar'});

ampar_reac_54a = addreaction(pfpc_model, 'amparp_ez + PP1 -> amparp_ez_pp1');
ampar_reac_54a_k = addkineticlaw(ampar_reac_54a,'MassAction');
set(ampar_reac_54a_k, 'ParameterVariableNames', {'kon_ampar_pp1'});

ampar_reac_55a = addreaction(pfpc_model, 'amparp_ez_pp1 -> amparp_ez + PP1');
ampar_reac_55a_k = addkineticlaw(ampar_reac_55a,'MassAction');
set(ampar_reac_55a_k, 'ParameterVariableNames', {'koff_ampar_pp1'});

ampar_reac_56a = addreaction(pfpc_model, 'amparp_ez_pp1 -> ampar_ez + PP1');
ampar_reac_56a_k = addkineticlaw(ampar_reac_56a,'MassAction');
set(ampar_reac_56a_k, 'ParameterVariableNames', {'kcat_pp1_ampar'});

ampar_reac_57a = addreaction(pfpc_model, 'amparp_ez_pick + PP1 -> amparp_ez_pick_pp1');
ampar_reac_57a_k = addkineticlaw(ampar_reac_57a,'MassAction');
set(ampar_reac_57a_k, 'ParameterVariableNames', {'kon_ampar_pp1'});

ampar_reac_58a = addreaction(pfpc_model, 'amparp_ez_pick_pp1 -> amparp_ez_pick + PP1');
ampar_reac_58a_k = addkineticlaw(ampar_reac_58a,'MassAction');
set(ampar_reac_58a_k, 'ParameterVariableNames', {'koff_ampar_pp1'});

ampar_reac_59a = addreaction(pfpc_model, 'amparp_ez_pick_pp1 -> ampar_ez_pick + PP1');
ampar_reac_59a_k = addkineticlaw(ampar_reac_59a,'MassAction');
set(ampar_reac_59a_k, 'ParameterVariableNames', {'kcat_pp1_ampar'});


ampar_reac_60 = addreaction(pfpc_model, 'ampar -> ampar_s');
ampar_reac_60_k = addkineticlaw(ampar_reac_60,'MassAction');
set(ampar_reac_60_k, 'ParameterVariableNames', {'kdiff_psd_syn'});

ampar_reac_61 = addreaction(pfpc_model, 'ampar_s -> ampar');
ampar_reac_61_k = addkineticlaw(ampar_reac_61,'MassAction');
set(ampar_reac_61_k, 'ParameterVariableNames', {'kdiff_syn_psd'});

ampar_reac_62 = addreaction(pfpc_model, 'amparp -> amparp_s');
ampar_reac_62_k = addkineticlaw(ampar_reac_62,'MassAction');
set(ampar_reac_62_k, 'ParameterVariableNames', {'kdiff_psd_syn'});

ampar_reac_63 = addreaction(pfpc_model, 'amparp_s -> amparp');
ampar_reac_63_k = addkineticlaw(ampar_reac_63,'MassAction');
set(ampar_reac_63_k, 'ParameterVariableNames', {'kdiff_syn_psd'});

ampar_reac_64 = addreaction(pfpc_model, 'ampar_s -> ampar_ez');
ampar_reac_64_k = addkineticlaw(ampar_reac_64,'MassAction');
set(ampar_reac_64_k, 'ParameterVariableNames', {'kdiff_syn_ez'});

ampar_reac_65 = addreaction(pfpc_model, 'ampar_ez -> ampar_s');
ampar_reac_65_k = addkineticlaw(ampar_reac_65,'MassAction');
set(ampar_reac_65_k, 'ParameterVariableNames', {'kdiff_ez_syn'});

ampar_reac_66 = addreaction(pfpc_model, 'amparp_s -> amparp_ez');
ampar_reac_66_k = addkineticlaw(ampar_reac_66,'MassAction');
set(ampar_reac_66_k, 'ParameterVariableNames', {'kdiff_syn_ez'});

ampar_reac_67 = addreaction(pfpc_model, 'amparp_ez -> amparp_s');
ampar_reac_67_k = addkineticlaw(ampar_reac_67,'MassAction');
set(ampar_reac_67_k, 'ParameterVariableNames', {'kdiff_ez_syn'});

ampar_reac_68 = addreaction(pfpc_model, 'ampar_c_pick -> ampar_c + pick');
ampar_reac_68_k = addkineticlaw(ampar_reac_68,'MassAction');
set(ampar_reac_68_k, 'ParameterVariableNames', {'koff_ampar_pick'});

ampar_reac_69 = addreaction(pfpc_model, 'amparp_c_pick -> amparp_c + pick');
ampar_reac_69_k = addkineticlaw(ampar_reac_69,'MassAction');
set(ampar_reac_69_k, 'ParameterVariableNames', {'koff_ampar_pick'});

ampar_reac_70 = addreaction(pfpc_model, 'ampar_pick -> ampar_s_pick');
ampar_reac_70_k = addkineticlaw(ampar_reac_70,'MassAction');
set(ampar_reac_70_k, 'ParameterVariableNames', {'kdiff_psd_syn'});

ampar_reac_71 = addreaction(pfpc_model, 'ampar_s_pick -> ampar_pick');
ampar_reac_71_k = addkineticlaw(ampar_reac_71,'MassAction');
set(ampar_reac_71_k, 'ParameterVariableNames', {'kdiff_syn_psd'});

ampar_reac_72 = addreaction(pfpc_model, 'amparp_pick -> amparp_s_pick');
ampar_reac_72_k = addkineticlaw(ampar_reac_72,'MassAction');
set(ampar_reac_72_k, 'ParameterVariableNames', {'kdiff_psd_syn'});

ampar_reac_73 = addreaction(pfpc_model, 'amparp_s_pick -> amparp_pick');
ampar_reac_73_k = addkineticlaw(ampar_reac_73,'MassAction');
set(ampar_reac_73_k, 'ParameterVariableNames', {'kdiff_syn_psd'});

ampar_reac_74 = addreaction(pfpc_model, 'ampar_s_pick -> ampar_ez_pick');
ampar_reac_74_k = addkineticlaw(ampar_reac_74,'MassAction');
set(ampar_reac_74_k, 'ParameterVariableNames', {'kdiff_syn_ez'});

ampar_reac_75 = addreaction(pfpc_model, 'ampar_ez_pick -> ampar_s_pick');
ampar_reac_75_k = addkineticlaw(ampar_reac_75,'MassAction');
set(ampar_reac_75_k, 'ParameterVariableNames', {'kdiff_ez_syn'});

ampar_reac_76 = addreaction(pfpc_model, 'amparp_s_pick -> amparp_ez_pick');
ampar_reac_76_k = addkineticlaw(ampar_reac_76,'MassAction');
set(ampar_reac_76_k, 'ParameterVariableNames', {'kdiff_syn_ez'});

ampar_reac_77 = addreaction(pfpc_model, 'amparp_ez_pick -> amparp_s_pick');
ampar_reac_77_k = addkineticlaw(ampar_reac_77,'MassAction');
set(ampar_reac_77_k, 'ParameterVariableNames', {'kdiff_ez_syn'});



ampar_reac_78 = addreaction(pfpc_model, 'ampar_pick + nsf -> ampar_pick_nsf');
ampar_reac_78_k = addkineticlaw(ampar_reac_78,'MassAction');
set(ampar_reac_78_k, 'ParameterVariableNames', {'kon_ampar_nsf'});

ampar_reac_78a = addreaction(pfpc_model, 'ampar_pick + nsf_no -> ampar_pick_nsf_no');
ampar_reac_78a_k = addkineticlaw(ampar_reac_78a,'MassAction');
set(ampar_reac_78a_k, 'ParameterVariableNames', {'kon_ampar_nsf'});

ampar_reac_78b = addreaction(pfpc_model, 'ampar_pick + nsf_2no -> ampar_pick_nsf_2no');
ampar_reac_78b_k = addkineticlaw(ampar_reac_78b,'MassAction');
set(ampar_reac_78b_k, 'ParameterVariableNames', {'kon_ampar_nsf'});

ampar_reac_78c = addreaction(pfpc_model, 'ampar_pick + nsf_3no -> ampar_pick_nsf_3no');
ampar_reac_78c_k = addkineticlaw(ampar_reac_78c,'MassAction');
set(ampar_reac_78c_k, 'ParameterVariableNames', {'kon_ampar_nsf_no'});

ampar_reac_79 = addreaction(pfpc_model, 'ampar_pick_nsf -> ampar_pick + nsf');
ampar_reac_79_k = addkineticlaw(ampar_reac_79,'MassAction');
set(ampar_reac_79_k, 'ParameterVariableNames', {'koff_ampar_nsf'});

ampar_reac_79a = addreaction(pfpc_model, 'ampar_pick_nsf_no -> ampar_pick + nsf_no');
ampar_reac_79a_k = addkineticlaw(ampar_reac_79a,'MassAction');
set(ampar_reac_79a_k, 'ParameterVariableNames', {'koff_ampar_nsf'});

ampar_reac_79b = addreaction(pfpc_model, 'ampar_pick_nsf_2no -> ampar_pick + nsf_2no');
ampar_reac_79b_k = addkineticlaw(ampar_reac_79b,'MassAction');
set(ampar_reac_79b_k, 'ParameterVariableNames', {'koff_ampar_nsf'});

ampar_reac_79c = addreaction(pfpc_model, 'ampar_pick_nsf_3no -> ampar_pick + nsf_3no');
ampar_reac_79c_k = addkineticlaw(ampar_reac_79c,'MassAction');
set(ampar_reac_79c_k, 'ParameterVariableNames', {'koff_ampar_nsf'});

ampar_reac_80 = addreaction(pfpc_model, 'ampar_pick_nsf -> ampar + pick + nsf');
ampar_reac_80_k = addkineticlaw(ampar_reac_80,'MassAction');
set(ampar_reac_80_k, 'ParameterVariableNames', {'kcat_nsf'});

ampar_reac_80a = addreaction(pfpc_model, 'ampar_pick_nsf_no -> ampar + pick + nsf_no');
ampar_reac_80a_k = addkineticlaw(ampar_reac_80a,'MassAction');
set(ampar_reac_80a_k, 'ParameterVariableNames', {'kcat_nsf'});

ampar_reac_80b = addreaction(pfpc_model, 'ampar_pick_nsf_2no -> ampar + pick + nsf_2no');
ampar_reac_80b_k = addkineticlaw(ampar_reac_80b,'MassAction');
set(ampar_reac_80b_k, 'ParameterVariableNames', {'kcat_nsf'});

ampar_reac_80c = addreaction(pfpc_model, 'ampar_pick_nsf_3no -> ampar + pick + nsf_3no');
ampar_reac_80c_k = addkineticlaw(ampar_reac_80c,'MassAction');
set(ampar_reac_80c_k, 'ParameterVariableNames', {'kcat_nsf_no'});



ampar_reac_81 = addreaction(pfpc_model, 'amparp_pick + nsf -> amparp_pick_nsf');
ampar_reac_81_k = addkineticlaw(ampar_reac_81,'MassAction');
set(ampar_reac_81_k, 'ParameterVariableNames', {'kon_ampar_nsf'});

ampar_reac_81a = addreaction(pfpc_model, 'amparp_pick + nsf_no -> amparp_pick_nsf_no');
ampar_reac_81a_k = addkineticlaw(ampar_reac_81a,'MassAction');
set(ampar_reac_81a_k, 'ParameterVariableNames', {'kon_ampar_nsf'});

ampar_reac_81b = addreaction(pfpc_model, 'amparp_pick + nsf_2no -> amparp_pick_nsf_2no');
ampar_reac_81b_k = addkineticlaw(ampar_reac_81b,'MassAction');
set(ampar_reac_81b_k, 'ParameterVariableNames', {'kon_ampar_nsf'});

ampar_reac_81c = addreaction(pfpc_model, 'amparp_pick + nsf_3no -> amparp_pick_nsf_3no');
ampar_reac_81c_k = addkineticlaw(ampar_reac_81c,'MassAction');
set(ampar_reac_81c_k, 'ParameterVariableNames', {'kon_ampar_nsf_no'});

ampar_reac_82 = addreaction(pfpc_model, 'amparp_pick_nsf -> amparp_pick + nsf');
ampar_reac_82_k = addkineticlaw(ampar_reac_82,'MassAction');
set(ampar_reac_82_k, 'ParameterVariableNames', {'koff_ampar_nsf'});

ampar_reac_82a = addreaction(pfpc_model, 'amparp_pick_nsf_no -> amparp_pick + nsf_no');
ampar_reac_82a_k = addkineticlaw(ampar_reac_82a,'MassAction');
set(ampar_reac_82a_k, 'ParameterVariableNames', {'koff_ampar_nsf'});

ampar_reac_82b = addreaction(pfpc_model, 'amparp_pick_nsf_2no -> amparp_pick + nsf_2no');
ampar_reac_82b_k = addkineticlaw(ampar_reac_82b,'MassAction');
set(ampar_reac_82b_k, 'ParameterVariableNames', {'koff_ampar_nsf'});

ampar_reac_82c = addreaction(pfpc_model, 'amparp_pick_nsf_3no -> amparp_pick + nsf_3no');
ampar_reac_82c_k = addkineticlaw(ampar_reac_82c,'MassAction');
set(ampar_reac_82c_k, 'ParameterVariableNames', {'koff_ampar_nsf'});

ampar_reac_83 = addreaction(pfpc_model, 'amparp_pick_nsf -> amparp + pick + nsf');
ampar_reac_83_k = addkineticlaw(ampar_reac_83,'MassAction');
set(ampar_reac_83_k, 'ParameterVariableNames', {'kcat_nsf'});

ampar_reac_83a = addreaction(pfpc_model, 'amparp_pick_nsf_no -> amparp + pick + nsf_no');
ampar_reac_83a_k = addkineticlaw(ampar_reac_83a,'MassAction');
set(ampar_reac_83a_k, 'ParameterVariableNames', {'kcat_nsf'});

ampar_reac_83b = addreaction(pfpc_model, 'amparp_pick_nsf_2no -> amparp + pick + nsf_2no');
ampar_reac_83b_k = addkineticlaw(ampar_reac_83b,'MassAction');
set(ampar_reac_83b_k, 'ParameterVariableNames', {'kcat_nsf'});

ampar_reac_83c = addreaction(pfpc_model, 'amparp_pick_nsf_3no -> amparp + pick + nsf_3no');
ampar_reac_83c_k = addkineticlaw(ampar_reac_83c,'MassAction');
set(ampar_reac_83c_k, 'ParameterVariableNames', {'kcat_nsf_no'});



ampar_reac_84 = addreaction(pfpc_model, 'ampar_s_pick + nsf -> ampar_s_pick_nsf');
ampar_reac_84_k = addkineticlaw(ampar_reac_84,'MassAction');
set(ampar_reac_84_k, 'ParameterVariableNames', {'kon_ampar_nsf'});

ampar_reac_84a = addreaction(pfpc_model, 'ampar_s_pick + nsf_no -> ampar_s_pick_nsf_no');
ampar_reac_84a_k = addkineticlaw(ampar_reac_84a,'MassAction');
set(ampar_reac_84a_k, 'ParameterVariableNames', {'kon_ampar_nsf'});

ampar_reac_84b = addreaction(pfpc_model, 'ampar_s_pick + nsf_2no -> ampar_s_pick_nsf_2no');
ampar_reac_84b_k = addkineticlaw(ampar_reac_84b,'MassAction');
set(ampar_reac_84b_k, 'ParameterVariableNames', {'kon_ampar_nsf'});

ampar_reac_84c = addreaction(pfpc_model, 'ampar_s_pick + nsf_3no -> ampar_s_pick_nsf_3no');
ampar_reac_84c_k = addkineticlaw(ampar_reac_84c,'MassAction');
set(ampar_reac_84c_k, 'ParameterVariableNames', {'kon_ampar_nsf_no'});

ampar_reac_85 = addreaction(pfpc_model, 'ampar_s_pick_nsf -> ampar_s_pick + nsf');
ampar_reac_85_k = addkineticlaw(ampar_reac_85,'MassAction');
set(ampar_reac_85_k, 'ParameterVariableNames', {'koff_ampar_nsf'});

ampar_reac_85a = addreaction(pfpc_model, 'ampar_s_pick_nsf_3no -> ampar_s_pick + nsf_3no');
ampar_reac_85a_k = addkineticlaw(ampar_reac_85a,'MassAction');
set(ampar_reac_85a_k, 'ParameterVariableNames', {'koff_ampar_nsf'});

ampar_reac_85b = addreaction(pfpc_model, 'ampar_s_pick_nsf_no -> ampar_s_pick + nsf_no');
ampar_reac_85b_k = addkineticlaw(ampar_reac_85b,'MassAction');
set(ampar_reac_85b_k, 'ParameterVariableNames', {'koff_ampar_nsf'});

ampar_reac_85c = addreaction(pfpc_model, 'ampar_s_pick_nsf_2no -> ampar_s_pick + nsf_2no');
ampar_reac_85c_k = addkineticlaw(ampar_reac_85c,'MassAction');
set(ampar_reac_85c_k, 'ParameterVariableNames', {'koff_ampar_nsf'});

ampar_reac_86 = addreaction(pfpc_model, 'ampar_s_pick_nsf -> ampar_s + pick + nsf');
ampar_reac_86_k = addkineticlaw(ampar_reac_86,'MassAction');
set(ampar_reac_86_k, 'ParameterVariableNames', {'kcat_nsf'});

ampar_reac_86a = addreaction(pfpc_model, 'ampar_s_pick_nsf_no -> ampar_s + pick + nsf_no');
ampar_reac_86a_k = addkineticlaw(ampar_reac_86a,'MassAction');
set(ampar_reac_86a_k, 'ParameterVariableNames', {'kcat_nsf'});

ampar_reac_86b = addreaction(pfpc_model, 'ampar_s_pick_nsf_2no -> ampar_s + pick + nsf_2no');
ampar_reac_86b_k = addkineticlaw(ampar_reac_86b,'MassAction');
set(ampar_reac_86b_k, 'ParameterVariableNames', {'kcat_nsf'});

ampar_reac_86c = addreaction(pfpc_model, 'ampar_s_pick_nsf_3no -> ampar_s + pick + nsf_3no');
ampar_reac_86c_k = addkineticlaw(ampar_reac_86c,'MassAction');
set(ampar_reac_86c_k, 'ParameterVariableNames', {'kcat_nsf_no'});



ampar_reac_87 = addreaction(pfpc_model, 'amparp_s_pick + nsf -> amparp_s_pick_nsf');
ampar_reac_87_k = addkineticlaw(ampar_reac_87,'MassAction');
set(ampar_reac_87_k, 'ParameterVariableNames', {'kon_ampar_nsf'});

ampar_reac_87a = addreaction(pfpc_model, 'amparp_s_pick + nsf_no -> amparp_s_pick_nsf_no');
ampar_reac_87a_k = addkineticlaw(ampar_reac_87a,'MassAction');
set(ampar_reac_87a_k, 'ParameterVariableNames', {'kon_ampar_nsf'});

ampar_reac_87b = addreaction(pfpc_model, 'amparp_s_pick + nsf_2no -> amparp_s_pick_nsf_2no');
ampar_reac_87b_k = addkineticlaw(ampar_reac_87b,'MassAction');
set(ampar_reac_87b_k, 'ParameterVariableNames', {'kon_ampar_nsf'});

ampar_reac_87c = addreaction(pfpc_model, 'amparp_s_pick + nsf_3no -> amparp_s_pick_nsf_3no');
ampar_reac_87c_k = addkineticlaw(ampar_reac_87c,'MassAction');
set(ampar_reac_87c_k, 'ParameterVariableNames', {'kon_ampar_nsf_no'});

ampar_reac_88 = addreaction(pfpc_model, 'amparp_s_pick_nsf -> amparp_s_pick + nsf');
ampar_reac_88_k = addkineticlaw(ampar_reac_88,'MassAction');
set(ampar_reac_88_k, 'ParameterVariableNames', {'koff_ampar_nsf'});

ampar_reac_88a = addreaction(pfpc_model, 'amparp_s_pick_nsf_no -> amparp_s_pick + nsf_no');
ampar_reac_88a_k = addkineticlaw(ampar_reac_88a,'MassAction');
set(ampar_reac_88a_k, 'ParameterVariableNames', {'koff_ampar_nsf'});

ampar_reac_88b = addreaction(pfpc_model, 'amparp_s_pick_nsf_2no -> amparp_s_pick + nsf_2no');
ampar_reac_88b_k = addkineticlaw(ampar_reac_88b,'MassAction');
set(ampar_reac_88b_k, 'ParameterVariableNames', {'koff_ampar_nsf'});

ampar_reac_88c = addreaction(pfpc_model, 'amparp_s_pick_nsf_3no -> amparp_s_pick + nsf_3no');
ampar_reac_88c_k = addkineticlaw(ampar_reac_88c,'MassAction');
set(ampar_reac_88c_k, 'ParameterVariableNames', {'koff_ampar_nsf'});

ampar_reac_89 = addreaction(pfpc_model, 'amparp_s_pick_nsf -> amparp_s + pick + nsf');
ampar_reac_89_k = addkineticlaw(ampar_reac_89,'MassAction');
set(ampar_reac_89_k, 'ParameterVariableNames', {'kcat_nsf'});

ampar_reac_89a = addreaction(pfpc_model, 'amparp_s_pick_nsf_no -> amparp_s + pick + nsf_no');
ampar_reac_89a_k = addkineticlaw(ampar_reac_89a,'MassAction');
set(ampar_reac_89a_k, 'ParameterVariableNames', {'kcat_nsf'});

ampar_reac_89b = addreaction(pfpc_model, 'amparp_s_pick_nsf_2no -> amparp_s + pick + nsf_2no');
ampar_reac_89b_k = addkineticlaw(ampar_reac_89b,'MassAction');
set(ampar_reac_89b_k, 'ParameterVariableNames', {'kcat_nsf'});

ampar_reac_89c = addreaction(pfpc_model, 'amparp_s_pick_nsf_3no -> amparp_s + pick + nsf_3no');
ampar_reac_89c_k = addkineticlaw(ampar_reac_89c,'MassAction');
set(ampar_reac_89c_k, 'ParameterVariableNames', {'kcat_nsf_no'});


ampar_reac_90 = addreaction(pfpc_model, 'ampar_ez_pick + nsf -> ampar_ez_pick_nsf');
ampar_reac_90_k = addkineticlaw(ampar_reac_90,'MassAction');
set(ampar_reac_90_k, 'ParameterVariableNames', {'kon_ampar_nsf'});

ampar_reac_90a = addreaction(pfpc_model, 'ampar_ez_pick + nsf_no -> ampar_ez_pick_nsf_no');
ampar_reac_90a_k = addkineticlaw(ampar_reac_90a,'MassAction');
set(ampar_reac_90a_k, 'ParameterVariableNames', {'kon_ampar_nsf'});

ampar_reac_90b = addreaction(pfpc_model, 'ampar_ez_pick + nsf_2no -> ampar_ez_pick_nsf_2no');
ampar_reac_90b_k = addkineticlaw(ampar_reac_90b,'MassAction');
set(ampar_reac_90b_k, 'ParameterVariableNames', {'kon_ampar_nsf'});

ampar_reac_90c = addreaction(pfpc_model, 'ampar_ez_pick + nsf_3no -> ampar_ez_pick_nsf_3no');
ampar_reac_90c_k = addkineticlaw(ampar_reac_90c,'MassAction');
set(ampar_reac_90c_k, 'ParameterVariableNames', {'kon_ampar_nsf_no'});

ampar_reac_91 = addreaction(pfpc_model, 'ampar_ez_pick_nsf -> ampar_ez_pick + nsf');
ampar_reac_91_k = addkineticlaw(ampar_reac_91,'MassAction');
set(ampar_reac_91_k, 'ParameterVariableNames', {'koff_ampar_nsf'});

ampar_reac_91a = addreaction(pfpc_model, 'ampar_ez_pick_nsf_no -> ampar_ez_pick + nsf_no');
ampar_reac_91a_k = addkineticlaw(ampar_reac_91a,'MassAction');
set(ampar_reac_91a_k, 'ParameterVariableNames', {'koff_ampar_nsf'});

ampar_reac_91b = addreaction(pfpc_model, 'ampar_ez_pick_nsf_2no -> ampar_ez_pick + nsf_2no');
ampar_reac_91b_k = addkineticlaw(ampar_reac_91b,'MassAction');
set(ampar_reac_91b_k, 'ParameterVariableNames', {'koff_ampar_nsf'});

ampar_reac_91c = addreaction(pfpc_model, 'ampar_ez_pick_nsf_3no -> ampar_ez_pick + nsf_3no');
ampar_reac_91c_k = addkineticlaw(ampar_reac_91c,'MassAction');
set(ampar_reac_91c_k, 'ParameterVariableNames', {'koff_ampar_nsf'});

ampar_reac_92 = addreaction(pfpc_model, 'ampar_ez_pick_nsf -> ampar_ez + pick + nsf');
ampar_reac_92_k = addkineticlaw(ampar_reac_92,'MassAction');
set(ampar_reac_92_k, 'ParameterVariableNames', {'kcat_nsf'});

ampar_reac_92a = addreaction(pfpc_model, 'ampar_ez_pick_nsf_no -> ampar_ez + pick + nsf_no');
ampar_reac_92a_k = addkineticlaw(ampar_reac_92a,'MassAction');
set(ampar_reac_92a_k, 'ParameterVariableNames', {'kcat_nsf'});

ampar_reac_92b = addreaction(pfpc_model, 'ampar_ez_pick_nsf_2no -> ampar_ez + pick + nsf_2no');
ampar_reac_92b_k = addkineticlaw(ampar_reac_92b,'MassAction');
set(ampar_reac_92b_k, 'ParameterVariableNames', {'kcat_nsf'});

ampar_reac_92c = addreaction(pfpc_model, 'ampar_ez_pick_nsf_3no -> ampar_ez + pick + nsf_3no');
ampar_reac_92c_k = addkineticlaw(ampar_reac_92c,'MassAction');
set(ampar_reac_92c_k, 'ParameterVariableNames', {'kcat_nsf_no'});



ampar_reac_93 = addreaction(pfpc_model, 'amparp_ez_pick + nsf -> amparp_ez_pick_nsf');
ampar_reac_93_k = addkineticlaw(ampar_reac_93,'MassAction');
set(ampar_reac_93_k, 'ParameterVariableNames', {'kon_ampar_nsf'});

ampar_reac_93a = addreaction(pfpc_model, 'amparp_ez_pick + nsf_no -> amparp_ez_pick_nsf_no');
ampar_reac_93a_k = addkineticlaw(ampar_reac_93a,'MassAction');
set(ampar_reac_93a_k, 'ParameterVariableNames', {'kon_ampar_nsf'});

ampar_reac_93b = addreaction(pfpc_model, 'amparp_ez_pick + nsf_2no -> amparp_ez_pick_nsf_2no');
ampar_reac_93b_k = addkineticlaw(ampar_reac_93b,'MassAction');
set(ampar_reac_93b_k, 'ParameterVariableNames', {'kon_ampar_nsf'});

ampar_reac_93c = addreaction(pfpc_model, 'amparp_ez_pick + nsf_3no -> amparp_ez_pick_nsf_3no');
ampar_reac_93c_k = addkineticlaw(ampar_reac_93c,'MassAction');
set(ampar_reac_93c_k, 'ParameterVariableNames', {'kon_ampar_nsf_no'});

ampar_reac_94 = addreaction(pfpc_model, 'amparp_ez_pick_nsf -> amparp_ez_pick + nsf');
ampar_reac_94_k = addkineticlaw(ampar_reac_94,'MassAction');
set(ampar_reac_94_k, 'ParameterVariableNames', {'koff_ampar_nsf'});

ampar_reac_94a = addreaction(pfpc_model, 'amparp_ez_pick_nsf_no -> amparp_ez_pick + nsf_no');
ampar_reac_94a_k = addkineticlaw(ampar_reac_94a,'MassAction');
set(ampar_reac_94a_k, 'ParameterVariableNames', {'koff_ampar_nsf'});

ampar_reac_94b = addreaction(pfpc_model, 'amparp_ez_pick_nsf_2no -> amparp_ez_pick + nsf_2no');
ampar_reac_94b_k = addkineticlaw(ampar_reac_94b,'MassAction');
set(ampar_reac_94b_k, 'ParameterVariableNames', {'koff_ampar_nsf'});

ampar_reac_94c = addreaction(pfpc_model, 'amparp_ez_pick_nsf_3no -> amparp_ez_pick + nsf_3no');
ampar_reac_94c_k = addkineticlaw(ampar_reac_94c,'MassAction');
set(ampar_reac_94c_k, 'ParameterVariableNames', {'koff_ampar_nsf'});

ampar_reac_95 = addreaction(pfpc_model, 'amparp_ez_pick_nsf -> amparp_ez + pick + nsf');
ampar_reac_95_k = addkineticlaw(ampar_reac_95,'MassAction');
set(ampar_reac_95_k, 'ParameterVariableNames', {'kcat_nsf'});

ampar_reac_95a = addreaction(pfpc_model, 'amparp_ez_pick_nsf_no -> amparp_ez + pick + nsf_no');
ampar_reac_95a_k = addkineticlaw(ampar_reac_95a,'MassAction');
set(ampar_reac_95a_k, 'ParameterVariableNames', {'kcat_nsf'});

ampar_reac_95b = addreaction(pfpc_model, 'amparp_ez_pick_nsf_2no -> amparp_ez + pick + nsf_2no');
ampar_reac_95b_k = addkineticlaw(ampar_reac_95b,'MassAction');
set(ampar_reac_95b_k, 'ParameterVariableNames', {'kcat_nsf'});

ampar_reac_95c = addreaction(pfpc_model, 'amparp_ez_pick_nsf_3no -> amparp_ez + pick + nsf_3no');
ampar_reac_95c_k = addkineticlaw(ampar_reac_95c,'MassAction');
set(ampar_reac_95c_k, 'ParameterVariableNames', {'kcat_nsf_no'});



ampar_reac_96 = addreaction(pfpc_model, 'ampar_c_pick + nsf -> ampar_c_pick_nsf');
ampar_reac_96_k = addkineticlaw(ampar_reac_96,'MassAction');
set(ampar_reac_96_k, 'ParameterVariableNames', {'kon_ampar_nsf'});

ampar_reac_96a = addreaction(pfpc_model, 'ampar_c_pick + nsf_no -> ampar_c_pick_nsf_no');
ampar_reac_96a_k = addkineticlaw(ampar_reac_96a,'MassAction');
set(ampar_reac_96a_k, 'ParameterVariableNames', {'kon_ampar_nsf'});

ampar_reac_96b = addreaction(pfpc_model, 'ampar_c_pick + nsf_2no -> ampar_c_pick_nsf_2no');
ampar_reac_96b_k = addkineticlaw(ampar_reac_96b,'MassAction');
set(ampar_reac_96b_k, 'ParameterVariableNames', {'kon_ampar_nsf'});

ampar_reac_96c = addreaction(pfpc_model, 'ampar_c_pick + nsf_3no -> ampar_c_pick_nsf_3no');
ampar_reac_96c_k = addkineticlaw(ampar_reac_96c,'MassAction');
set(ampar_reac_96c_k, 'ParameterVariableNames', {'kon_ampar_nsf_no'});

ampar_reac_97 = addreaction(pfpc_model, 'ampar_c_pick_nsf -> ampar_c_pick + nsf');
ampar_reac_97_k = addkineticlaw(ampar_reac_97,'MassAction');
set(ampar_reac_97_k, 'ParameterVariableNames', {'koff_ampar_nsf'});

ampar_reac_97a = addreaction(pfpc_model, 'ampar_c_pick_nsf_no -> ampar_c_pick + nsf_no');
ampar_reac_97a_k = addkineticlaw(ampar_reac_97a,'MassAction');
set(ampar_reac_97a_k, 'ParameterVariableNames', {'koff_ampar_nsf'});

ampar_reac_97b = addreaction(pfpc_model, 'ampar_c_pick_nsf_2no -> ampar_c_pick + nsf_2no');
ampar_reac_97b_k = addkineticlaw(ampar_reac_97b,'MassAction');
set(ampar_reac_97b_k, 'ParameterVariableNames', {'koff_ampar_nsf'});

ampar_reac_97c = addreaction(pfpc_model, 'ampar_c_pick_nsf_3no -> ampar_c_pick + nsf_3no');
ampar_reac_97c_k = addkineticlaw(ampar_reac_97c,'MassAction');
set(ampar_reac_97c_k, 'ParameterVariableNames', {'koff_ampar_nsf'});

ampar_reac_98 = addreaction(pfpc_model, 'ampar_c_pick_nsf -> ampar_c + pick + nsf');
ampar_reac_98_k = addkineticlaw(ampar_reac_98,'MassAction');
set(ampar_reac_98_k, 'ParameterVariableNames', {'kcat_nsf'});

ampar_reac_98a = addreaction(pfpc_model, 'ampar_c_pick_nsf_no -> ampar_c + pick + nsf_no');
ampar_reac_98a_k = addkineticlaw(ampar_reac_98a,'MassAction');
set(ampar_reac_98a_k, 'ParameterVariableNames', {'kcat_nsf'});

ampar_reac_98b = addreaction(pfpc_model, 'ampar_c_pick_nsf_2no -> ampar_c + pick + nsf_2no');
ampar_reac_98b_k = addkineticlaw(ampar_reac_98b,'MassAction');
set(ampar_reac_98b_k, 'ParameterVariableNames', {'kcat_nsf'});

ampar_reac_98c = addreaction(pfpc_model, 'ampar_c_pick_nsf_3no -> ampar_c + pick + nsf_3no');
ampar_reac_98c_k = addkineticlaw(ampar_reac_98c,'MassAction');
set(ampar_reac_98c_k, 'ParameterVariableNames', {'kcat_nsf_no'});



ampar_reac_99 = addreaction(pfpc_model, 'amparp_c_pick + nsf -> amparp_c_pick_nsf');
ampar_reac_99_k = addkineticlaw(ampar_reac_99,'MassAction');
set(ampar_reac_99_k, 'ParameterVariableNames', {'kon_ampar_nsf'});

ampar_reac_99a = addreaction(pfpc_model, 'amparp_c_pick + nsf_no -> amparp_c_pick_nsf_no');
ampar_reac_99a_k = addkineticlaw(ampar_reac_99a,'MassAction');
set(ampar_reac_99a_k, 'ParameterVariableNames', {'kon_ampar_nsf'});

ampar_reac_99b = addreaction(pfpc_model, 'amparp_c_pick + nsf_2no -> amparp_c_pick_nsf_2no');
ampar_reac_99b_k = addkineticlaw(ampar_reac_99b,'MassAction');
set(ampar_reac_99b_k, 'ParameterVariableNames', {'kon_ampar_nsf'});

ampar_reac_99c = addreaction(pfpc_model, 'amparp_c_pick + nsf_3no -> amparp_c_pick_nsf_3no');
ampar_reac_99c_k = addkineticlaw(ampar_reac_99c,'MassAction');
set(ampar_reac_99c_k, 'ParameterVariableNames', {'kon_ampar_nsf_no'});

ampar_reac_100 = addreaction(pfpc_model, 'amparp_c_pick_nsf -> amparp_c_pick + nsf');
ampar_reac_100_k = addkineticlaw(ampar_reac_100,'MassAction');
set(ampar_reac_100_k, 'ParameterVariableNames', {'koff_ampar_nsf'});

ampar_reac_100a = addreaction(pfpc_model, 'amparp_c_pick_nsf_no -> amparp_c_pick + nsf_no');
ampar_reac_100a_k = addkineticlaw(ampar_reac_100a,'MassAction');
set(ampar_reac_100a_k, 'ParameterVariableNames', {'koff_ampar_nsf'});

ampar_reac_100b = addreaction(pfpc_model, 'amparp_c_pick_nsf_2no -> amparp_c_pick + nsf_2no');
ampar_reac_100b_k = addkineticlaw(ampar_reac_100b,'MassAction');
set(ampar_reac_100b_k, 'ParameterVariableNames', {'koff_ampar_nsf'});

ampar_reac_100c = addreaction(pfpc_model, 'amparp_c_pick_nsf_3no -> amparp_c_pick + nsf_3no');
ampar_reac_100c_k = addkineticlaw(ampar_reac_100c,'MassAction');
set(ampar_reac_100c_k, 'ParameterVariableNames', {'koff_ampar_nsf'});

ampar_reac_101 = addreaction(pfpc_model, 'amparp_c_pick_nsf -> amparp_c + pick + nsf');
ampar_reac_101_k = addkineticlaw(ampar_reac_101,'MassAction');
set(ampar_reac_101_k, 'ParameterVariableNames', {'kcat_nsf'});

ampar_reac_101a = addreaction(pfpc_model, 'amparp_c_pick_nsf_no -> amparp_c + pick + nsf_no');
ampar_reac_101a_k = addkineticlaw(ampar_reac_101a,'MassAction');
set(ampar_reac_101a_k, 'ParameterVariableNames', {'kcat_nsf'});

ampar_reac_101b = addreaction(pfpc_model, 'amparp_c_pick_nsf_2no -> amparp_c + pick + nsf_2no');
ampar_reac_101b_k = addkineticlaw(ampar_reac_101b,'MassAction');
set(ampar_reac_101b_k, 'ParameterVariableNames', {'kcat_nsf'});

ampar_reac_101c = addreaction(pfpc_model, 'amparp_c_pick_nsf_3no -> amparp_c + pick + nsf_3no');
ampar_reac_101c_k = addkineticlaw(ampar_reac_101c,'MassAction');
set(ampar_reac_101c_k, 'ParameterVariableNames', {'kcat_nsf_no'});


%%PKC phosphorylation of NSF
%Species

NSFBPpase  = addspecies(cytosol, 'NSFBPpase', 'InitialAmount', 0.3, 'InitialAmountUnits', 'micromolarity');
NSFBP = addspecies(cytosol, 'NSFBP', 'InitialAmount', 4, 'InitialAmountUnits', 'micromolarity');
NSFBPp = addspecies(cytosol, 'NSFBPp', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
NSFBPp_NSFBPpase = addspecies(cytosol, 'NSFBPp_NSFBPpase', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
NSFBPpp_NSFBPpase = addspecies(cytosol, 'NSFBPpp_NSFBPpase', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');

CaPKC_act_NSFBP = addspecies(cytosol, 'CaPKC_act_NSFBP', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
Ca2PKC_act_NSFBP = addspecies(cytosol, 'Ca2PKC_act_NSFBP', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
CaPKC_AA_act_NSFBP = addspecies(cytosol, 'CaPKC_AA_act_NSFBP', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
Ca2PKC_AA_act_NSFBP = addspecies(cytosol, 'Ca2PKC_AA_act_NSFBP', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
DAG_CaPKC_act_NSFBP = addspecies(cytosol, 'DAG_CaPKC_act_NSFBP', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
DAG_Ca2PKC_act_NSFBP = addspecies(cytosol, 'DAG_Ca2PKC_act_NSFBP', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
DAG_CaPKC_AA_act_NSFBP = addspecies(cytosol, 'DAG_CaPKC_AA_act_NSFBP', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
DAG_Ca2PKC_AA_act_NSFBP = addspecies(cytosol, 'DAG_Ca2PKC_AA_act_NSFBP', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');

CaPKC_act_NSFBPp = addspecies(cytosol, 'CaPKC_act_NSFBPp', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
Ca2PKC_act_NSFBPp = addspecies(cytosol, 'Ca2PKC_act_NSFBPp', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
CaPKC_AA_act_NSFBPp = addspecies(cytosol, 'CaPKC_AA_act_NSFBPp', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
Ca2PKC_AA_act_NSFBPp = addspecies(cytosol, 'Ca2PKC_AA_act_NSFBPp', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
DAG_CaPKC_act_NSFBPp = addspecies(cytosol, 'DAG_CaPKC_act_NSFBPp', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
DAG_Ca2PKC_act_NSFBPp = addspecies(cytosol, 'DAG_Ca2PKC_act_NSFBPp', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
DAG_CaPKC_AA_act_NSFBPp = addspecies(cytosol, 'DAG_CaPKC_AA_act_NSFBPp', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
DAG_Ca2PKC_AA_act_NSFBPp = addspecies(cytosol, 'DAG_Ca2PKC_AA_act_NSFBPp', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');


nsf_no_NSFBPpp  = addspecies(cytosol, 'nsf_no_NSFBPpp', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
nsf_2no_NSFBPpp  = addspecies(cytosol, 'nsf_2no_NSFBPpp', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
nsf_3no_NSFBPpp  = addspecies(cytosol, 'nsf_3no_NSFBPpp', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');

NSFBPp_tot = addspecies(cytosol, 'NSFBPp_tot', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
NSFBPpp_tot = addspecies(cytosol, 'NSFBPpp_tot', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');

%Parameters
kon_nsf_nsfbp = addparameter(pfpc_model,'kon_nsf_nsfbp', 'Value', 100, 'ValueUnits', '1/(micromolarity*second)');
kon_nsf_nsfbp1 = addparameter(pfpc_model,'kon_nsf_nsfbp1', 'Value', 0, 'ValueUnits', '1/(micromolarity*second)');
koff_nsf_nsfbp = addparameter(pfpc_model,'koff_nsf_nsfbp', 'Value', 2, 'ValueUnits', '1/second');

kon_pkc_nsfbp = addparameter(pfpc_model,'kon_pkc_nsfbp', 'Value', 1, 'ValueUnits', '1/(micromolarity*second)');
koff_pkc_nsfbp = addparameter(pfpc_model,'koff_pkc_nsfbp', 'Value', 40, 'ValueUnits', '1/second');

kon_pkc_nsfbp2 = addparameter(pfpc_model,'kon_pkc_nsfbp2', 'Value', 100, 'ValueUnits', '1/(micromolarity*second)');
koff_pkc_nsfbp2 = addparameter(pfpc_model,'koff_pkc_nsfbp2', 'Value', 0.1, 'ValueUnits', '1/second');

kcat_pkc_nsfbp = addparameter(pfpc_model,'kcat_pkc_nsfbp', 'Value', 0.05, 'ValueUnits', '1/second');
kcat_pkc_nsfbp2 = addparameter(pfpc_model,'kcat_pkc_nsfbp2', 'Value', 0.5, 'ValueUnits', '1/second');


kon_nsfbppase_nsfbp = addparameter(pfpc_model,'kon_nsfbppase_nsfbp', 'Value', 20, 'ValueUnits', '1/(micromolarity*second)');
koff_nsfbppase_nsfbp = addparameter(pfpc_model,'koff_nsfbppase_nsfbp', 'Value', 2, 'ValueUnits', '1/second');
kcat_nsfbppase = addparameter(pfpc_model,'kcat_nsfbppase', 'Value', 10, 'ValueUnits', '1/second');

%Reactions

pkc_nsfbp_reac_1a = addreaction(pfpc_model, 'CaPKC_act + NSFBP -> CaPKC_act_NSFBP');
pkc_nsfbp_reac_1a_k = addkineticlaw(pkc_nsfbp_reac_1a,'MassAction');
set(pkc_nsfbp_reac_1a_k, 'ParameterVariableNames', {'kon_pkc_nsfbp'});

pkc_nsfbp_reac_1b = addreaction(pfpc_model, 'CaPKC_act_NSFBP -> CaPKC_act + NSFBP');
pkc_nsfbp_reac_1b_k = addkineticlaw(pkc_nsfbp_reac_1b,'MassAction');
set(pkc_nsfbp_reac_1b_k, 'ParameterVariableNames', {'koff_pkc_nsfbp'});

pkc_nsfbp_reac_1c = addreaction(pfpc_model, 'CaPKC_act_NSFBP -> CaPKC_act + NSFBPp');
pkc_nsfbp_reac_1c_k = addkineticlaw(pkc_nsfbp_reac_1c,'MassAction');
set(pkc_nsfbp_reac_1c_k, 'ParameterVariableNames', {'kcat_pkc_nsfbp'});


pkc_nsfbp_reac_2a = addreaction(pfpc_model, 'Ca2PKC_act + NSFBP -> Ca2PKC_act_NSFBP');
pkc_nsfbp_reac_2a_k = addkineticlaw(pkc_nsfbp_reac_2a,'MassAction');
set(pkc_nsfbp_reac_2a_k, 'ParameterVariableNames', {'kon_pkc_nsfbp'});

pkc_nsfbp_reac_2b = addreaction(pfpc_model, 'Ca2PKC_act_NSFBP -> Ca2PKC_act + NSFBP');
pkc_nsfbp_reac_2b_k = addkineticlaw(pkc_nsfbp_reac_2b,'MassAction');
set(pkc_nsfbp_reac_2b_k, 'ParameterVariableNames', {'koff_pkc_nsfbp'});

pkc_nsfbp_reac_2c = addreaction(pfpc_model, 'Ca2PKC_act_NSFBP -> Ca2PKC_act + NSFBPp');
pkc_nsfbp_reac_2c_k = addkineticlaw(pkc_nsfbp_reac_2c,'MassAction');
set(pkc_nsfbp_reac_2c_k, 'ParameterVariableNames', {'kcat_pkc_nsfbp'});


pkc_nsfbp_reac_3a = addreaction(pfpc_model, 'CaPKC_AA_act + NSFBP -> CaPKC_AA_act_NSFBP');
pkc_nsfbp_reac_3a_k = addkineticlaw(pkc_nsfbp_reac_3a,'MassAction');
set(pkc_nsfbp_reac_3a_k, 'ParameterVariableNames', {'kon_pkc_nsfbp'});

pkc_nsfbp_reac_3b = addreaction(pfpc_model, 'CaPKC_AA_act_NSFBP -> CaPKC_AA_act + NSFBP');
pkc_nsfbp_reac_3b_k = addkineticlaw(pkc_nsfbp_reac_3b,'MassAction');
set(pkc_nsfbp_reac_3b_k, 'ParameterVariableNames', {'koff_pkc_nsfbp'});

pkc_nsfbp_reac_3c = addreaction(pfpc_model, 'CaPKC_AA_act_NSFBP -> CaPKC_AA_act + NSFBPp');
pkc_nsfbp_reac_3c_k = addkineticlaw(pkc_nsfbp_reac_3c,'MassAction');
set(pkc_nsfbp_reac_3c_k, 'ParameterVariableNames', {'kcat_pkc_nsfbp'});


pkc_nsfbp_reac_4a = addreaction(pfpc_model, 'Ca2PKC_AA_act + NSFBP -> Ca2PKC_AA_act_NSFBP');
pkc_nsfbp_reac_4a_k = addkineticlaw(pkc_nsfbp_reac_4a,'MassAction');
set(pkc_nsfbp_reac_4a_k, 'ParameterVariableNames', {'kon_pkc_nsfbp'});

pkc_nsfbp_reac_4b = addreaction(pfpc_model, 'Ca2PKC_AA_act_NSFBP -> Ca2PKC_AA_act + NSFBP');
pkc_nsfbp_reac_4b_k = addkineticlaw(pkc_nsfbp_reac_4b,'MassAction');
set(pkc_nsfbp_reac_4b_k, 'ParameterVariableNames', {'koff_pkc_nsfbp'});

pkc_nsfbp_reac_4c = addreaction(pfpc_model, 'Ca2PKC_AA_act_NSFBP -> Ca2PKC_AA_act + NSFBPp');
pkc_nsfbp_reac_4c_k = addkineticlaw(pkc_nsfbp_reac_4c,'MassAction');
set(pkc_nsfbp_reac_4c_k, 'ParameterVariableNames', {'kcat_pkc_nsfbp'});


pkc_nsfbp_reac_5a = addreaction(pfpc_model, 'DAG_CaPKC_act + NSFBP -> DAG_CaPKC_act_NSFBP');
pkc_nsfbp_reac_5a_k = addkineticlaw(pkc_nsfbp_reac_5a,'MassAction');
set(pkc_nsfbp_reac_5a_k, 'ParameterVariableNames', {'kon_pkc_nsfbp'});

pkc_nsfbp_reac_5b = addreaction(pfpc_model, 'DAG_CaPKC_act_NSFBP -> DAG_CaPKC_act + NSFBP');
pkc_nsfbp_reac_5b_k = addkineticlaw(pkc_nsfbp_reac_5b,'MassAction');
set(pkc_nsfbp_reac_5b_k, 'ParameterVariableNames', {'koff_pkc_nsfbp'});

pkc_nsfbp_reac_5c = addreaction(pfpc_model, 'DAG_CaPKC_act_NSFBP -> DAG_CaPKC_act + NSFBPp');
pkc_nsfbp_reac_5c_k = addkineticlaw(pkc_nsfbp_reac_5c,'MassAction');
set(pkc_nsfbp_reac_5c_k, 'ParameterVariableNames', {'kcat_pkc_nsfbp'});


pkc_nsfbp_reac_6a = addreaction(pfpc_model, 'DAG_Ca2PKC_act + NSFBP -> DAG_Ca2PKC_act_NSFBP');
pkc_nsfbp_reac_6a_k = addkineticlaw(pkc_nsfbp_reac_6a,'MassAction');
set(pkc_nsfbp_reac_6a_k, 'ParameterVariableNames', {'kon_pkc_nsfbp'});

pkc_nsfbp_reac_6b = addreaction(pfpc_model, 'DAG_Ca2PKC_act_NSFBP -> DAG_Ca2PKC_act + NSFBP');
pkc_nsfbp_reac_6b_k = addkineticlaw(pkc_nsfbp_reac_6b,'MassAction');
set(pkc_nsfbp_reac_6b_k, 'ParameterVariableNames', {'koff_pkc_nsfbp'});

pkc_nsfbp_reac_6c = addreaction(pfpc_model, 'DAG_Ca2PKC_act_NSFBP -> DAG_Ca2PKC_act + NSFBPp');
pkc_nsfbp_reac_6c_k = addkineticlaw(pkc_nsfbp_reac_6c,'MassAction');
set(pkc_nsfbp_reac_6c_k, 'ParameterVariableNames', {'kcat_pkc_nsfbp'});


pkc_nsfbp_reac_7a = addreaction(pfpc_model, 'DAG_CaPKC_AA_act + NSFBP -> DAG_CaPKC_AA_act_NSFBP');
pkc_nsfbp_reac_7a_k = addkineticlaw(pkc_nsfbp_reac_7a,'MassAction');
set(pkc_nsfbp_reac_7a_k, 'ParameterVariableNames', {'kon_pkc_nsfbp'});

pkc_nsfbp_reac_7b = addreaction(pfpc_model, 'DAG_CaPKC_AA_act_NSFBP -> DAG_CaPKC_AA_act + NSFBP');
pkc_nsfbp_reac_7b_k = addkineticlaw(pkc_nsfbp_reac_7b,'MassAction');
set(pkc_nsfbp_reac_7b_k, 'ParameterVariableNames', {'koff_pkc_nsfbp'});

pkc_nsfbp_reac_7c = addreaction(pfpc_model, 'DAG_CaPKC_AA_act_NSFBP -> DAG_CaPKC_AA_act + NSFBPp');
pkc_nsfbp_reac_7c_k = addkineticlaw(pkc_nsfbp_reac_7c,'MassAction');
set(pkc_nsfbp_reac_7c_k, 'ParameterVariableNames', {'kcat_pkc_nsfbp'});


pkc_nsfbp_reac_8a = addreaction(pfpc_model, 'DAG_Ca2PKC_AA_act + NSFBP -> DAG_Ca2PKC_AA_act_NSFBP');
pkc_nsfbp_reac_8a_k = addkineticlaw(pkc_nsfbp_reac_8a,'MassAction');
set(pkc_nsfbp_reac_8a_k, 'ParameterVariableNames', {'kon_pkc_nsfbp'});

pkc_nsfbp_reac_8b = addreaction(pfpc_model, 'DAG_Ca2PKC_AA_act_NSFBP -> DAG_Ca2PKC_AA_act + NSFBP');
pkc_nsfbp_reac_8b_k = addkineticlaw(pkc_nsfbp_reac_8b,'MassAction');
set(pkc_nsfbp_reac_8b_k, 'ParameterVariableNames', {'koff_pkc_nsfbp'});

pkc_nsfbp_reac_8c = addreaction(pfpc_model, 'DAG_Ca2PKC_AA_act_NSFBP -> DAG_Ca2PKC_AA_act + NSFBPp');
pkc_nsfbp_reac_8c_k = addkineticlaw(pkc_nsfbp_reac_8c,'MassAction');
set(pkc_nsfbp_reac_8c_k, 'ParameterVariableNames', {'kcat_pkc_nsfbp'});


pkc_nsfbp_reac_11a = addreaction(pfpc_model, 'CaPKC_act + NSFBPp -> CaPKC_act_NSFBPp');
pkc_nsfbp_reac_11a_k = addkineticlaw(pkc_nsfbp_reac_11a,'MassAction');
set(pkc_nsfbp_reac_11a_k, 'ParameterVariableNames', {'kon_pkc_nsfbp2'});

pkc_nsfbp_reac_11b = addreaction(pfpc_model, 'CaPKC_act_NSFBPp -> CaPKC_act + NSFBPp');
pkc_nsfbp_reac_11b_k = addkineticlaw(pkc_nsfbp_reac_11b,'MassAction');
set(pkc_nsfbp_reac_11b_k, 'ParameterVariableNames', {'koff_pkc_nsfbp2'});

pkc_nsfbp_reac_11c = addreaction(pfpc_model, 'CaPKC_act_NSFBPp -> CaPKC_act + NSFBPpp');
pkc_nsfbp_reac_11c_k = addkineticlaw(pkc_nsfbp_reac_11c,'MassAction');
set(pkc_nsfbp_reac_11c_k, 'ParameterVariableNames', {'kcat_pkc_nsfbp2'});


pkc_nsfbp_reac_12a = addreaction(pfpc_model, 'Ca2PKC_act + NSFBPp -> Ca2PKC_act_NSFBPp');
pkc_nsfbp_reac_12a_k = addkineticlaw(pkc_nsfbp_reac_12a,'MassAction');
set(pkc_nsfbp_reac_12a_k, 'ParameterVariableNames', {'kon_pkc_nsfbp2'});

pkc_nsfbp_reac_12b = addreaction(pfpc_model, 'Ca2PKC_act_NSFBPp -> Ca2PKC_act + NSFBPp');
pkc_nsfbp_reac_12b_k = addkineticlaw(pkc_nsfbp_reac_12b,'MassAction');
set(pkc_nsfbp_reac_12b_k, 'ParameterVariableNames', {'koff_pkc_nsfbp2'});

pkc_nsfbp_reac_12c = addreaction(pfpc_model, 'Ca2PKC_act_NSFBPp -> Ca2PKC_act + NSFBPpp');
pkc_nsfbp_reac_12c_k = addkineticlaw(pkc_nsfbp_reac_12c,'MassAction');
set(pkc_nsfbp_reac_12c_k, 'ParameterVariableNames', {'kcat_pkc_nsfbp2'});


pkc_nsfbp_reac_13a = addreaction(pfpc_model, 'CaPKC_AA_act + NSFBPp -> CaPKC_AA_act_NSFBPp');
pkc_nsfbp_reac_13a_k = addkineticlaw(pkc_nsfbp_reac_13a,'MassAction');
set(pkc_nsfbp_reac_13a_k, 'ParameterVariableNames', {'kon_pkc_nsfbp2'});

pkc_nsfbp_reac_13b = addreaction(pfpc_model, 'CaPKC_AA_act_NSFBPp -> CaPKC_AA_act + NSFBPp');
pkc_nsfbp_reac_13b_k = addkineticlaw(pkc_nsfbp_reac_13b,'MassAction');
set(pkc_nsfbp_reac_13b_k, 'ParameterVariableNames', {'koff_pkc_nsfbp2'});

pkc_nsfbp_reac_13c = addreaction(pfpc_model, 'CaPKC_AA_act_NSFBPp -> CaPKC_AA_act + NSFBPpp');
pkc_nsfbp_reac_13c_k = addkineticlaw(pkc_nsfbp_reac_13c,'MassAction');
set(pkc_nsfbp_reac_13c_k, 'ParameterVariableNames', {'kcat_pkc_nsfbp2'});


pkc_nsfbp_reac_14a = addreaction(pfpc_model, 'Ca2PKC_AA_act + NSFBPp -> Ca2PKC_AA_act_NSFBPp');
pkc_nsfbp_reac_14a_k = addkineticlaw(pkc_nsfbp_reac_14a,'MassAction');
set(pkc_nsfbp_reac_14a_k, 'ParameterVariableNames', {'kon_pkc_nsfbp2'});

pkc_nsfbp_reac_14b = addreaction(pfpc_model, 'Ca2PKC_AA_act_NSFBPp -> Ca2PKC_AA_act + NSFBPp');
pkc_nsfbp_reac_14b_k = addkineticlaw(pkc_nsfbp_reac_14b,'MassAction');
set(pkc_nsfbp_reac_14b_k, 'ParameterVariableNames', {'koff_pkc_nsfbp2'});

pkc_nsfbp_reac_14c = addreaction(pfpc_model, 'Ca2PKC_AA_act_NSFBPp -> Ca2PKC_AA_act + NSFBPpp');
pkc_nsfbp_reac_14c_k = addkineticlaw(pkc_nsfbp_reac_14c,'MassAction');
set(pkc_nsfbp_reac_14c_k, 'ParameterVariableNames', {'kcat_pkc_nsfbp2'});


pkc_nsfbp_reac_15a = addreaction(pfpc_model, 'DAG_CaPKC_act + NSFBPp -> DAG_CaPKC_act_NSFBPp');
pkc_nsfbp_reac_15a_k = addkineticlaw(pkc_nsfbp_reac_15a,'MassAction');
set(pkc_nsfbp_reac_15a_k, 'ParameterVariableNames', {'kon_pkc_nsfbp2'});

pkc_nsfbp_reac_15b = addreaction(pfpc_model, 'DAG_CaPKC_act_NSFBPp -> DAG_CaPKC_act + NSFBPp');
pkc_nsfbp_reac_15b_k = addkineticlaw(pkc_nsfbp_reac_15b,'MassAction');
set(pkc_nsfbp_reac_15b_k, 'ParameterVariableNames', {'koff_pkc_nsfbp2'});

pkc_nsfbp_reac_15c = addreaction(pfpc_model, 'DAG_CaPKC_act_NSFBPp -> DAG_CaPKC_act + NSFBPpp');
pkc_nsfbp_reac_15c_k = addkineticlaw(pkc_nsfbp_reac_15c,'MassAction');
set(pkc_nsfbp_reac_15c_k, 'ParameterVariableNames', {'kcat_pkc_nsfbp2'});


pkc_nsfbp_reac_16a = addreaction(pfpc_model, 'DAG_Ca2PKC_act + NSFBPp -> DAG_Ca2PKC_act_NSFBPp');
pkc_nsfbp_reac_16a_k = addkineticlaw(pkc_nsfbp_reac_16a,'MassAction');
set(pkc_nsfbp_reac_16a_k, 'ParameterVariableNames', {'kon_pkc_nsfbp2'});

pkc_nsfbp_reac_16b = addreaction(pfpc_model, 'DAG_Ca2PKC_act_NSFBPp -> DAG_Ca2PKC_act + NSFBPp');
pkc_nsfbp_reac_16b_k = addkineticlaw(pkc_nsfbp_reac_16b,'MassAction');
set(pkc_nsfbp_reac_16b_k, 'ParameterVariableNames', {'koff_pkc_nsfbp2'});

pkc_nsfbp_reac_16c = addreaction(pfpc_model, 'DAG_Ca2PKC_act_NSFBPp -> DAG_Ca2PKC_act + NSFBPpp');
pkc_nsfbp_reac_16c_k = addkineticlaw(pkc_nsfbp_reac_16c,'MassAction');
set(pkc_nsfbp_reac_16c_k, 'ParameterVariableNames', {'kcat_pkc_nsfbp2'});


pkc_nsfbp_reac_17a = addreaction(pfpc_model, 'DAG_CaPKC_AA_act + NSFBPp -> DAG_CaPKC_AA_act_NSFBPp');
pkc_nsfbp_reac_17a_k = addkineticlaw(pkc_nsfbp_reac_17a,'MassAction');
set(pkc_nsfbp_reac_17a_k, 'ParameterVariableNames', {'kon_pkc_nsfbp2'});

pkc_nsfbp_reac_17b = addreaction(pfpc_model, 'DAG_CaPKC_AA_act_NSFBPp -> DAG_CaPKC_AA_act + NSFBPp');
pkc_nsfbp_reac_17b_k = addkineticlaw(pkc_nsfbp_reac_17b,'MassAction');
set(pkc_nsfbp_reac_17b_k, 'ParameterVariableNames', {'koff_pkc_nsfbp2'});

pkc_nsfbp_reac_17c = addreaction(pfpc_model, 'DAG_CaPKC_AA_act_NSFBPp -> DAG_CaPKC_AA_act + NSFBPpp');
pkc_nsfbp_reac_17c_k = addkineticlaw(pkc_nsfbp_reac_17c,'MassAction');
set(pkc_nsfbp_reac_17c_k, 'ParameterVariableNames', {'kcat_pkc_nsfbp2'});


pkc_nsfbp_reac_18a = addreaction(pfpc_model, 'DAG_Ca2PKC_AA_act + NSFBPp -> DAG_Ca2PKC_AA_act_NSFBPp');
pkc_nsfbp_reac_18a_k = addkineticlaw(pkc_nsfbp_reac_18a,'MassAction');
set(pkc_nsfbp_reac_18a_k, 'ParameterVariableNames', {'kon_pkc_nsfbp2'});

pkc_nsfbp_reac_18b = addreaction(pfpc_model, 'DAG_Ca2PKC_AA_act_NSFBPp -> DAG_Ca2PKC_AA_act + NSFBPp');
pkc_nsfbp_reac_18b_k = addkineticlaw(pkc_nsfbp_reac_18b,'MassAction');
set(pkc_nsfbp_reac_18b_k, 'ParameterVariableNames', {'koff_pkc_nsfbp2'});

pkc_nsfbp_reac_18c = addreaction(pfpc_model, 'DAG_Ca2PKC_AA_act_NSFBPp -> DAG_Ca2PKC_AA_act + NSFBPpp');
pkc_nsfbp_reac_18c_k = addkineticlaw(pkc_nsfbp_reac_18c,'MassAction');
set(pkc_nsfbp_reac_18c_k, 'ParameterVariableNames', {'kcat_pkc_nsfbp2'});


pkc_nsfbp_reac_9a = addreaction(pfpc_model, 'NSFBPp + NSFBPpase -> NSFBPp_NSFBPpase');
pkc_nsfbp_reac_9a_k = addkineticlaw(pkc_nsfbp_reac_9a,'MassAction');
set(pkc_nsfbp_reac_9a_k, 'ParameterVariableNames', {'kon_nsfbppase_nsfbp'});

pkc_nsfbp_reac_9b = addreaction(pfpc_model, 'NSFBPp_NSFBPpase -> NSFBPp + NSFBPpase');
pkc_nsfbp_reac_9b_k = addkineticlaw(pkc_nsfbp_reac_9b,'MassAction');
set(pkc_nsfbp_reac_9b_k, 'ParameterVariableNames', {'koff_nsfbppase_nsfbp'});

pkc_nsfbp_reac_9c = addreaction(pfpc_model, 'NSFBPp_NSFBPpase -> NSFBP + NSFBPpase');
pkc_nsfbp_reac_9c_k = addkineticlaw(pkc_nsfbp_reac_9c,'MassAction');
set(pkc_nsfbp_reac_9c_k, 'ParameterVariableNames', {'kcat_nsfbppase'});


pkc_nsfbp_reac_10a = addreaction(pfpc_model, 'NSFBPpp + NSFBPpase -> NSFBPpp_NSFBPpase');
pkc_nsfbp_reac_10a_k = addkineticlaw(pkc_nsfbp_reac_10a,'MassAction');
set(pkc_nsfbp_reac_10a_k, 'ParameterVariableNames', {'kon_nsfbppase_nsfbp'});

pkc_nsfbp_reac_10b = addreaction(pfpc_model, 'NSFBPpp_NSFBPpase -> NSFBPpp + NSFBPpase');
pkc_nsfbp_reac_10b_k = addkineticlaw(pkc_nsfbp_reac_10b,'MassAction');
set(pkc_nsfbp_reac_10b_k, 'ParameterVariableNames', {'koff_nsfbppase_nsfbp'});

pkc_nsfbp_reac_10c = addreaction(pfpc_model, 'NSFBPpp_NSFBPpase -> NSFBPp + NSFBPpase');
pkc_nsfbp_reac_10c_k = addkineticlaw(pkc_nsfbp_reac_10c,'MassAction');
set(pkc_nsfbp_reac_10c_k, 'ParameterVariableNames', {'kcat_nsfbppase'});


%%%%%%%%%%%%%Binding of nsf to nsf-binding-protein
pkc_nsfbp_reac_13a = addreaction(pfpc_model, 'nsf_3no + NSFBPpp -> nsf_3no_NSFBPpp');
pkc_nsfbp_reac_13a_k = addkineticlaw(pkc_nsfbp_reac_13a,'MassAction');
set(pkc_nsfbp_reac_13a_k, 'ParameterVariableNames', {'kon_nsf_nsfbp'});

pkc_nsfbp_reac_13b = addreaction(pfpc_model, 'nsf_3no_NSFBPpp -> nsf_3no + NSFBPpp');
pkc_nsfbp_reac_13b_k = addkineticlaw(pkc_nsfbp_reac_13b,'MassAction');
set(pkc_nsfbp_reac_13b_k, 'ParameterVariableNames', {'koff_nsf_nsfbp'});

%%%%PKA Activation Model%%%%
%Species
PKA = addspecies(cytosol, 'PKA', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
R2C2 = addspecies(cytosol, 'R2C2', 'InitialAmount', 0.07, 'InitialAmountUnits', 'micromolarity');
R2C2_cAMP = addspecies(cytosol, 'R2C2_cAMP', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
R2C2_2cAMPbb = addspecies(cytosol, 'R2C2_2cAMPbb', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
R2C2_2cAMPab = addspecies(cytosol, 'R2C2_2cAMPab', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
R2C2_3cAMP = addspecies(cytosol, 'R2C2_3cAMP', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
R2C2_4cAMP = addspecies(cytosol, 'R2C2_4cAMP', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
R2C_2cAMPab = addspecies(cytosol, 'R2C_2cAMPab', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
R2C_3cAMP = addspecies(cytosol, 'R2C_3cAMP', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
R2C_4cAMP = addspecies(cytosol, 'R2C_4cAMP', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
R2_4cAMP = addspecies(cytosol, 'R2_4cAMP', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');

%Parameters
kon_r2c2_camp1 = addparameter(pfpc_model,'kon_r2c2_camp1', 'Value', 2, 'ValueUnits', '1/(micromolarity*second)');
koff_r2c2_camp1 = addparameter(pfpc_model,'koff_r2c2_camp1', 'Value', 0.75, 'ValueUnits', '1/second');

kon_r2c2_camp2bb = addparameter(pfpc_model,'kon_r2c2_camp2bb', 'Value', 1, 'ValueUnits', '1/(micromolarity*second)');
koff_r2c2_camp2bb = addparameter(pfpc_model,'koff_r2c2_camp2bb', 'Value', 1.5, 'ValueUnits', '1/second');

kon_r2c2_camp2ab = addparameter(pfpc_model,'kon_r2c2_camp2ab', 'Value', 10, 'ValueUnits', '1/(micromolarity*second)');
koff_r2c2_camp2ab = addparameter(pfpc_model,'koff_r2c2_camp2ab', 'Value', 7.5, 'ValueUnits', '1/second');

kon_r2c2_camp3bb = addparameter(pfpc_model,'kon_r2c2_camp3bb', 'Value', 20, 'ValueUnits', '1/(micromolarity*second)');
koff_r2c2_camp3bb = addparameter(pfpc_model,'koff_r2c2_camp3bb', 'Value', 7.5, 'ValueUnits', '1/second');

kon_r2c2_camp3ab = addparameter(pfpc_model,'kon_r2c2_camp3ab', 'Value', 1, 'ValueUnits', '1/(micromolarity*second)');
koff_r2c2_camp3ab = addparameter(pfpc_model,'koff_r2c2_camp3ab', 'Value', 0.75, 'ValueUnits', '1/second');

kon_r2c2_camp4 = addparameter(pfpc_model,'kon_r2c2_camp4', 'Value', 10, 'ValueUnits', '1/(micromolarity*second)');
koff_r2c2_camp4 = addparameter(pfpc_model,'koff_r2c2_camp4', 'Value', 15, 'ValueUnits', '1/second');

kact_pka_23 = addparameter(pfpc_model,'kact_pka_23', 'Value', 0.005, 'ValueUnits', '1/second');
kdeact_pka_23 = addparameter(pfpc_model,'kdeact_pka_23', 'Value', 5, 'ValueUnits', '1/(micromolarity*second)');

kact_pka_4 = addparameter(pfpc_model,'kact_pka_4', 'Value', 6, 'ValueUnits', '1/second');
kdeact_pka_4 = addparameter(pfpc_model,'kdeact_pka_4', 'Value', 5, 'ValueUnits', '1/(micromolarity*second)');

kon_r2c_camp3 = addparameter(pfpc_model,'kon_r2c_camp3', 'Value', 1, 'ValueUnits', '1/(micromolarity*second)');
koff_r2c_camp3 = addparameter(pfpc_model,'koff_r2c_camp3', 'Value', 0.75, 'ValueUnits', '1/second');

kon_r2c_camp4 = addparameter(pfpc_model,'kon_r2c_camp4', 'Value', 10, 'ValueUnits', '1/(micromolarity*second)');
koff_r2c_camp4 = addparameter(pfpc_model,'koff_r2c_camp4', 'Value', 7.5, 'ValueUnits', '1/second');

kact_pka_r2c = addparameter(pfpc_model,'kact_pka_r2c', 'Value', 3, 'ValueUnits', '1/second');
kdeact_pka_r2c = addparameter(pfpc_model,'kdeact_pka_r2c', 'Value', 10, 'ValueUnits', '1/(micromolarity*second)');


%Reactions
pka_reac_1 = addreaction(pfpc_model, 'R2C2 + cAMP -> R2C2_cAMP');
pka_reac_1_k = addkineticlaw(pka_reac_1,'MassAction');
set(pka_reac_1_k, 'ParameterVariableNames', {'kon_r2c2_camp1'});

pka_reac_2 = addreaction(pfpc_model, 'R2C2_cAMP -> R2C2 + cAMP');
pka_reac_2_k = addkineticlaw(pka_reac_2,'MassAction');
set(pka_reac_2_k, 'ParameterVariableNames', {'koff_r2c2_camp1'});

pka_reac_3 = addreaction(pfpc_model, 'R2C2_cAMP + cAMP -> R2C2_2cAMPbb');
pka_reac_3_k = addkineticlaw(pka_reac_3,'MassAction');
set(pka_reac_3_k, 'ParameterVariableNames', {'kon_r2c2_camp2bb'});

pka_reac_4 = addreaction(pfpc_model, 'R2C2_2cAMPbb -> R2C2_cAMP + cAMP');
pka_reac_4_k = addkineticlaw(pka_reac_4,'MassAction');
set(pka_reac_4_k, 'ParameterVariableNames', {'koff_r2c2_camp2bb'});

pka_reac_5 = addreaction(pfpc_model, 'R2C2_cAMP + cAMP -> R2C2_2cAMPab');
pka_reac_5_k = addkineticlaw(pka_reac_5,'MassAction');
set(pka_reac_5_k, 'ParameterVariableNames', {'kon_r2c2_camp2ab'});

pka_reac_6 = addreaction(pfpc_model, 'R2C2_2cAMPab -> R2C2_cAMP + cAMP');
pka_reac_6_k = addkineticlaw(pka_reac_6,'MassAction');
set(pka_reac_6_k, 'ParameterVariableNames', {'koff_r2c2_camp2ab'});

pka_reac_7 = addreaction(pfpc_model, 'R2C2_2cAMPbb + cAMP -> R2C2_3cAMP');
pka_reac_7_k = addkineticlaw(pka_reac_7,'MassAction');
set(pka_reac_7_k, 'ParameterVariableNames', {'kon_r2c2_camp3bb'});

pka_reac_8 = addreaction(pfpc_model, 'R2C2_3cAMP -> R2C2_2cAMPbb + cAMP');
pka_reac_8_k = addkineticlaw(pka_reac_8,'MassAction');
set(pka_reac_8_k, 'ParameterVariableNames', {'koff_r2c2_camp3bb'});

pka_reac_9 = addreaction(pfpc_model, 'R2C2_2cAMPab + cAMP -> R2C2_3cAMP');
pka_reac_9_k = addkineticlaw(pka_reac_9,'MassAction');
set(pka_reac_9_k, 'ParameterVariableNames', {'kon_r2c2_camp3ab'});

pka_reac_10 = addreaction(pfpc_model, 'R2C2_3cAMP -> R2C2_2cAMPab');
pka_reac_10_k = addkineticlaw(pka_reac_10,'MassAction');
set(pka_reac_10_k, 'ParameterVariableNames', {'koff_r2c2_camp3ab'});

pka_reac_11 = addreaction(pfpc_model, 'R2C2_3cAMP + cAMP -> R2C2_4cAMP');
pka_reac_11_k = addkineticlaw(pka_reac_11,'MassAction');
set(pka_reac_11_k, 'ParameterVariableNames', {'kon_r2c2_camp4'});

pka_reac_12 = addreaction(pfpc_model, 'R2C2_4cAMP -> R2C2_3cAMP + cAMP');
pka_reac_12_k = addkineticlaw(pka_reac_12,'MassAction');
set(pka_reac_12_k, 'ParameterVariableNames', {'koff_r2c2_camp4'});

pka_reac_13 = addreaction(pfpc_model, 'R2C2_2cAMPab -> R2C_2cAMPab + PKA');
pka_reac_13_k = addkineticlaw(pka_reac_13,'MassAction');
set(pka_reac_13_k, 'ParameterVariableNames', {'kact_pka_23'});

pka_reac_14 = addreaction(pfpc_model, 'R2C_2cAMPab + PKA -> R2C2_2cAMPab');
pka_reac_14_k = addkineticlaw(pka_reac_14,'MassAction');
set(pka_reac_14_k, 'ParameterVariableNames', {'kdeact_pka_23'});

pka_reac_15 = addreaction(pfpc_model, 'R2C2_3cAMP -> R2C_3cAMP + PKA');
pka_reac_15_k = addkineticlaw(pka_reac_15,'MassAction');
set(pka_reac_15_k, 'ParameterVariableNames', {'kact_pka_23'});

pka_reac_16 = addreaction(pfpc_model, 'R2C_3cAMP + PKA -> R2C2_3cAMP');
pka_reac_16_k = addkineticlaw(pka_reac_16,'MassAction');
set(pka_reac_16_k, 'ParameterVariableNames', {'kdeact_pka_23'});

pka_reac_17 = addreaction(pfpc_model, 'R2C2_4cAMP -> R2C_4cAMP + PKA');
pka_reac_17_k = addkineticlaw(pka_reac_17,'MassAction');
set(pka_reac_17_k, 'ParameterVariableNames', {'kact_pka_4'});

pka_reac_18 = addreaction(pfpc_model, 'R2C_4cAMP + PKA -> R2C2_4cAMP');
pka_reac_18_k = addkineticlaw(pka_reac_18,'MassAction');
set(pka_reac_18_k, 'ParameterVariableNames', {'kdeact_pka_4'});

pka_reac_19 = addreaction(pfpc_model, 'R2C_2cAMPab + cAMP -> R2C_3cAMP');
pka_reac_19_k = addkineticlaw(pka_reac_19,'MassAction');
set(pka_reac_19_k, 'ParameterVariableNames', {'kon_r2c_camp3'});

pka_reac_20 = addreaction(pfpc_model, 'R2C_3cAMP -> R2C_2cAMPab + cAMP');
pka_reac_20_k = addkineticlaw(pka_reac_20,'MassAction');
set(pka_reac_20_k, 'ParameterVariableNames', {'koff_r2c_camp3'});

pka_reac_21 = addreaction(pfpc_model, 'R2C_3cAMP + cAMP -> R2C_4cAMP');
pka_reac_21_k = addkineticlaw(pka_reac_21,'MassAction');
set(pka_reac_21_k, 'ParameterVariableNames', {'kon_r2c_camp4'});

pka_reac_22 = addreaction(pfpc_model, 'R2C_4cAMP -> R2C_3cAMP + cAMP');
pka_reac_22_k = addkineticlaw(pka_reac_22,'MassAction');
set(pka_reac_22_k, 'ParameterVariableNames', {'koff_r2c_camp4'});

pka_reac_23 = addreaction(pfpc_model, 'R2C_4cAMP -> R2_4cAMP + PKA');
pka_reac_23_k = addkineticlaw(pka_reac_23,'MassAction');
set(pka_reac_23_k, 'ParameterVariableNames', {'kact_pka_r2c'});

pka_reac_24 = addreaction(pfpc_model, 'R2_4cAMP + PKA -> R2C_4cAMP');
pka_reac_24_k = addkineticlaw(pka_reac_24,'MassAction');
set(pka_reac_24_k, 'ParameterVariableNames', {'kdeact_pka_r2c'});

%%%%%%DARPP MODEL%%%%%

%Species
D32 = addspecies(cytosol, 'D32', 'InitialAmount', 1.8, 'InitialAmountUnits', 'micromolarity');
D32p = addspecies(cytosol, 'D32p', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
PKA_D32 = addspecies(cytosol, 'PKA_D32', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
CaN_D32p = addspecies(cytosol, 'CaN_D32p', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
D32p_PP1 = addspecies(cytosol, 'D32p_PP1', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');


%Parameters
kon_pka_d32 = addparameter(pfpc_model,'kon_pka_d32', 'Value', 5.6, 'ValueUnits', '1/(micromolarity*second)');
koff_pka_d32 = addparameter(pfpc_model,'koff_pka_d32', 'Value', 10.8, 'ValueUnits', '1/second');
kcat_pka = addparameter(pfpc_model,'kcat_pka', 'Value', 2.7, 'ValueUnits', '1/second');

kon_can_d32 = addparameter(pfpc_model,'kon_can_d32', 'Value', 6.3, 'ValueUnits', '1/(micromolarity*second)', 'ConstantValue', false);
koff_can_d32 = addparameter(pfpc_model,'koff_can_d32', 'Value', 0.8, 'ValueUnits', '1/second');
kcat_can_d32 = addparameter(pfpc_model,'kcat_can_d32', 'Value', 2, 'ValueUnits', '1/second');

kon_d32_pp1 = addparameter(pfpc_model,'kon_d32_pp1', 'Value', 2, 'ValueUnits', '1/(micromolarity*second)');
koff_d32_pp1 = addparameter(pfpc_model,'koff_d32_pp1', 'Value', 0.5, 'ValueUnits', '1/second');

%Reactions

pka_reac_25 = addreaction(pfpc_model, 'PKA + D32 -> PKA_D32');
pka_reac_25_k = addkineticlaw(pka_reac_25,'MassAction');
set(pka_reac_25_k, 'ParameterVariableNames', {'kon_pka_d32'});

pka_reac_26 = addreaction(pfpc_model, 'PKA_D32 -> PKA + D32');
pka_reac_26_k = addkineticlaw(pka_reac_26,'MassAction');
set(pka_reac_26_k, 'ParameterVariableNames', {'koff_pka_d32'});

pka_reac_27 = addreaction(pfpc_model, 'PKA_D32 -> PKA + D32p');
pka_reac_27_k = addkineticlaw(pka_reac_27,'MassAction');
set(pka_reac_27_k, 'ParameterVariableNames', {'kcat_pka'});

pka_reac_28 = addreaction(pfpc_model, 'CaN_CaM + D32p -> CaN_D32p');
pka_reac_28_k = addkineticlaw(pka_reac_28,'MassAction');
set(pka_reac_28_k, 'ParameterVariableNames', {'kon_can_d32'});

pka_reac_29 = addreaction(pfpc_model, 'CaN_D32p -> CaN_CaM + D32p');
pka_reac_29_k = addkineticlaw(pka_reac_29,'MassAction');
set(pka_reac_29_k, 'ParameterVariableNames', {'koff_can_d32'});

pka_reac_30 = addreaction(pfpc_model, 'CaN_D32p -> CaN_CaM + D32');
pka_reac_30_k = addkineticlaw(pka_reac_30,'MassAction');
set(pka_reac_30_k, 'ParameterVariableNames', {'kcat_can_d32'});

pka_reac_31 = addreaction(pfpc_model, 'D32p + PP1 -> D32p_PP1');
pka_reac_31_k = addkineticlaw(pka_reac_31,'MassAction');
set(pka_reac_31_k, 'ParameterVariableNames', {'kon_d32_pp1'});

pka_reac_32 = addreaction(pfpc_model, 'D32p_PP1 -> D32p + PP1');
pka_reac_32_k = addkineticlaw(pka_reac_32,'MassAction');
set(pka_reac_32_k, 'ParameterVariableNames', {'koff_d32_pp1'});

kdeg_no2 = addparameter(pfpc_model,'kdeg_no2', 'Value', 0.1, 'ValueUnits', '1/second');
no_deg = addreaction(pfpc_model, 'NO -> null');
no_deg_k = addkineticlaw(no_deg,'MassAction');
set(no_deg_k, 'ParameterVariableNames', {'kdeg_no2'});

%%%%Repeated Assignments
ERK_act_Rule = addrule(pfpc_model,'ERK_act = ppERK','repeatedAssignment');
Ras_act_Rule = addrule(pfpc_model,'Ras_act = ras_gtp + ras_RAF + ras_RAF_MEK + ras_RAF_pMEK','repeatedAssignment');
RAF_act_Rule = addrule(pfpc_model,'RAF_act = ras_RAF_act+ras_RAF_MEK+ras_RAF_pMEK+ras_RAF_RAFK+ras_pRAF_act','repeatedAssignment');
MEK_phospho_Rule = addrule(pfpc_model,'MEK_phospho = pMEK + ras_RAF_pMEK + ppMEK + ppMEK_ERK + ppMEK_pERK + ppMEK_PP2A + pMEK_PP2A','repeatedAssignment');
PLA2_active_Rule = addrule(pfpc_model,'PLA2_active = Ca2pPLA2_act + DAG_Ca2pPLA2_act + Ca2PLA2_act + DAG_Ca2PLA2_act + Ca2pPLA2_act_APC + DAG_Ca2pPLA2_act_APC + Ca2PLA2_act_APC + DAG_Ca2PLA2_act_APC','repeatedAssignment');
PKC_active_Rule = addrule(pfpc_model,'PKC_active = CaPKC_act + Ca2PKC_act + CaPKC_AA_act + Ca2PKC_AA_act + DAG_CaPKC_act + DAG_Ca2PKC_act + DAG_CaPKC_AA_act + DAG_Ca2PKC_AA_act','repeatedAssignment');
AMPAR_PSD_Rule = addrule(pfpc_model,'ampar_psd = amparp_pick_nsf + ampar_pick_nsf + amparp_pick_nsf_no + amparp_pick_nsf_2no + amparp_pick_nsf_3no + ampar_pick_nsf_no + ampar_pick_nsf_2no + ampar_pick_nsf_3no +ampar + ampar_x + amparp + amparp_x + ampar_pick + amparp_pick + amparp_pp2a + amparp_pick_pp2a + ampar_CaPKC_act + ampar_x_CaPKC_act + ampar_pick_CaPKC_act + ampar_Ca2PKC_act + ampar_x_Ca2PKC_act + ampar_pick_Ca2PKC_act + ampar_CaPKC_AA_act + ampar_x_CaPKC_AA_act+ ampar_pick_CaPKC_AA_act + ampar_Ca2PKC_AA_act + ampar_x_Ca2PKC_AA_act + ampar_pick_Ca2PKC_AA_act + ampar_DAG_CaPKC_act + ampar_x_DAG_CaPKC_act + ampar_pick_DAG_CaPKC_act + ampar_DAG_Ca2PKC_act + ampar_x_DAG_Ca2PKC_act + ampar_pick_DAG_Ca2PKC_act + ampar_DAG_CaPKC_AA_act + ampar_x_DAG_CaPKC_AA_act + ampar_pick_DAG_CaPKC_AA_act + ampar_DAG_Ca2PKC_AA_act + ampar_x_DAG_Ca2PKC_AA_act + ampar_pick_DAG_Ca2PKC_AA_act','repeatedAssignment');
AMPAR_cytosol_Rule = addrule(pfpc_model,'ampar_cytosol = ampar_c + amparp_c + ampar_c_pick + amparp_c_pick','repeatedAssignment');

NSF_NO3_Rule = addrule(pfpc_model,'NSF_NO3_tot = nsf_3no + nsf_3no_NSFBPpp','repeatedAssignment');

NSFBPp_Rule = addrule(pfpc_model,'NSFBPp_tot = NSFBPp + CaPKC_act_NSFBPp + Ca2PKC_act_NSFBPp + CaPKC_AA_act_NSFBPp + Ca2PKC_AA_act_NSFBPp + DAG_CaPKC_act_NSFBPp + DAG_Ca2PKC_act_NSFBPp + DAG_CaPKC_AA_act_NSFBPp + DAG_Ca2PKC_AA_act_NSFBPp + NSFBPp_NSFBPpase','repeatedAssignment');

NSFBPpp_Rule = addrule(pfpc_model,'NSFBPpp_tot = NSFBPpp + nsf_no_NSFBPpp + nsf_2no_NSFBPpp + nsf_3no_NSFBPpp','repeatedAssignment');




%%%%Influx parameter to vary Ca spike%%%%%
kinflux_v = 0;

%Parameters
kflux_ip3r = addparameter(pfpc_model,'kflux_ip3r', 'Value', 8, 'ValueUnits', '1/(micromolarity*second)');%450

%Reactions
%E1
ca_reac_1 = addreaction(pfpc_model, 'Ca_ER + IP3R_open -> Ca + IP3R_open + Ca_ER')
ca_reac_1_k = addkineticlaw(ca_reac_1,'MassAction');
set(ca_reac_1_k, 'ParameterVariableNames', {'kflux_ip3r'});
ca_reac_2 = addreaction(pfpc_model, 'Ca + IP3R_open -> IP3R_open')
ca_reac_2_k = addkineticlaw(ca_reac_2,'MassAction');
set(ca_reac_2_k, 'ParameterVariableNames', {'kflux_ip3r'});

%%%%Influx parameter to vary Ca spike%%%%%
%Ca influx
ca_reac_31 = addreaction(pfpc_model, 'null -> Ca');
ca_reac_31_k = addkineticlaw(ca_reac_31,'MassAction');
set(ca_reac_31_k, 'ParameterVariableNames', {'kinflux'});

%Ca influx (PF AMAPR)
ca_reac_pf = addreaction(pfpc_model, 'null -> Ca');
ca_reac_pf_k = addkineticlaw(ca_reac_pf,'MassAction');
set(ca_reac_pf_k, 'ParameterVariableNames', {'kinflux_pf'});

%Glutamate pulses (PF)
glu_pulse_1 = addreaction(pfpc_model, 'null -> glu');
glu_pulse_1_k = addkineticlaw(glu_pulse_1,'MassAction');
set(glu_pulse_1_k, 'ParameterVariableNames', {'kinflux_glu'});


%%%%%%%% AG2 Model %%%%%%%%%%
%%%%%%%% AG2 Model %%%%%%%%%%
% Add species
AG2 = addspecies(cytosol, 'AG2', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
DAGL = addspecies(cytosol, 'DAGL', 'InitialAmount', 0.5, 'InitialAmountUnits', 'micromolarity');
DAGL_DAG = addspecies(cytosol, 'DAGL_DAG', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
MGAT = addspecies(cytosol, 'MGAT', 'InitialAmount', 1, 'InitialAmountUnits', 'micromolarity');
MGAT_AG2 = addspecies(cytosol, 'MGAT_AG2', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');

AG2_buff = addspecies(cytosol, 'AG2_buff', 'InitialAmount', 5, 'InitialAmountUnits', 'micromolarity');
AG2_AG2_buff = addspecies(cytosol, 'AG2_AG2_buff', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');

AG2T = addspecies(cytosol, 'AG2T', 'InitialAmount', 0.1, 'InitialAmountUnits', 'micromolarity');
AG2T_AG2 = addspecies(cytosol, 'AG2T_AG2', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
AG2T_AG2_2 = addspecies(cytosol, 'AG2T_AG2_2', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
AG2_ex = addspecies(cytosol, 'AG2_ex', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
AG2_deg = addspecies(cytosol, 'AG2_deg', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');

MAGL = addspecies(cytosol, 'MAGL', 'InitialAmount', 0.5, 'InitialAmountUnits', 'micromolarity');
MAGL_AG2 = addspecies(cytosol, 'MAGL_AG2', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
GAT = addspecies(cytosol, 'GAT', 'InitialAmount', 0.5, 'InitialAmountUnits', 'micromolarity');
GAT_AA = addspecies(cytosol, 'GAT_AA', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');


% Add parameters
kon_dagl_dag = addparameter(pfpc_model,'kon_dagl_dag', 'Value', 1, 'ValueUnits', '1/(micromolarity*second)');
koff_dagl_dag = addparameter(pfpc_model,'koff_dagl_dag', 'Value', 5, 'ValueUnits', '1/second');
kcat_dagl = addparameter(pfpc_model,'kcat_dagl', 'Value', 1, 'ValueUnits', '1/second');

kon_mgat_ag2 = addparameter(pfpc_model,'kon_mgat_ag2', 'Value', 1, 'ValueUnits', '1/(micromolarity*second)');
koff_mgat_ag2 = addparameter(pfpc_model,'koff_mgat_ag2', 'Value', 5, 'ValueUnits', '1/second');
kcat_mgat = addparameter(pfpc_model,'kcat_mgat', 'Value', 0.5, 'ValueUnits', '1/second');

kon_ag2_buff = addparameter(pfpc_model,'kon_ag2_buff', 'Value', 100, 'ValueUnits', '1/(micromolarity*second)');
koff_ag2_buff = addparameter(pfpc_model,'koff_ag2_buff', 'Value', 1, 'ValueUnits', '1/second');

kon_ag2t_ag2 = addparameter(pfpc_model,'kon_ag2t_ag2', 'Value', 1, 'ValueUnits', '1/(micromolarity*second)');
kon_ag2t_ag2_2 = addparameter(pfpc_model,'kon_ag2t_ag2_2', 'Value', 20, 'ValueUnits', '1/(micromolarity*second)');

koff_ag2t_ag2 = addparameter(pfpc_model,'koff_ag2t_ag2', 'Value', 10, 'ValueUnits', '1/second');
kflux_ag2t = addparameter(pfpc_model,'kflux_ag2t', 'Value', 1, 'ValueUnits', '1/second');

kdiff_ag2 = addparameter(pfpc_model,'kdiff_ag2', 'Value', 0.05, 'ValueUnits', '1/second');

kon_magl_ag2 = addparameter(pfpc_model,'kon_magl_ag2', 'Value', 1, 'ValueUnits', '1/(micromolarity*second)');
koff_magl_ag2 = addparameter(pfpc_model,'koff_magl_ag2', 'Value', 5, 'ValueUnits', '1/second');
kcat_magl = addparameter(pfpc_model,'kcat_magl', 'Value', 1, 'ValueUnits', '1/second');

kon_gat_aa = addparameter(pfpc_model,'kon_gat_aa', 'Value', 1, 'ValueUnits', '1/(micromolarity*second)');
koff_gat_aa = addparameter(pfpc_model,'koff_gat_aa', 'Value', 3, 'ValueUnits', '1/second');
kcat_gat = addparameter(pfpc_model,'kcat_gat', 'Value', 10, 'ValueUnits', '1/second');

% Add reactions
ag2_reac_1 = addreaction(pfpc_model, 'DAGL + DAG -> DAGL_DAG')
ag2_reac_1_k = addkineticlaw(ag2_reac_1,'MassAction');
set(ag2_reac_1_k, 'ParameterVariableNames', {'kon_dagl_dag'});

ag2_reac_2 = addreaction(pfpc_model, 'DAGL_DAG -> DAGL + DAG')
ag2_reac_2_k = addkineticlaw(ag2_reac_2,'MassAction');
set(ag2_reac_2_k, 'ParameterVariableNames', {'koff_dagl_dag'});

ag2_reac_3 = addreaction(pfpc_model, 'DAGL_DAG -> DAGL + AG2')
ag2_reac_3_k = addkineticlaw(ag2_reac_3,'MassAction');
set(ag2_reac_3_k, 'ParameterVariableNames', {'kcat_dagl'});


ag2_reac_1a = addreaction(pfpc_model, 'MAGL + AG2 -> MAGL_AG2')
ag2_reac_1a_k = addkineticlaw(ag2_reac_1a,'MassAction');
set(ag2_reac_1a_k, 'ParameterVariableNames', {'kon_magl_ag2'});

ag2_reac_2a = addreaction(pfpc_model, 'MAGL_AG2 -> MAGL + AG2')
ag2_reac_2a_k = addkineticlaw(ag2_reac_2a,'MassAction');
set(ag2_reac_2a_k, 'ParameterVariableNames', {'koff_magl_ag2'});

ag2_reac_3a = addreaction(pfpc_model, 'MAGL_AG2 -> MAGL + AA')
ag2_reac_3a_k = addkineticlaw(ag2_reac_3a,'MassAction');
set(ag2_reac_3a_k, 'ParameterVariableNames', {'kcat_magl'});


ag2_reac_4 = addreaction(pfpc_model, 'MGAT + AG2 -> MGAT_AG2')
ag2_reac_4_k = addkineticlaw(ag2_reac_4,'MassAction');
set(ag2_reac_4_k, 'ParameterVariableNames', {'kon_mgat_ag2'});

ag2_reac_5 = addreaction(pfpc_model, 'MGAT_AG2 -> MGAT + AG2')
ag2_reac_5_k = addkineticlaw(ag2_reac_5,'MassAction');
set(ag2_reac_5_k, 'ParameterVariableNames', {'koff_mgat_ag2'});

ag2_reac_6 = addreaction(pfpc_model, 'MGAT_AG2 -> MGAT + DAG')
ag2_reac_6_k = addkineticlaw(ag2_reac_6,'MassAction');
set(ag2_reac_6_k, 'ParameterVariableNames', {'kcat_mgat'});


ag2_reac_4a = addreaction(pfpc_model, 'GAT + AA -> GAT_AA')
ag2_reac_4a_k = addkineticlaw(ag2_reac_4a,'MassAction');
set(ag2_reac_4a_k, 'ParameterVariableNames', {'kon_gat_aa'});

ag2_reac_5a = addreaction(pfpc_model, 'GAT_AA -> GAT + AA')
ag2_reac_5a_k = addkineticlaw(ag2_reac_5a,'MassAction');
set(ag2_reac_5a_k, 'ParameterVariableNames', {'koff_gat_aa'});

ag2_reac_6a = addreaction(pfpc_model, 'GAT_AA -> GAT + AG2')
ag2_reac_6a_k = addkineticlaw(ag2_reac_6a,'MassAction');
set(ag2_reac_6a_k, 'ParameterVariableNames', {'kcat_gat'});


ag2_reac_7 = addreaction(pfpc_model, 'AG2_buff + AG2 -> AG2_AG2_buff')
ag2_reac_7_k = addkineticlaw(ag2_reac_7,'MassAction');
set(ag2_reac_7_k, 'ParameterVariableNames', {'kon_ag2_buff'});

ag2_reac_8 = addreaction(pfpc_model, '	AG2_AG2_buff -> AG2_buff + AG2')
ag2_reac_8_k = addkineticlaw(ag2_reac_8,'MassAction');
set(ag2_reac_8_k, 'ParameterVariableNames', {'koff_ag2_buff'});

%AG2 transport
ag2_reac_9 = addreaction(pfpc_model, 'AG2 + AG2T -> AG2T_AG2')
ag2_reac_9_k = addkineticlaw(ag2_reac_9,'MassAction');
set(ag2_reac_9_k, 'ParameterVariableNames', {'kon_ag2t_ag2'});

ag2_reac_10 = addreaction(pfpc_model, 'AG2T_AG2 -> AG2 + AG2T')
ag2_reac_10_k = addkineticlaw(ag2_reac_10,'MassAction');
set(ag2_reac_10_k, 'ParameterVariableNames', {'koff_ag2t_ag2'});

ag2_reac_9a = addreaction(pfpc_model, 'AG2 + AG2T_AG2 -> AG2T_AG2_2')
ag2_reac_9a_k = addkineticlaw(ag2_reac_9a,'MassAction');
set(ag2_reac_9a_k, 'ParameterVariableNames', {'kon_ag2t_ag2_2'});

ag2_reac_10a = addreaction(pfpc_model, 'AG2T_AG2_2 -> AG2 + AG2T_AG2')
ag2_reac_10a_k = addkineticlaw(ag2_reac_10a,'MassAction');
set(ag2_reac_10a_k, 'ParameterVariableNames', {'koff_ag2t_ag2'});

ag2_reac_11 = addreaction(pfpc_model, 'AG2T_AG2_2 -> AG2_ex + AG2_ex + AG2T')
ag2_reac_11_k = addkineticlaw(ag2_reac_11,'MassAction');
set(ag2_reac_11_k, 'ParameterVariableNames', {'kflux_ag2t'});

ag2_reac_12 = addreaction(pfpc_model, 'AG2_ex -> AG2_deg')
ag2_reac_12_k = addkineticlaw(ag2_reac_12,'MassAction');
set(ag2_reac_12_k, 'ParameterVariableNames', {'kdiff_ag2'});

%%%CB1R Model%%%
%%%%%CB1R Model%%%%%
%Species (uM)
AG2_ex_deg = addspecies(cytosol, 'AG2_ex_deg', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
CB1R = addspecies(cytosol, 'CB1R', 'InitialAmount', 8.3333, 'InitialAmountUnits', 'micromolarity');
CB1R_gq = addspecies(cytosol, 'CB1R_gq', 'InitialAmount', 6.6667, 'InitialAmountUnits', 'micromolarity');
CB1R_AG2_ex = addspecies(cytosol, 'CB1R_AG2_ex', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
CB1R_gq_AG2_ex = addspecies(cytosol, 'CB1R_gq_AG2_ex', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
gq_gdp_p = addspecies(cytosol, 'gq_gdp_p', 'InitialAmount', 43.333, 'InitialAmountUnits', 'micromolarity');
gi_gdp_p = addspecies(cytosol, 'gi_gdp_p', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
gi_gtp_p = addspecies(cytosol, 'gi_gtp_p', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
gbc_p = addspecies(cytosol, 'gbc_p', 'InitialAmount', 0.0, 'InitialAmountUnits', 'micromolarity');
CB1R_tot = addspecies(cytosol, 'CB1R_tot', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
CB1R_tot_Rule = addrule(pfpc_model,'CB1R_tot = CB1R + CB1R_gq + CB1R_AG2_ex + CB1R_gq_AG2_ex','repeatedAssignment');

%Parameters
kon_cb1r_ag2 = addparameter(pfpc_model,'kon_cb1r_ag2', 'Value', 0.02, 'ValueUnits', '1/(micromolarity*second)');
koff_cb1r_ag2 = addparameter(pfpc_model,'koff_cb1r_ag2', 'Value', 100, 'ValueUnits', '1/second');
kon_cb1r_gq = addparameter(pfpc_model,'kon_cb1r_gq', 'Value', 2, 'ValueUnits', '1/(micromolarity*second)');
koff_cb1r_gq = addparameter(pfpc_model,'koff_cb1r_gq', 'Value', 100, 'ValueUnits', '1/second');
kact_gi = addparameter(pfpc_model,'kact_gi', 'Value', 40, 'ValueUnits', '1/second');
kdeact_gi = addparameter(pfpc_model,'kdeact_gi', 'Value', 0.04, 'ValueUnits', '1/second');
kdeg_ag2 = addparameter(pfpc_model,'kdeg_ag2', 'Value', 0.02, 'ValueUnits', '1/second', 'ConstantValue', false);%9000

%Reactions
%A1
cb1r_reac_1 = addreaction(pfpc_model, 'CB1R + AG2_ex -> CB1R_AG2_ex')
cb1r_reac_1_k = addkineticlaw(cb1r_reac_1,'MassAction');
set(cb1r_reac_1_k, 'ParameterVariableNames', {'kon_cb1r_ag2'});

cb1r_reac_2 = addreaction(pfpc_model, 'CB1R_AG2_ex -> CB1R + AG2_ex')
cb1r_reac_2_k = addkineticlaw(cb1r_reac_2,'MassAction');
set(cb1r_reac_2_k, 'ParameterVariableNames', {'koff_cb1r_ag2'});
%A2
cb1r_reac_3 = addreaction(pfpc_model, 'CB1R_gq + AG2_ex -> CB1R_gq_AG2_ex')
cb1r_reac_3_k = addkineticlaw(cb1r_reac_3,'MassAction');
set(cb1r_reac_3_k, 'ParameterVariableNames', {'kon_cb1r_ag2'});

cb1r_reac_4 = addreaction(pfpc_model, 'CB1R_gq_AG2_ex -> CB1R_gq + AG2_ex')
cb1r_reac_4_k = addkineticlaw(cb1r_reac_4,'MassAction');
set(cb1r_reac_4_k, 'ParameterVariableNames', {'koff_cb1r_ag2'});

%A3
cb1r_reac_5 = addreaction(pfpc_model, 'CB1R + gq_gdp_p -> CB1R_gq')
cb1r_reac_5_k = addkineticlaw(cb1r_reac_5,'MassAction');
set(cb1r_reac_5_k, 'ParameterVariableNames', {'kon_cb1r_gq'});

cb1r_reac_6 = addreaction(pfpc_model, 'CB1R_gq -> CB1R + gq_gdp_p')
cb1r_reac_6_k = addkineticlaw(cb1r_reac_6,'MassAction');
set(cb1r_reac_6_k, 'ParameterVariableNames', {'koff_cb1r_gq'});

%A4

cb1r_reac_7 = addreaction(pfpc_model, 'CB1R_AG2_ex + gq_gdp_p -> CB1R_gq_AG2_ex')
cb1r_reac_7_k = addkineticlaw(cb1r_reac_7,'MassAction');
set(cb1r_reac_7_k, 'ParameterVariableNames', {'kon_cb1r_gq'});

cb1r_reac_8 = addreaction(pfpc_model, 'CB1R_gq_AG2_ex -> CB1R_AG2_ex + gq_gdp_p')
cb1r_reac_8_k = addkineticlaw(cb1r_reac_8,'MassAction');
set(cb1r_reac_8_k, 'ParameterVariableNames', {'koff_cb1r_gq'});

%A5
cb1r_reac_9 = addreaction(pfpc_model, 'CB1R_gq_AG2_ex -> CB1R_AG2_ex + gi_gtp_p + gbc_p')
cb1r_reac_9_k = addkineticlaw(cb1r_reac_9,'MassAction');
set(cb1r_reac_9_k, 'ParameterVariableNames', {'kact_gi'});


%A7
cb1r_reac_11 = addreaction(pfpc_model, 'gi_gtp_p -> gi_gdp_p')
cb1r_reac_11_k = addkineticlaw(cb1r_reac_11,'MassAction');
set(cb1r_reac_11_k, 'ParameterVariableNames', {'kdeact_gi'});

%A8
cb1r_reac_12 = addreaction(pfpc_model, 'gi_gdp_p + gbc_p -> gq_gdp_p')
cb1r_reac_12_k = addkineticlaw(cb1r_reac_12,'MassAction');
set(cb1r_reac_12_k, 'ParameterVariableNames', {'ktrimer'});


%A9 CLEAR 2-AG
cb1r_reac_13 = addreaction(pfpc_model, 'AG2_ex -> AG2_ex_deg')
cb1r_reac_13_k = addkineticlaw(cb1r_reac_13,'MassAction');
set(cb1r_reac_13_k, 'ParameterVariableNames', {'kdeg_ag2'});


%%%%%PI3K Model%%%%%
%Species (uM)
PI3K = addspecies(cytosol, 'PI3K', 'InitialAmount', 1, 'InitialAmountUnits', 'micromolarity');
PI3K_gbc = addspecies(cytosol, 'PI3K_gbc', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
pip2_p = addspecies(cytosol, 'pip2_p', 'InitialAmount', 4166.7, 'InitialAmountUnits', 'micromolarity');
pip3 = addspecies(cytosol, 'pip3', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
PI3K_gbc_pip2 = addspecies(cytosol, 'PI3K_gbc_pip2', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');

%Parameters
kon_pi3k_gbc = addparameter(pfpc_model,'kon_pi3k_gbc', 'Value', 32, 'ValueUnits', '1/(micromolarity*second)');
koff_pi3k_gbc = addparameter(pfpc_model,'koff_pi3k_gbc', 'Value', 20, 'ValueUnits', '1/second');
kcat_pi3k = addparameter(pfpc_model,'kcat_pi3k', 'Value', 160, 'ValueUnits', '1/second');
kon_pi3k_pip2 = addparameter(pfpc_model,'kon_pi3k_pip2', 'Value', 1, 'ValueUnits', '1/(micromolarity*second)');
koff_pi3k_pip2 = addparameter(pfpc_model,'koff_pi3k_pip2', 'Value', 170, 'ValueUnits', '1/second');
kdeg_pip3 = addparameter(pfpc_model,'kdeg_pip3', 'Value', 10, 'ValueUnits', '1/second');
%Reactions
%B3
pi3k_reac_1 = addreaction(pfpc_model, 'PI3K + gbc_p -> PI3K_gbc')
pi3k_reac_1_k = addkineticlaw(pi3k_reac_1,'MassAction');
set(pi3k_reac_1_k, 'ParameterVariableNames', {'kon_pi3k_gbc'});

pi3k_reac_2 = addreaction(pfpc_model, 'PI3K_gbc -> PI3K + gbc_p')
pi3k_reac_2_k = addkineticlaw(pi3k_reac_2,'MassAction');
set(pi3k_reac_2_k, 'ParameterVariableNames', {'koff_pi3k_gbc'});


%B8
pi3k_reac_3 = addreaction(pfpc_model, 'PI3K_gbc + pip2_p -> PI3K_gbc_pip2')
pi3k_reac_3_k = addkineticlaw(pi3k_reac_3,'MassAction');
set(pi3k_reac_3_k, 'ParameterVariableNames', {'kon_pi3k_pip2'});

pi3k_reac_4 = addreaction(pfpc_model, 'PI3K_gbc_pip2 -> PI3K_gbc + pip2_p')
pi3k_reac_4_k = addkineticlaw(pi3k_reac_4,'MassAction');
set(pi3k_reac_4_k, 'ParameterVariableNames', {'koff_pi3k_pip2'});

%B7
pi3k_reac_5 = addreaction(pfpc_model, 'PI3K_gbc_pip2 -> PI3K + gbc_p + pip3')
pi3k_reac_5_k = addkineticlaw(pi3k_reac_5,'MassAction');
set(pi3k_reac_5_k, 'ParameterVariableNames', {'kcat_pi3k'});

%B8
pi3k_reac_6 = addreaction(pfpc_model, 'pip3 -> pip2_p')
pi3k_reac_6_k = addkineticlaw(pi3k_reac_6,'MassAction');
set(pi3k_reac_6_k, 'ParameterVariableNames', {'kdeg_pip3'});

%%%%%%%%PDK1 Model%%%%%%%%%%
%PDK1 translocation

PDK1 = addspecies(cytosol, 'PDK1', 'InitialAmount', 1, 'InitialAmountUnits', 'micromolarity');
AKT = addspecies(cytosol, 'AKT', 'InitialAmount', 1, 'InitialAmountUnits', 'micromolarity');
AKT_PDK1 = addspecies(cytosol, 'AKT_PDK1', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
AKTp_PDK1 = addspecies(cytosol, 'AKTp_PDK1', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
PDK1_PIP3 = addspecies(cytosol, 'PDK1_PIP3', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
PP2A_p = addspecies(cytosol, 'PP2A_p', 'InitialAmount', 0.5, 'InitialAmountUnits', 'micromolarity');
AKTp = addspecies(cytosol, 'AKTp', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
AKTp2 = addspecies(cytosol, 'AKTp2', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
AKTp_PP2A = addspecies(cytosol, 'AKTp_PP2A', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
AKTp2_PP2A = addspecies(cytosol, 'AKTp2_PP2A', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');


kon_pdk_pip3 = addparameter(pfpc_model,'kon_pdk_pip3', 'Value', 1.4, 'ValueUnits', '1/(micromolarity*second)');
koff_pdk_pip3 = addparameter(pfpc_model,'koff_pdk_pip3', 'Value', 9.8, 'ValueUnits', '1/second');

%Km = 0.4
kon_akt_pdk1 = addparameter(pfpc_model,'kon_akt_pdk1', 'Value', 1, 'ValueUnits', '1/(micromolarity*second)');
koff_akt_pdk1 = addparameter(pfpc_model,'koff_akt_pdk1', 'Value', 40, 'ValueUnits', '1/second');
kcat_pdk1 = addparameter(pfpc_model,'kcat_pdk1', 'Value', 1, 'ValueUnits', '1/second');

%Km = 4.8
kon_akt_pp2a = addparameter(pfpc_model,'kon_akt_pp2a', 'Value', 1.875, 'ValueUnits', '1/(micromolarity*second)');
koff_akt_pp2a = addparameter(pfpc_model,'koff_akt_pp2a', 'Value', 7.2, 'ValueUnits', '1/second');
kcat_pp2a_akt = addparameter(pfpc_model,'kcat_pp2a_akt', 'Value', 1.8, 'ValueUnits', '1/second');

%Km = 0.8
kon_akt_pdk2 = addparameter(pfpc_model,'kon_akt_pdk2', 'Value', 125, 'ValueUnits', '1/(micromolarity*second)');
koff_akt_pdk2 = addparameter(pfpc_model,'koff_akt_pdk2', 'Value', 80, 'ValueUnits', '1/second');
kcat_pdk2 = addparameter(pfpc_model,'kcat_pdk2', 'Value', 20, 'ValueUnits', '1/second');



pdk1_reac_1 = addreaction(pfpc_model, 'PDK1 + pip3 -> PDK1_PIP3')
pdk1_reac_1_k = addkineticlaw(pdk1_reac_1,'MassAction');
set(pdk1_reac_1_k, 'ParameterVariableNames', {'kon_pdk_pip3'});

pdk1_reac_2 = addreaction(pfpc_model, 'PDK1_PIP3 -> PDK1 + pip3')
pdk1_reac_2_k = addkineticlaw(pdk1_reac_2,'MassAction');
set(pdk1_reac_2_k, 'ParameterVariableNames', {'koff_pdk_pip3'});


pdk1_reac_5 = addreaction(pfpc_model, 'AKT + PDK1_PIP3 -> AKT_PDK1')
pdk1_reac_5_k = addkineticlaw(pdk1_reac_5,'MassAction');
set(pdk1_reac_5_k, 'ParameterVariableNames', {'kon_akt_pdk1'});

pdk1_reac_6 = addreaction(pfpc_model, 'AKT_PDK1 -> AKT + PDK1_PIP3')
pdk1_reac_6_k = addkineticlaw(pdk1_reac_6,'MassAction');
set(pdk1_reac_6_k, 'ParameterVariableNames', {'koff_akt_pdk1'});

pdk1_reac_7 = addreaction(pfpc_model, 'AKT_PDK1 -> AKTp + PDK1_PIP3')
pdk1_reac_7_k = addkineticlaw(pdk1_reac_7,'MassAction');
set(pdk1_reac_7_k, 'ParameterVariableNames', {'kcat_pdk1'});


pdk1_reac_8 = addreaction(pfpc_model, 'AKTp + PP2A_p -> AKTp_PP2A')
pdk1_reac_8_k = addkineticlaw(pdk1_reac_8,'MassAction');
set(pdk1_reac_8_k, 'ParameterVariableNames', {'kon_akt_pp2a'});

pdk1_reac_9 = addreaction(pfpc_model, 'AKTp_PP2A -> AKTp + PP2A_p')
pdk1_reac_9_k = addkineticlaw(pdk1_reac_9,'MassAction');
set(pdk1_reac_9_k, 'ParameterVariableNames', {'koff_akt_pp2a'});

pdk1_reac_10 = addreaction(pfpc_model, 'AKTp_PP2A -> AKT + PP2A_p')
pdk1_reac_10_k = addkineticlaw(pdk1_reac_10,'MassAction');
set(pdk1_reac_10_k, 'ParameterVariableNames', {'kcat_pp2a_akt'});


pdk1_reac_11 = addreaction(pfpc_model, 'AKTp + PDK1_PIP3 -> AKTp_PDK1')
pdk1_reac_11_k = addkineticlaw(pdk1_reac_11,'MassAction');
set(pdk1_reac_11_k, 'ParameterVariableNames', {'kon_akt_pdk2'});

pdk1_reac_12 = addreaction(pfpc_model, 'AKTp_PDK1 -> AKTp + PDK1_PIP3')
pdk1_reac_12_k = addkineticlaw(pdk1_reac_12,'MassAction');
set(pdk1_reac_12_k, 'ParameterVariableNames', {'koff_akt_pdk2'});

pdk1_reac_13 = addreaction(pfpc_model, 'AKTp_PDK1 -> AKTp2 + PDK1_PIP3')
pdk1_reac_13_k = addkineticlaw(pdk1_reac_13,'MassAction');
set(pdk1_reac_13_k, 'ParameterVariableNames', {'kcat_pdk2'});


pdk1_reac_14 = addreaction(pfpc_model, 'AKTp2 + PP2A_p -> AKTp2_PP2A')
pdk1_reac_14_k = addkineticlaw(pdk1_reac_14,'MassAction');
set(pdk1_reac_14_k, 'ParameterVariableNames', {'kon_akt_pp2a'});

pdk1_reac_15 = addreaction(pfpc_model, 'AKTp2_PP2A -> AKTp2 + PP2A_p')
pdk1_reac_15_k = addkineticlaw(pdk1_reac_15,'MassAction');
set(pdk1_reac_15_k, 'ParameterVariableNames', {'koff_akt_pp2a'});

pdk1_reac_16 = addreaction(pfpc_model, 'AKTp2_PP2A -> AKTp + PP2A_p')
pdk1_reac_16_k = addkineticlaw(pdk1_reac_16,'MassAction');
set(pdk1_reac_16_k, 'ParameterVariableNames', {'kcat_pp2a_akt'});


%%%AKT_NOS Model
NOS = addspecies(cytosol, 'NOS', 'InitialAmount', 1, 'InitialAmountUnits', 'micromolarity');
NOSp = addspecies(cytosol, 'NOSp', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');
AKTp2_NOS = addspecies(cytosol, 'AKTp2_NOS', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');

NOSP = addspecies(cytosol, 'NOSP', 'InitialAmount', 0.4, 'InitialAmountUnits', 'micromolarity');
NOSp_NOSP = addspecies(cytosol, 'NOSp_NSOP', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');

NOp = addspecies(cytosol, 'NOp', 'InitialAmount', 0, 'InitialAmountUnits', 'micromolarity');

kon_nos_akt = addparameter(pfpc_model,'kon_nos_akt', 'Value', 250, 'ValueUnits', '1/(micromolarity*second)');
koff_nos_akt = addparameter(pfpc_model,'koff_nos_akt', 'Value', 4, 'ValueUnits', '1/second');
kcat_akt = addparameter(pfpc_model,'kcat_akt', 'Value', 1, 'ValueUnits', '1/second');

kon_nos_nosp = addparameter(pfpc_model,'kon_nos_nosp', 'Value', 250, 'ValueUnits', '1/(micromolarity*second)');
koff_nos_nosp = addparameter(pfpc_model,'koff_nos_nosp', 'Value', 4, 'ValueUnits', '1/second');
kcat_nosp = addparameter(pfpc_model,'kcat_nosp', 'Value', 1, 'ValueUnits', '1/second');

kcat_nos = addparameter(pfpc_model,'kcat_nos', 'Value', 0.2, 'ValueUnits', '1/second', 'ConstantValue', false);
kdiff_no = addparameter(pfpc_model,'kdiff_no', 'Value', 0.1, 'ValueUnits', '1/second');

%Km = 0.02
pdk1_reac_17 = addreaction(pfpc_model, 'NOS + AKTp2 -> AKTp2_NOS')
pdk1_reac_17_k = addkineticlaw(pdk1_reac_17,'MassAction');
set(pdk1_reac_17_k, 'ParameterVariableNames', {'kon_nos_akt'});

pdk1_reac_18 = addreaction(pfpc_model, 'AKTp2_NOS -> NOS + AKTp2')
pdk1_reac_18_k = addkineticlaw(pdk1_reac_18,'MassAction');
set(pdk1_reac_18_k, 'ParameterVariableNames', {'koff_nos_akt'});

pdk1_reac_19 = addreaction(pfpc_model, 'AKTp2_NOS -> NOSp + AKTp2')
pdk1_reac_19_k = addkineticlaw(pdk1_reac_19,'MassAction');
set(pdk1_reac_19_k, 'ParameterVariableNames', {'kcat_akt'});

%Km = 0.02
pdk1_reac_20 = addreaction(pfpc_model, 'NOSp + NOSP -> NOSp_NOSP')
pdk1_reac_20_k = addkineticlaw(pdk1_reac_20,'MassAction');
set(pdk1_reac_20_k, 'ParameterVariableNames', {'kon_nos_nosp'});

pdk1_reac_21 = addreaction(pfpc_model, 'NOSp_NOSP -> NOSp + NOSP')
pdk1_reac_21_k = addkineticlaw(pdk1_reac_21,'MassAction');
set(pdk1_reac_21_k, 'ParameterVariableNames', {'koff_nos_nosp'});

pdk1_reac_22 = addreaction(pfpc_model, 'NOSp_NOSP -> NOS + NOSP')
pdk1_reac_22_k = addkineticlaw(pdk1_reac_22,'MassAction');
set(pdk1_reac_22_k, 'ParameterVariableNames', {'kcat_nosp'});

%%NO synthesis by active NOS
pdk1_reac_23 = addreaction(pfpc_model, 'NOSp -> NOSp + NOp')
pdk1_reac_23_k = addkineticlaw(pdk1_reac_23,'MassAction');
set(pdk1_reac_23_k, 'ParameterVariableNames', {'kcat_nos'});

pdk1_reac_24 = addreaction(pfpc_model, 'NOp -> NO')
pdk1_reac_24_k = addkineticlaw(pdk1_reac_24,'MassAction');
set(pdk1_reac_24_k, 'ParameterVariableNames', {'kdiff_no'});



%%%%%%PULSES FOR CALCIUM AND GLUTAMATE FROM PF ACTIVITY%%%%%%%%%%%
glu_pulse = addparameter(pfpc_model,'glu_pulse', 'Value', 5300, 'ValueUnits', 'micromole/second', 'ConstantValue', false);
ca_pulse = addparameter(pfpc_model,'ca_pulse', 'Value', 2000, 'ValueUnits', 'micromole/second', 'ConstantValue', false);

%%%%Events for specific experiments. Do not uncomment without reason.
%%%%%%PF PULSES%%%%%%
%num_pulse = 200;
%pulse_times = [1:1:num_pulse]+10000;
%for i = 1:num_pulse
%	pulse_time = pulse_times(i);
%	pf_pulse_ltp(pulse_time)
%end

%removal of CaMKII kinase activity after 10 minutes
%eventObj7 = addevent(pfpc_model, 'time >= 10600', 'kon_camkii_pde = 0')
%removal of CaMKII kinase activity at beginning of induction
%eventObj7 = addevent(pfpc_model, 'time >= 10000', 'kon_camkii_pde = 0')
%eventObj7 = addevent(pfpc_model, 'time >= 10000', 'kon_ampar_pick = 56')
%eventObj8 = addevent(pfpc_model, 'time >= 10030', 'kon_ampar_pick = 7')

%%%%Inhibit PP2A at induction time
%eventObj7 = addevent(pfpc_model, 'time >= 10000', 'kon_mek_pp2a = 0')
%eventObj8 = addevent(pfpc_model, 'time >= 10000', 'kon_pla2_pp2a = 0')
%eventObj9 = addevent(pfpc_model, 'time >= 10000', 'kon_camkii_pp2a = 0')
%eventObj10 = addevent(pfpc_model, 'time >= 10000', 'kon_camkii_pp2a_2 = 0')
%eventObj11 = addevent(pfpc_model, 'time >= 10000', 'kon_ampar_pp2a = 0')

%%%%Inhibit PP1 at induction time
%eventObj8 = addevent(pfpc_model, 'time >= 10000', 'kon_pla2_pp1 = 0')
%eventObj9 = addevent(pfpc_model, 'time >= 10000', 'kon_camkii_pp1 = 0')
%eventObj10 = addevent(pfpc_model, 'time >= 10000', 'kon_camkii_pp1_2 = 0')
%eventObj11 = addevent(pfpc_model, 'time >= 10000', 'kon_ampar_pp1 = 0')

%%%%Inhibit calcineurin at induction time
%eventObj11 = addevent(pfpc_model, 'time >= 10000', 'kon_can_cam = 0')
%eventObj12 = addevent(pfpc_model, 'time >= 10000', 'kon_can_pgsub = 0')
%eventObj13 = addevent(pfpc_model, 'time >= 10000', 'kon_can_d32 = 0')


%configure simulation
config = addconfigset(pfpc_model, 'config');
config.StopTime = 20000;
config.SolverOptions.AbsoluteToleranceScaling = false;
%config.SolverType = 'ode45';
%set(config, 'MaximumNumberOfLogs', 100000)

%Species to log
ca_log = pfpc_model.Species(1);
pde_act_log = pfpc_model.Species(64);
no_log = pfpc_model.Species(94);
camkii_log = pfpc_model.Species(115);
pp2a_act_log = pfpc_model.Species(116);
erk_act_log = pfpc_model.Species(168);
mek_p_log = pfpc_model.Species(169);
aa_log = pfpc_model.Species(171);
pla2_act_log = pfpc_model.Species(198);
pkc_act_log = pfpc_model.Species(236);
ampar_psd_log = pfpc_model.Species(361);
ampar_cyt_log = pfpc_model.Species(362);
nsf_no_log = pfpc_model.Species(365);
nsfp_log = pfpc_model.Species(375);
cam_act_log = pfpc_model.Species(27);
camp_log = pfpc_model.Species(269);
pp1_log = pfpc_model.Species(23);
pka_log = pfpc_model.Species(396);
dag_log = pfpc_model.Species(203);
ip3_log = pfpc_model.Species(264);

%config.RuntimeOptions.StatesToLog = [ca_log pde_act_log no_log camkii_log pp2a_act_log erk_act_log mek_p_log aa_log pla2_act_log pkc_act_log ampar_psd_log ampar_cyt_log nsf_no_log nsfp_log cam_act_log camp_log pp1_log pka_log dag_log ip3_log];

setactiveconfigset(pfpc_model,config);

%config.RuntimeOptions.StatesToLog = {'Ca','ampar_psd', 'PDE_active', 'NO', 'CaMKII_Auton', 'PP2A_active', 'ERK_act', 'MEK_phospho', 'AA', 'PLA2_active', 'PKC_active','ampar_cytosol', 'nsf_no', 'nsfp', 'DAG', 'IP3', 'IP3R_open', 'CaMKII_305p'};

%%%%SIMULATION

%simDataObj = sbiosimulate(pfpc_model);
%sbioplot(simDataObj)
