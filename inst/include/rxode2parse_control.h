#pragma once
#ifndef __rxode2parse_control_H__
#define __rxode2parse_control_H__

#define RxMv_params 0
#define RxMv_lhs 1
#define RxMv_state 2
#define RxMv_trans 3
#define RxMv_model 4
#define RxMv_ini 5
#define RxMv_dfdy 6
#define RxMv_sens 7
#define RxMv_state_ignore 8
#define RxMv_version 9
#define RxMv_normal_state 10
#define RxMv_needSort 11
#define RxMv_nMtime 12
#define RxMv_extraCmt 13
#define RxMv_stateExtra 14
#define RxMv_dvid 15
#define RxMv_indLin 16
#define RxMv_flags 17
#define RxMv_slhs 18
#define RxMv_alag 19
#define RxMv_timeId 20
#define RxMv_md5 21
#define RxMvFlag_ncmt 0
#define RxMvFlag_ka 1
#define RxMvFlag_linB 2
#define RxMvFlag_maxeta 3
#define RxMvFlag_maxtheta 4
#define RxMvFlag_hasCmt 5
#define RxMvFlag_linCmt 6
#define RxMvFlag_linCmtFlg 7
#define RxMvFlag_nIndSim 8
#define RxMvFlag_simflg 9
#define RxMvFlag_thread 10
#define RxMvFlag_nLlik 11

#define RxMvTrans_lib_name 0
#define RxMvTrans_jac 1
#define RxMvTrans_prefix 2
#define RxMvTrans_dydt 3
#define RxMvTrans_calc_jac 4
#define RxMvTrans_calc_lhs 5
#define RxMvTrans_model_vars 6
#define RxMvTrans_theta 7
#define RxMvTrans_inis 8
#define RxMvTrans_dydt_lsoda 9
#define RxMvTrans_calc_jac_lsoda 10
#define RxMvTrans_ode_solver_solvedata 11
#define RxMvTrans_ode_solver_get_solvedata 12
#define RxMvTrans_dydt_liblsoda 13
#define RxMvTrans_F 14
#define RxMvTrans_Lag 15
#define RxMvTrans_Rate 16
#define RxMvTrans_Dur 17
#define RxMvTrans_mtime 18
#define RxMvTrans_assignFuns 19
#define RxMvTrans_ME 20
#define RxMvTrans_IndF 21
#define RxTrans_ndose 0
#define RxTrans_nobs 1
#define RxTrans_nid 2
#define RxTrans_cov1 3
#define RxTrans_covParPos 4
#define RxTrans_covParPosTV 5
#define RxTrans_sub0 6
#define RxTrans_baseSize 7
#define RxTrans_nTv 8
#define RxTrans_lst 9
#define RxTrans_nme 10
#define RxTrans_covParPos0 11
#define RxTrans_covUnits 12
#define RxTrans_pars 13
#define RxTrans_allBolus 14
#define RxTrans_allInf 15
#define RxTrans_mxCmt 16
#define RxTrans_lib_name 17
#define RxTrans_addCmt 18
#define RxTrans_cmtInfo 19
#define RxTrans_idLvl 20
#define RxTrans_allTimeVar 21
#define RxTrans_keepDosingOnly 22
#define RxTrans_censAdd 23
#define RxTrans_limitAdd 24
#define RxTrans_levelInfo 25
#define RxTrans_idInfo 26
#define RxTrans_maxShift 27
#define RxTrans_keepL 28
#define RxTransNames CharacterVector _en(29);_en[0]="ndose";_en[1]="nobs";_en[2]="nid";_en[3]="cov1";_en[4]="covParPos";_en[5]="covParPosTV";_en[6]="sub0";_en[7]="baseSize";_en[8]="nTv";_en[9]="lst";_en[10]="nme";_en[11]="covParPos0";_en[12]="covUnits";_en[13]="pars";_en[14]="allBolus";_en[15]="allInf";_en[16]="mxCmt";_en[17]="lib_name";_en[18]="addCmt";_en[19]="cmtInfo";_en[20]="idLvl";_en[21]="allTimeVar";_en[22]="keepDosingOnly";_en[23]="censAdd";_en[24]="limitAdd";_en[25]="levelInfo";_en[26]="idInfo";_en[27]="maxShift";_en[28]="keepL";e.names() = _en;

#endif // __rxode2parse_control_H__
