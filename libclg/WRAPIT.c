


#include <stdio.h>
/*
* The following are the required NCAR Graphics include files.
* They should be located in ${NCARG_ROOT}/include
*/
#include <ncarg/hlu/hlu.h>
#include <ncarg/hlu/NresDB.h>
#include <ncarg/ncl/defs.h>
#include <ncarg/ncl/NclDataDefs.h>
#include <ncarg/ncl/NclBuiltInSupport.h>
#include <ncarg/gks.h>
#include <ncarg/ncl/NclBuiltIns.h>



extern void NGCALLF(qsort,QSORT)();


NhlErrorTypes qsort_W( void ) {
	int i;
	float *arr;
	int *n_elem;
	long arr_dimsizes[NCL_MAX_DIMENSIONS];
	int arr_ndims;

	n_elem = (int*) NclGetArgValue(
		1,
		2,
		NULL,
		NULL,
		NULL,
		NULL,
		NULL,
		1);


	arr = (float*) NclGetArgValue(
		0,
		2,
		&arr_ndims,
		arr_dimsizes,
		NULL,
		NULL,
		NULL,
		1);

	if(*n_elem != (int)arr_dimsizes[0]) {
		NhlPError(NhlFATAL,NhlEUNKNOWN,"qsort: dimension size of dimension (0) of arr must be equal to the value of n_elem");
		return(NhlFATAL);
	}
	NGCALLF(qsort,QSORT)(arr,n_elem);

	return(NhlNOERROR);
}
extern void NGCALLF(search_index_binary,SEARCH_INDEX_BINARY)();


NhlErrorTypes search_index_binary_W( void ) {
	int i;
	float *sorted_arr;
	int *n_elem;
	float *key;
	int *res_idx;
	long sorted_arr_dimsizes[NCL_MAX_DIMENSIONS];
	int sorted_arr_ndims;

	res_idx = (int*) NclGetArgValue(
		3,
		4,
		NULL,
		NULL,
		NULL,
		NULL,
		NULL,
		1);


	key = (float*) NclGetArgValue(
		2,
		4,
		NULL,
		NULL,
		NULL,
		NULL,
		NULL,
		1);


	n_elem = (int*) NclGetArgValue(
		1,
		4,
		NULL,
		NULL,
		NULL,
		NULL,
		NULL,
		1);


	sorted_arr = (float*) NclGetArgValue(
		0,
		4,
		&sorted_arr_ndims,
		sorted_arr_dimsizes,
		NULL,
		NULL,
		NULL,
		1);

	if(*n_elem != (int)sorted_arr_dimsizes[0]) {
		NhlPError(NhlFATAL,NhlEUNKNOWN,"search_index_binary: dimension size of dimension (0) of sorted_arr must be equal to the value of n_elem");
		return(NhlFATAL);
	}
	NGCALLF(search_index_binary,SEARCH_INDEX_BINARY)(sorted_arr,n_elem,key,res_idx);

	return(NhlNOERROR);
}
extern void NGCALLF(compute_pr_ets_us_cwrf,COMPUTE_PR_ETS_US_CWRF)();


NhlErrorTypes compute_pr_ets_us_cwrf_W( void ) {
	int i;
	float *pr_obs;
	float *pr_mod;
	int *N_ETS;
	float *lcc;
	float *us_landmask;
	void *buffer_zone;
	void *ocean;
	void *us_land_only;
	int *NX;
	int *NY;
	float *ets;
	long pr_obs_dimsizes[NCL_MAX_DIMENSIONS];
	int pr_obs_ndims;
	long pr_mod_dimsizes[NCL_MAX_DIMENSIONS];
	int pr_mod_ndims;
	long lcc_dimsizes[NCL_MAX_DIMENSIONS];
	int lcc_ndims;
	long us_landmask_dimsizes[NCL_MAX_DIMENSIONS];
	int us_landmask_ndims;
	long ets_dimsizes[NCL_MAX_DIMENSIONS];
	int ets_ndims;

	ets = (float*) NclGetArgValue(
		10,
		11,
		&ets_ndims,
		ets_dimsizes,
		NULL,
		NULL,
		NULL,
		1);


	NY = (int*) NclGetArgValue(
		9,
		11,
		NULL,
		NULL,
		NULL,
		NULL,
		NULL,
		1);


	NX = (int*) NclGetArgValue(
		8,
		11,
		NULL,
		NULL,
		NULL,
		NULL,
		NULL,
		1);


	us_land_only = (void*) NclGetArgValue(
		7,
		11,
		NULL,
		NULL,
		NULL,
		NULL,
		NULL,
		1);


	ocean = (void*) NclGetArgValue(
		6,
		11,
		NULL,
		NULL,
		NULL,
		NULL,
		NULL,
		1);


	buffer_zone = (void*) NclGetArgValue(
		5,
		11,
		NULL,
		NULL,
		NULL,
		NULL,
		NULL,
		1);


	us_landmask = (float*) NclGetArgValue(
		4,
		11,
		&us_landmask_ndims,
		us_landmask_dimsizes,
		NULL,
		NULL,
		NULL,
		1);


	lcc = (float*) NclGetArgValue(
		3,
		11,
		&lcc_ndims,
		lcc_dimsizes,
		NULL,
		NULL,
		NULL,
		1);


	N_ETS = (int*) NclGetArgValue(
		2,
		11,
		NULL,
		NULL,
		NULL,
		NULL,
		NULL,
		1);


	pr_mod = (float*) NclGetArgValue(
		1,
		11,
		&pr_mod_ndims,
		pr_mod_dimsizes,
		NULL,
		NULL,
		NULL,
		1);


	pr_obs = (float*) NclGetArgValue(
		0,
		11,
		&pr_obs_ndims,
		pr_obs_dimsizes,
		NULL,
		NULL,
		NULL,
		1);

	if(*NX != (int)pr_obs_dimsizes[1]) {
		NhlPError(NhlFATAL,NhlEUNKNOWN,"compute_pr_ets_us_cwrf: dimension size of dimension (1) of pr_obs must be equal to the value of NX");
		return(NhlFATAL);
	}
	if(*NY != (int)pr_obs_dimsizes[0]) {
		NhlPError(NhlFATAL,NhlEUNKNOWN,"compute_pr_ets_us_cwrf: dimension size of dimension (0) of pr_obs must be equal to the value of NY");
		return(NhlFATAL);
	}
	if(*NX != (int)pr_mod_dimsizes[1]) {
		NhlPError(NhlFATAL,NhlEUNKNOWN,"compute_pr_ets_us_cwrf: dimension size of dimension (1) of pr_mod must be equal to the value of NX");
		return(NhlFATAL);
	}
	if(*NY != (int)pr_mod_dimsizes[0]) {
		NhlPError(NhlFATAL,NhlEUNKNOWN,"compute_pr_ets_us_cwrf: dimension size of dimension (0) of pr_mod must be equal to the value of NY");
		return(NhlFATAL);
	}
	if(*NX != (int)lcc_dimsizes[1]) {
		NhlPError(NhlFATAL,NhlEUNKNOWN,"compute_pr_ets_us_cwrf: dimension size of dimension (1) of lcc must be equal to the value of NX");
		return(NhlFATAL);
	}
	if(*NY != (int)lcc_dimsizes[0]) {
		NhlPError(NhlFATAL,NhlEUNKNOWN,"compute_pr_ets_us_cwrf: dimension size of dimension (0) of lcc must be equal to the value of NY");
		return(NhlFATAL);
	}
	if(*NX != (int)us_landmask_dimsizes[1]) {
		NhlPError(NhlFATAL,NhlEUNKNOWN,"compute_pr_ets_us_cwrf: dimension size of dimension (1) of us_landmask must be equal to the value of NX");
		return(NhlFATAL);
	}
	if(*NY != (int)us_landmask_dimsizes[0]) {
		NhlPError(NhlFATAL,NhlEUNKNOWN,"compute_pr_ets_us_cwrf: dimension size of dimension (0) of us_landmask must be equal to the value of NY");
		return(NhlFATAL);
	}
	if(*N_ETS != (int)ets_dimsizes[0]) {
		NhlPError(NhlFATAL,NhlEUNKNOWN,"compute_pr_ets_us_cwrf: dimension size of dimension (0) of ets must be equal to the value of N_ETS");
		return(NhlFATAL);
	}
	NGCALLF(compute_pr_ets_us_cwrf,COMPUTE_PR_ETS_US_CWRF)(pr_obs,pr_mod,N_ETS,lcc,us_landmask,buffer_zone,ocean,us_land_only,NX,NY,ets);

	return(NhlNOERROR);
}
extern void NGCALLF(write_out_5dayavg2monthlyavg_soda_mom,WRITE_OUT_5DAYAVG2MONTHLYAVG_SODA_MOM)();


NhlErrorTypes write_out_5dayAvg2monthlyAvg_SODA_MOM_W( void ) {
	int i;
	int *curr_yr;
	int *curr_mon;
	float *avg_temp;
	float *avg_salt;
	float *avg_u;
	float *avg_v;
	float *avg_sea_level;
	float *avg_wt;
	float *avg_tau_x;
	float *avg_tau_y;
	int *XT_OCEAN;
	int *YT_OCEAN;
	int *ST_OCEAN;
	long avg_temp_dimsizes[NCL_MAX_DIMENSIONS];
	int avg_temp_ndims;
	long avg_salt_dimsizes[NCL_MAX_DIMENSIONS];
	int avg_salt_ndims;
	long avg_u_dimsizes[NCL_MAX_DIMENSIONS];
	int avg_u_ndims;
	long avg_v_dimsizes[NCL_MAX_DIMENSIONS];
	int avg_v_ndims;
	long avg_sea_level_dimsizes[NCL_MAX_DIMENSIONS];
	int avg_sea_level_ndims;
	long avg_wt_dimsizes[NCL_MAX_DIMENSIONS];
	int avg_wt_ndims;
	long avg_tau_x_dimsizes[NCL_MAX_DIMENSIONS];
	int avg_tau_x_ndims;
	long avg_tau_y_dimsizes[NCL_MAX_DIMENSIONS];
	int avg_tau_y_ndims;

	ST_OCEAN = (int*) NclGetArgValue(
		12,
		13,
		NULL,
		NULL,
		NULL,
		NULL,
		NULL,
		1);


	YT_OCEAN = (int*) NclGetArgValue(
		11,
		13,
		NULL,
		NULL,
		NULL,
		NULL,
		NULL,
		1);


	XT_OCEAN = (int*) NclGetArgValue(
		10,
		13,
		NULL,
		NULL,
		NULL,
		NULL,
		NULL,
		1);


	avg_tau_y = (float*) NclGetArgValue(
		9,
		13,
		&avg_tau_y_ndims,
		avg_tau_y_dimsizes,
		NULL,
		NULL,
		NULL,
		1);


	avg_tau_x = (float*) NclGetArgValue(
		8,
		13,
		&avg_tau_x_ndims,
		avg_tau_x_dimsizes,
		NULL,
		NULL,
		NULL,
		1);


	avg_wt = (float*) NclGetArgValue(
		7,
		13,
		&avg_wt_ndims,
		avg_wt_dimsizes,
		NULL,
		NULL,
		NULL,
		1);


	avg_sea_level = (float*) NclGetArgValue(
		6,
		13,
		&avg_sea_level_ndims,
		avg_sea_level_dimsizes,
		NULL,
		NULL,
		NULL,
		1);


	avg_v = (float*) NclGetArgValue(
		5,
		13,
		&avg_v_ndims,
		avg_v_dimsizes,
		NULL,
		NULL,
		NULL,
		1);


	avg_u = (float*) NclGetArgValue(
		4,
		13,
		&avg_u_ndims,
		avg_u_dimsizes,
		NULL,
		NULL,
		NULL,
		1);


	avg_salt = (float*) NclGetArgValue(
		3,
		13,
		&avg_salt_ndims,
		avg_salt_dimsizes,
		NULL,
		NULL,
		NULL,
		1);


	avg_temp = (float*) NclGetArgValue(
		2,
		13,
		&avg_temp_ndims,
		avg_temp_dimsizes,
		NULL,
		NULL,
		NULL,
		1);


	curr_mon = (int*) NclGetArgValue(
		1,
		13,
		NULL,
		NULL,
		NULL,
		NULL,
		NULL,
		1);


	curr_yr = (int*) NclGetArgValue(
		0,
		13,
		NULL,
		NULL,
		NULL,
		NULL,
		NULL,
		1);

	if(*XT_OCEAN != (int)avg_temp_dimsizes[3]) {
		NhlPError(NhlFATAL,NhlEUNKNOWN,"write_out_5dayAvg2monthlyAvg_SODA_MOM: dimension size of dimension (3) of avg_temp must be equal to the value of XT_OCEAN");
		return(NhlFATAL);
	}
	if(*YT_OCEAN != (int)avg_temp_dimsizes[2]) {
		NhlPError(NhlFATAL,NhlEUNKNOWN,"write_out_5dayAvg2monthlyAvg_SODA_MOM: dimension size of dimension (2) of avg_temp must be equal to the value of YT_OCEAN");
		return(NhlFATAL);
	}
	if(*ST_OCEAN != (int)avg_temp_dimsizes[1]) {
		NhlPError(NhlFATAL,NhlEUNKNOWN,"write_out_5dayAvg2monthlyAvg_SODA_MOM: dimension size of dimension (1) of avg_temp must be equal to the value of ST_OCEAN");
		return(NhlFATAL);
	}
	if(*XT_OCEAN != (int)avg_salt_dimsizes[3]) {
		NhlPError(NhlFATAL,NhlEUNKNOWN,"write_out_5dayAvg2monthlyAvg_SODA_MOM: dimension size of dimension (3) of avg_salt must be equal to the value of XT_OCEAN");
		return(NhlFATAL);
	}
	if(*YT_OCEAN != (int)avg_salt_dimsizes[2]) {
		NhlPError(NhlFATAL,NhlEUNKNOWN,"write_out_5dayAvg2monthlyAvg_SODA_MOM: dimension size of dimension (2) of avg_salt must be equal to the value of YT_OCEAN");
		return(NhlFATAL);
	}
	if(*ST_OCEAN != (int)avg_salt_dimsizes[1]) {
		NhlPError(NhlFATAL,NhlEUNKNOWN,"write_out_5dayAvg2monthlyAvg_SODA_MOM: dimension size of dimension (1) of avg_salt must be equal to the value of ST_OCEAN");
		return(NhlFATAL);
	}
	if(*XT_OCEAN != (int)avg_u_dimsizes[3]) {
		NhlPError(NhlFATAL,NhlEUNKNOWN,"write_out_5dayAvg2monthlyAvg_SODA_MOM: dimension size of dimension (3) of avg_u must be equal to the value of XT_OCEAN");
		return(NhlFATAL);
	}
	if(*YT_OCEAN != (int)avg_u_dimsizes[2]) {
		NhlPError(NhlFATAL,NhlEUNKNOWN,"write_out_5dayAvg2monthlyAvg_SODA_MOM: dimension size of dimension (2) of avg_u must be equal to the value of YT_OCEAN");
		return(NhlFATAL);
	}
	if(*ST_OCEAN != (int)avg_u_dimsizes[1]) {
		NhlPError(NhlFATAL,NhlEUNKNOWN,"write_out_5dayAvg2monthlyAvg_SODA_MOM: dimension size of dimension (1) of avg_u must be equal to the value of ST_OCEAN");
		return(NhlFATAL);
	}
	if(*XT_OCEAN != (int)avg_v_dimsizes[3]) {
		NhlPError(NhlFATAL,NhlEUNKNOWN,"write_out_5dayAvg2monthlyAvg_SODA_MOM: dimension size of dimension (3) of avg_v must be equal to the value of XT_OCEAN");
		return(NhlFATAL);
	}
	if(*YT_OCEAN != (int)avg_v_dimsizes[2]) {
		NhlPError(NhlFATAL,NhlEUNKNOWN,"write_out_5dayAvg2monthlyAvg_SODA_MOM: dimension size of dimension (2) of avg_v must be equal to the value of YT_OCEAN");
		return(NhlFATAL);
	}
	if(*ST_OCEAN != (int)avg_v_dimsizes[1]) {
		NhlPError(NhlFATAL,NhlEUNKNOWN,"write_out_5dayAvg2monthlyAvg_SODA_MOM: dimension size of dimension (1) of avg_v must be equal to the value of ST_OCEAN");
		return(NhlFATAL);
	}
	if(*XT_OCEAN != (int)avg_sea_level_dimsizes[2]) {
		NhlPError(NhlFATAL,NhlEUNKNOWN,"write_out_5dayAvg2monthlyAvg_SODA_MOM: dimension size of dimension (2) of avg_sea_level must be equal to the value of XT_OCEAN");
		return(NhlFATAL);
	}
	if(*YT_OCEAN != (int)avg_sea_level_dimsizes[1]) {
		NhlPError(NhlFATAL,NhlEUNKNOWN,"write_out_5dayAvg2monthlyAvg_SODA_MOM: dimension size of dimension (1) of avg_sea_level must be equal to the value of YT_OCEAN");
		return(NhlFATAL);
	}
	if(*XT_OCEAN != (int)avg_wt_dimsizes[3]) {
		NhlPError(NhlFATAL,NhlEUNKNOWN,"write_out_5dayAvg2monthlyAvg_SODA_MOM: dimension size of dimension (3) of avg_wt must be equal to the value of XT_OCEAN");
		return(NhlFATAL);
	}
	if(*YT_OCEAN != (int)avg_wt_dimsizes[2]) {
		NhlPError(NhlFATAL,NhlEUNKNOWN,"write_out_5dayAvg2monthlyAvg_SODA_MOM: dimension size of dimension (2) of avg_wt must be equal to the value of YT_OCEAN");
		return(NhlFATAL);
	}
	if(*ST_OCEAN != (int)avg_wt_dimsizes[1]) {
		NhlPError(NhlFATAL,NhlEUNKNOWN,"write_out_5dayAvg2monthlyAvg_SODA_MOM: dimension size of dimension (1) of avg_wt must be equal to the value of ST_OCEAN");
		return(NhlFATAL);
	}
	if(*XT_OCEAN != (int)avg_tau_x_dimsizes[2]) {
		NhlPError(NhlFATAL,NhlEUNKNOWN,"write_out_5dayAvg2monthlyAvg_SODA_MOM: dimension size of dimension (2) of avg_tau_x must be equal to the value of XT_OCEAN");
		return(NhlFATAL);
	}
	if(*YT_OCEAN != (int)avg_tau_x_dimsizes[1]) {
		NhlPError(NhlFATAL,NhlEUNKNOWN,"write_out_5dayAvg2monthlyAvg_SODA_MOM: dimension size of dimension (1) of avg_tau_x must be equal to the value of YT_OCEAN");
		return(NhlFATAL);
	}
	if(*XT_OCEAN != (int)avg_tau_y_dimsizes[2]) {
		NhlPError(NhlFATAL,NhlEUNKNOWN,"write_out_5dayAvg2monthlyAvg_SODA_MOM: dimension size of dimension (2) of avg_tau_y must be equal to the value of XT_OCEAN");
		return(NhlFATAL);
	}
	if(*YT_OCEAN != (int)avg_tau_y_dimsizes[1]) {
		NhlPError(NhlFATAL,NhlEUNKNOWN,"write_out_5dayAvg2monthlyAvg_SODA_MOM: dimension size of dimension (1) of avg_tau_y must be equal to the value of YT_OCEAN");
		return(NhlFATAL);
	}
	NGCALLF(write_out_5dayavg2monthlyavg_soda_mom,WRITE_OUT_5DAYAVG2MONTHLYAVG_SODA_MOM)(curr_yr,curr_mon,avg_temp,avg_salt,avg_u,avg_v,avg_sea_level,avg_wt,avg_tau_x,avg_tau_y,XT_OCEAN,YT_OCEAN,ST_OCEAN);

	return(NhlNOERROR);
}
extern void NGCALLF(write_out_5dayavg2monthlyavg_soda_mom_ice,WRITE_OUT_5DAYAVG2MONTHLYAVG_SODA_MOM_ICE)();


NhlErrorTypes write_out_5dayAvg2monthlyAvg_SODA_MOM_ice_W( void ) {
	int i;
	int *curr_yr;
	int *curr_mon;
	float *avg_CN;
	float *avg_HI;
	int *xt;
	int *yt;
	int *ct;
	long avg_CN_dimsizes[NCL_MAX_DIMENSIONS];
	int avg_CN_ndims;
	long avg_HI_dimsizes[NCL_MAX_DIMENSIONS];
	int avg_HI_ndims;

	ct = (int*) NclGetArgValue(
		6,
		7,
		NULL,
		NULL,
		NULL,
		NULL,
		NULL,
		1);


	yt = (int*) NclGetArgValue(
		5,
		7,
		NULL,
		NULL,
		NULL,
		NULL,
		NULL,
		1);


	xt = (int*) NclGetArgValue(
		4,
		7,
		NULL,
		NULL,
		NULL,
		NULL,
		NULL,
		1);


	avg_HI = (float*) NclGetArgValue(
		3,
		7,
		&avg_HI_ndims,
		avg_HI_dimsizes,
		NULL,
		NULL,
		NULL,
		1);


	avg_CN = (float*) NclGetArgValue(
		2,
		7,
		&avg_CN_ndims,
		avg_CN_dimsizes,
		NULL,
		NULL,
		NULL,
		1);


	curr_mon = (int*) NclGetArgValue(
		1,
		7,
		NULL,
		NULL,
		NULL,
		NULL,
		NULL,
		1);


	curr_yr = (int*) NclGetArgValue(
		0,
		7,
		NULL,
		NULL,
		NULL,
		NULL,
		NULL,
		1);

	if(*xt != (int)avg_CN_dimsizes[3]) {
		NhlPError(NhlFATAL,NhlEUNKNOWN,"write_out_5dayAvg2monthlyAvg_SODA_MOM_ice: dimension size of dimension (3) of avg_CN must be equal to the value of xt");
		return(NhlFATAL);
	}
	if(*yt != (int)avg_CN_dimsizes[2]) {
		NhlPError(NhlFATAL,NhlEUNKNOWN,"write_out_5dayAvg2monthlyAvg_SODA_MOM_ice: dimension size of dimension (2) of avg_CN must be equal to the value of yt");
		return(NhlFATAL);
	}
	if(*ct != (int)avg_CN_dimsizes[1]) {
		NhlPError(NhlFATAL,NhlEUNKNOWN,"write_out_5dayAvg2monthlyAvg_SODA_MOM_ice: dimension size of dimension (1) of avg_CN must be equal to the value of ct");
		return(NhlFATAL);
	}
	if(*xt != (int)avg_HI_dimsizes[2]) {
		NhlPError(NhlFATAL,NhlEUNKNOWN,"write_out_5dayAvg2monthlyAvg_SODA_MOM_ice: dimension size of dimension (2) of avg_HI must be equal to the value of xt");
		return(NhlFATAL);
	}
	if(*yt != (int)avg_HI_dimsizes[1]) {
		NhlPError(NhlFATAL,NhlEUNKNOWN,"write_out_5dayAvg2monthlyAvg_SODA_MOM_ice: dimension size of dimension (1) of avg_HI must be equal to the value of yt");
		return(NhlFATAL);
	}
	NGCALLF(write_out_5dayavg2monthlyavg_soda_mom_ice,WRITE_OUT_5DAYAVG2MONTHLYAVG_SODA_MOM_ICE)(curr_yr,curr_mon,avg_CN,avg_HI,xt,yt,ct);

	return(NhlNOERROR);
}
extern void NGCALLF(write_out_5dayavg2monthlyavg_soda_mom_1979_1981,WRITE_OUT_5DAYAVG2MONTHLYAVG_SODA_MOM_1979_1981)();


NhlErrorTypes write_out_5dayAvg2monthlyAvg_SODA_MOM_1979_1981_W( void ) {
	int i;
	int *curr_yr;
	int *curr_mon;
	float *avg_temp;
	float *avg_salt;
	float *avg_u;
	float *avg_v;
	float *avg_ssh;
	float *avg_mlt;
	float *avg_mlp;
	float *avg_mls;
	float *avg_pbot;
	float *avg_wt;
	float *avg_prho;
	float *avg_tau_x;
	float *avg_tau_y;
	int *XT_OCEAN;
	int *YT_OCEAN;
	int *ST_OCEAN;
	long avg_temp_dimsizes[NCL_MAX_DIMENSIONS];
	int avg_temp_ndims;
	long avg_salt_dimsizes[NCL_MAX_DIMENSIONS];
	int avg_salt_ndims;
	long avg_u_dimsizes[NCL_MAX_DIMENSIONS];
	int avg_u_ndims;
	long avg_v_dimsizes[NCL_MAX_DIMENSIONS];
	int avg_v_ndims;
	long avg_ssh_dimsizes[NCL_MAX_DIMENSIONS];
	int avg_ssh_ndims;
	long avg_mlt_dimsizes[NCL_MAX_DIMENSIONS];
	int avg_mlt_ndims;
	long avg_mlp_dimsizes[NCL_MAX_DIMENSIONS];
	int avg_mlp_ndims;
	long avg_mls_dimsizes[NCL_MAX_DIMENSIONS];
	int avg_mls_ndims;
	long avg_pbot_dimsizes[NCL_MAX_DIMENSIONS];
	int avg_pbot_ndims;
	long avg_wt_dimsizes[NCL_MAX_DIMENSIONS];
	int avg_wt_ndims;
	long avg_prho_dimsizes[NCL_MAX_DIMENSIONS];
	int avg_prho_ndims;
	long avg_tau_x_dimsizes[NCL_MAX_DIMENSIONS];
	int avg_tau_x_ndims;
	long avg_tau_y_dimsizes[NCL_MAX_DIMENSIONS];
	int avg_tau_y_ndims;

	ST_OCEAN = (int*) NclGetArgValue(
		17,
		18,
		NULL,
		NULL,
		NULL,
		NULL,
		NULL,
		1);


	YT_OCEAN = (int*) NclGetArgValue(
		16,
		18,
		NULL,
		NULL,
		NULL,
		NULL,
		NULL,
		1);


	XT_OCEAN = (int*) NclGetArgValue(
		15,
		18,
		NULL,
		NULL,
		NULL,
		NULL,
		NULL,
		1);


	avg_tau_y = (float*) NclGetArgValue(
		14,
		18,
		&avg_tau_y_ndims,
		avg_tau_y_dimsizes,
		NULL,
		NULL,
		NULL,
		1);


	avg_tau_x = (float*) NclGetArgValue(
		13,
		18,
		&avg_tau_x_ndims,
		avg_tau_x_dimsizes,
		NULL,
		NULL,
		NULL,
		1);


	avg_prho = (float*) NclGetArgValue(
		12,
		18,
		&avg_prho_ndims,
		avg_prho_dimsizes,
		NULL,
		NULL,
		NULL,
		1);


	avg_wt = (float*) NclGetArgValue(
		11,
		18,
		&avg_wt_ndims,
		avg_wt_dimsizes,
		NULL,
		NULL,
		NULL,
		1);


	avg_pbot = (float*) NclGetArgValue(
		10,
		18,
		&avg_pbot_ndims,
		avg_pbot_dimsizes,
		NULL,
		NULL,
		NULL,
		1);


	avg_mls = (float*) NclGetArgValue(
		9,
		18,
		&avg_mls_ndims,
		avg_mls_dimsizes,
		NULL,
		NULL,
		NULL,
		1);


	avg_mlp = (float*) NclGetArgValue(
		8,
		18,
		&avg_mlp_ndims,
		avg_mlp_dimsizes,
		NULL,
		NULL,
		NULL,
		1);


	avg_mlt = (float*) NclGetArgValue(
		7,
		18,
		&avg_mlt_ndims,
		avg_mlt_dimsizes,
		NULL,
		NULL,
		NULL,
		1);


	avg_ssh = (float*) NclGetArgValue(
		6,
		18,
		&avg_ssh_ndims,
		avg_ssh_dimsizes,
		NULL,
		NULL,
		NULL,
		1);


	avg_v = (float*) NclGetArgValue(
		5,
		18,
		&avg_v_ndims,
		avg_v_dimsizes,
		NULL,
		NULL,
		NULL,
		1);


	avg_u = (float*) NclGetArgValue(
		4,
		18,
		&avg_u_ndims,
		avg_u_dimsizes,
		NULL,
		NULL,
		NULL,
		1);


	avg_salt = (float*) NclGetArgValue(
		3,
		18,
		&avg_salt_ndims,
		avg_salt_dimsizes,
		NULL,
		NULL,
		NULL,
		1);


	avg_temp = (float*) NclGetArgValue(
		2,
		18,
		&avg_temp_ndims,
		avg_temp_dimsizes,
		NULL,
		NULL,
		NULL,
		1);


	curr_mon = (int*) NclGetArgValue(
		1,
		18,
		NULL,
		NULL,
		NULL,
		NULL,
		NULL,
		1);


	curr_yr = (int*) NclGetArgValue(
		0,
		18,
		NULL,
		NULL,
		NULL,
		NULL,
		NULL,
		1);

	if(*XT_OCEAN != (int)avg_temp_dimsizes[3]) {
		NhlPError(NhlFATAL,NhlEUNKNOWN,"write_out_5dayAvg2monthlyAvg_SODA_MOM_1979_1981: dimension size of dimension (3) of avg_temp must be equal to the value of XT_OCEAN");
		return(NhlFATAL);
	}
	if(*YT_OCEAN != (int)avg_temp_dimsizes[2]) {
		NhlPError(NhlFATAL,NhlEUNKNOWN,"write_out_5dayAvg2monthlyAvg_SODA_MOM_1979_1981: dimension size of dimension (2) of avg_temp must be equal to the value of YT_OCEAN");
		return(NhlFATAL);
	}
	if(*ST_OCEAN != (int)avg_temp_dimsizes[1]) {
		NhlPError(NhlFATAL,NhlEUNKNOWN,"write_out_5dayAvg2monthlyAvg_SODA_MOM_1979_1981: dimension size of dimension (1) of avg_temp must be equal to the value of ST_OCEAN");
		return(NhlFATAL);
	}
	if(*XT_OCEAN != (int)avg_salt_dimsizes[3]) {
		NhlPError(NhlFATAL,NhlEUNKNOWN,"write_out_5dayAvg2monthlyAvg_SODA_MOM_1979_1981: dimension size of dimension (3) of avg_salt must be equal to the value of XT_OCEAN");
		return(NhlFATAL);
	}
	if(*YT_OCEAN != (int)avg_salt_dimsizes[2]) {
		NhlPError(NhlFATAL,NhlEUNKNOWN,"write_out_5dayAvg2monthlyAvg_SODA_MOM_1979_1981: dimension size of dimension (2) of avg_salt must be equal to the value of YT_OCEAN");
		return(NhlFATAL);
	}
	if(*ST_OCEAN != (int)avg_salt_dimsizes[1]) {
		NhlPError(NhlFATAL,NhlEUNKNOWN,"write_out_5dayAvg2monthlyAvg_SODA_MOM_1979_1981: dimension size of dimension (1) of avg_salt must be equal to the value of ST_OCEAN");
		return(NhlFATAL);
	}
	if(*XT_OCEAN != (int)avg_u_dimsizes[3]) {
		NhlPError(NhlFATAL,NhlEUNKNOWN,"write_out_5dayAvg2monthlyAvg_SODA_MOM_1979_1981: dimension size of dimension (3) of avg_u must be equal to the value of XT_OCEAN");
		return(NhlFATAL);
	}
	if(*YT_OCEAN != (int)avg_u_dimsizes[2]) {
		NhlPError(NhlFATAL,NhlEUNKNOWN,"write_out_5dayAvg2monthlyAvg_SODA_MOM_1979_1981: dimension size of dimension (2) of avg_u must be equal to the value of YT_OCEAN");
		return(NhlFATAL);
	}
	if(*ST_OCEAN != (int)avg_u_dimsizes[1]) {
		NhlPError(NhlFATAL,NhlEUNKNOWN,"write_out_5dayAvg2monthlyAvg_SODA_MOM_1979_1981: dimension size of dimension (1) of avg_u must be equal to the value of ST_OCEAN");
		return(NhlFATAL);
	}
	if(*XT_OCEAN != (int)avg_v_dimsizes[3]) {
		NhlPError(NhlFATAL,NhlEUNKNOWN,"write_out_5dayAvg2monthlyAvg_SODA_MOM_1979_1981: dimension size of dimension (3) of avg_v must be equal to the value of XT_OCEAN");
		return(NhlFATAL);
	}
	if(*YT_OCEAN != (int)avg_v_dimsizes[2]) {
		NhlPError(NhlFATAL,NhlEUNKNOWN,"write_out_5dayAvg2monthlyAvg_SODA_MOM_1979_1981: dimension size of dimension (2) of avg_v must be equal to the value of YT_OCEAN");
		return(NhlFATAL);
	}
	if(*ST_OCEAN != (int)avg_v_dimsizes[1]) {
		NhlPError(NhlFATAL,NhlEUNKNOWN,"write_out_5dayAvg2monthlyAvg_SODA_MOM_1979_1981: dimension size of dimension (1) of avg_v must be equal to the value of ST_OCEAN");
		return(NhlFATAL);
	}
	if(*XT_OCEAN != (int)avg_ssh_dimsizes[2]) {
		NhlPError(NhlFATAL,NhlEUNKNOWN,"write_out_5dayAvg2monthlyAvg_SODA_MOM_1979_1981: dimension size of dimension (2) of avg_ssh must be equal to the value of XT_OCEAN");
		return(NhlFATAL);
	}
	if(*YT_OCEAN != (int)avg_ssh_dimsizes[1]) {
		NhlPError(NhlFATAL,NhlEUNKNOWN,"write_out_5dayAvg2monthlyAvg_SODA_MOM_1979_1981: dimension size of dimension (1) of avg_ssh must be equal to the value of YT_OCEAN");
		return(NhlFATAL);
	}
	if(*XT_OCEAN != (int)avg_mlt_dimsizes[2]) {
		NhlPError(NhlFATAL,NhlEUNKNOWN,"write_out_5dayAvg2monthlyAvg_SODA_MOM_1979_1981: dimension size of dimension (2) of avg_mlt must be equal to the value of XT_OCEAN");
		return(NhlFATAL);
	}
	if(*YT_OCEAN != (int)avg_mlt_dimsizes[1]) {
		NhlPError(NhlFATAL,NhlEUNKNOWN,"write_out_5dayAvg2monthlyAvg_SODA_MOM_1979_1981: dimension size of dimension (1) of avg_mlt must be equal to the value of YT_OCEAN");
		return(NhlFATAL);
	}
	if(*XT_OCEAN != (int)avg_mlp_dimsizes[2]) {
		NhlPError(NhlFATAL,NhlEUNKNOWN,"write_out_5dayAvg2monthlyAvg_SODA_MOM_1979_1981: dimension size of dimension (2) of avg_mlp must be equal to the value of XT_OCEAN");
		return(NhlFATAL);
	}
	if(*YT_OCEAN != (int)avg_mlp_dimsizes[1]) {
		NhlPError(NhlFATAL,NhlEUNKNOWN,"write_out_5dayAvg2monthlyAvg_SODA_MOM_1979_1981: dimension size of dimension (1) of avg_mlp must be equal to the value of YT_OCEAN");
		return(NhlFATAL);
	}
	if(*XT_OCEAN != (int)avg_mls_dimsizes[2]) {
		NhlPError(NhlFATAL,NhlEUNKNOWN,"write_out_5dayAvg2monthlyAvg_SODA_MOM_1979_1981: dimension size of dimension (2) of avg_mls must be equal to the value of XT_OCEAN");
		return(NhlFATAL);
	}
	if(*YT_OCEAN != (int)avg_mls_dimsizes[1]) {
		NhlPError(NhlFATAL,NhlEUNKNOWN,"write_out_5dayAvg2monthlyAvg_SODA_MOM_1979_1981: dimension size of dimension (1) of avg_mls must be equal to the value of YT_OCEAN");
		return(NhlFATAL);
	}
	if(*XT_OCEAN != (int)avg_pbot_dimsizes[2]) {
		NhlPError(NhlFATAL,NhlEUNKNOWN,"write_out_5dayAvg2monthlyAvg_SODA_MOM_1979_1981: dimension size of dimension (2) of avg_pbot must be equal to the value of XT_OCEAN");
		return(NhlFATAL);
	}
	if(*YT_OCEAN != (int)avg_pbot_dimsizes[1]) {
		NhlPError(NhlFATAL,NhlEUNKNOWN,"write_out_5dayAvg2monthlyAvg_SODA_MOM_1979_1981: dimension size of dimension (1) of avg_pbot must be equal to the value of YT_OCEAN");
		return(NhlFATAL);
	}
	if(*XT_OCEAN != (int)avg_wt_dimsizes[3]) {
		NhlPError(NhlFATAL,NhlEUNKNOWN,"write_out_5dayAvg2monthlyAvg_SODA_MOM_1979_1981: dimension size of dimension (3) of avg_wt must be equal to the value of XT_OCEAN");
		return(NhlFATAL);
	}
	if(*YT_OCEAN != (int)avg_wt_dimsizes[2]) {
		NhlPError(NhlFATAL,NhlEUNKNOWN,"write_out_5dayAvg2monthlyAvg_SODA_MOM_1979_1981: dimension size of dimension (2) of avg_wt must be equal to the value of YT_OCEAN");
		return(NhlFATAL);
	}
	if(*ST_OCEAN != (int)avg_wt_dimsizes[1]) {
		NhlPError(NhlFATAL,NhlEUNKNOWN,"write_out_5dayAvg2monthlyAvg_SODA_MOM_1979_1981: dimension size of dimension (1) of avg_wt must be equal to the value of ST_OCEAN");
		return(NhlFATAL);
	}
	if(*XT_OCEAN != (int)avg_prho_dimsizes[3]) {
		NhlPError(NhlFATAL,NhlEUNKNOWN,"write_out_5dayAvg2monthlyAvg_SODA_MOM_1979_1981: dimension size of dimension (3) of avg_prho must be equal to the value of XT_OCEAN");
		return(NhlFATAL);
	}
	if(*YT_OCEAN != (int)avg_prho_dimsizes[2]) {
		NhlPError(NhlFATAL,NhlEUNKNOWN,"write_out_5dayAvg2monthlyAvg_SODA_MOM_1979_1981: dimension size of dimension (2) of avg_prho must be equal to the value of YT_OCEAN");
		return(NhlFATAL);
	}
	if(*ST_OCEAN != (int)avg_prho_dimsizes[1]) {
		NhlPError(NhlFATAL,NhlEUNKNOWN,"write_out_5dayAvg2monthlyAvg_SODA_MOM_1979_1981: dimension size of dimension (1) of avg_prho must be equal to the value of ST_OCEAN");
		return(NhlFATAL);
	}
	if(*XT_OCEAN != (int)avg_tau_x_dimsizes[2]) {
		NhlPError(NhlFATAL,NhlEUNKNOWN,"write_out_5dayAvg2monthlyAvg_SODA_MOM_1979_1981: dimension size of dimension (2) of avg_tau_x must be equal to the value of XT_OCEAN");
		return(NhlFATAL);
	}
	if(*YT_OCEAN != (int)avg_tau_x_dimsizes[1]) {
		NhlPError(NhlFATAL,NhlEUNKNOWN,"write_out_5dayAvg2monthlyAvg_SODA_MOM_1979_1981: dimension size of dimension (1) of avg_tau_x must be equal to the value of YT_OCEAN");
		return(NhlFATAL);
	}
	if(*XT_OCEAN != (int)avg_tau_y_dimsizes[2]) {
		NhlPError(NhlFATAL,NhlEUNKNOWN,"write_out_5dayAvg2monthlyAvg_SODA_MOM_1979_1981: dimension size of dimension (2) of avg_tau_y must be equal to the value of XT_OCEAN");
		return(NhlFATAL);
	}
	if(*YT_OCEAN != (int)avg_tau_y_dimsizes[1]) {
		NhlPError(NhlFATAL,NhlEUNKNOWN,"write_out_5dayAvg2monthlyAvg_SODA_MOM_1979_1981: dimension size of dimension (1) of avg_tau_y must be equal to the value of YT_OCEAN");
		return(NhlFATAL);
	}
	NGCALLF(write_out_5dayavg2monthlyavg_soda_mom_1979_1981,WRITE_OUT_5DAYAVG2MONTHLYAVG_SODA_MOM_1979_1981)(curr_yr,curr_mon,avg_temp,avg_salt,avg_u,avg_v,avg_ssh,avg_mlt,avg_mlp,avg_mls,avg_pbot,avg_wt,avg_prho,avg_tau_x,avg_tau_y,XT_OCEAN,YT_OCEAN,ST_OCEAN);

	return(NhlNOERROR);
}
extern void NGCALLF(cal_temp_fcs_obs_diff_dist,CAL_TEMP_FCS_OBS_DIFF_DIST)();


NhlErrorTypes cal_temp_fcs_obs_diff_dist_W( void ) {
	int i;
	int *yr;
	int *mon;
	int *day;
	int *jul_day;
	int *temp_tmp;

	temp_tmp = (int*) NclGetArgValue(
		4,
		5,
		NULL,
		NULL,
		NULL,
		NULL,
		NULL,
		1);


	jul_day = (int*) NclGetArgValue(
		3,
		5,
		NULL,
		NULL,
		NULL,
		NULL,
		NULL,
		1);


	day = (int*) NclGetArgValue(
		2,
		5,
		NULL,
		NULL,
		NULL,
		NULL,
		NULL,
		1);


	mon = (int*) NclGetArgValue(
		1,
		5,
		NULL,
		NULL,
		NULL,
		NULL,
		NULL,
		1);


	yr = (int*) NclGetArgValue(
		0,
		5,
		NULL,
		NULL,
		NULL,
		NULL,
		NULL,
		1);

	NGCALLF(cal_temp_fcs_obs_diff_dist,CAL_TEMP_FCS_OBS_DIFF_DIST)(yr,mon,day,jul_day,temp_tmp);

	return(NhlNOERROR);
}
extern void NGCALLF(cal_salt_fcs_obs_diff_dist,CAL_SALT_FCS_OBS_DIFF_DIST)();


NhlErrorTypes cal_salt_fcs_obs_diff_dist_W( void ) {
	int i;
	int *yr;
	int *mon;
	int *day;
	int *jul_day;
	int *salt_tmp;

	salt_tmp = (int*) NclGetArgValue(
		4,
		5,
		NULL,
		NULL,
		NULL,
		NULL,
		NULL,
		1);


	jul_day = (int*) NclGetArgValue(
		3,
		5,
		NULL,
		NULL,
		NULL,
		NULL,
		NULL,
		1);


	day = (int*) NclGetArgValue(
		2,
		5,
		NULL,
		NULL,
		NULL,
		NULL,
		NULL,
		1);


	mon = (int*) NclGetArgValue(
		1,
		5,
		NULL,
		NULL,
		NULL,
		NULL,
		NULL,
		1);


	yr = (int*) NclGetArgValue(
		0,
		5,
		NULL,
		NULL,
		NULL,
		NULL,
		NULL,
		1);

	NGCALLF(cal_salt_fcs_obs_diff_dist,CAL_SALT_FCS_OBS_DIFF_DIST)(yr,mon,day,jul_day,salt_tmp);

	return(NhlNOERROR);
}
extern void NGCALLF(sst_noaa2soda_1hr,SST_NOAA2SODA_1HR)();


NhlErrorTypes sst_noaa2soda_1hr_W( void ) {
	int i;
	int *ni;
	int *nj;
	int *time;
	float *lon_noaa;
	float *lat_noaa;
	float *sst_noaa;
	float *sses_bias;
	short *l2p_flags;
	byte *quality_level;
	float *sst_soda;
	int *n_sst_soda;
	long lon_noaa_dimsizes[NCL_MAX_DIMENSIONS];
	int lon_noaa_ndims;
	long lat_noaa_dimsizes[NCL_MAX_DIMENSIONS];
	int lat_noaa_ndims;
	long sst_noaa_dimsizes[NCL_MAX_DIMENSIONS];
	int sst_noaa_ndims;
	long sses_bias_dimsizes[NCL_MAX_DIMENSIONS];
	int sses_bias_ndims;
	long l2p_flags_dimsizes[NCL_MAX_DIMENSIONS];
	int l2p_flags_ndims;
	long quality_level_dimsizes[NCL_MAX_DIMENSIONS];
	int quality_level_ndims;

	n_sst_soda = (int*) NclGetArgValue(
		10,
		11,
		NULL,
		NULL,
		NULL,
		NULL,
		NULL,
		1);


	sst_soda = (float*) NclGetArgValue(
		9,
		11,
		NULL,
		NULL,
		NULL,
		NULL,
		NULL,
		1);


	quality_level = (byte*) NclGetArgValue(
		8,
		11,
		&quality_level_ndims,
		quality_level_dimsizes,
		NULL,
		NULL,
		NULL,
		1);


	l2p_flags = (short*) NclGetArgValue(
		7,
		11,
		&l2p_flags_ndims,
		l2p_flags_dimsizes,
		NULL,
		NULL,
		NULL,
		1);


	sses_bias = (float*) NclGetArgValue(
		6,
		11,
		&sses_bias_ndims,
		sses_bias_dimsizes,
		NULL,
		NULL,
		NULL,
		1);


	sst_noaa = (float*) NclGetArgValue(
		5,
		11,
		&sst_noaa_ndims,
		sst_noaa_dimsizes,
		NULL,
		NULL,
		NULL,
		1);


	lat_noaa = (float*) NclGetArgValue(
		4,
		11,
		&lat_noaa_ndims,
		lat_noaa_dimsizes,
		NULL,
		NULL,
		NULL,
		1);


	lon_noaa = (float*) NclGetArgValue(
		3,
		11,
		&lon_noaa_ndims,
		lon_noaa_dimsizes,
		NULL,
		NULL,
		NULL,
		1);


	time = (int*) NclGetArgValue(
		2,
		11,
		NULL,
		NULL,
		NULL,
		NULL,
		NULL,
		1);


	nj = (int*) NclGetArgValue(
		1,
		11,
		NULL,
		NULL,
		NULL,
		NULL,
		NULL,
		1);


	ni = (int*) NclGetArgValue(
		0,
		11,
		NULL,
		NULL,
		NULL,
		NULL,
		NULL,
		1);

	if(*ni != (int)lon_noaa_dimsizes[1]) {
		NhlPError(NhlFATAL,NhlEUNKNOWN,"sst_noaa2soda_1hr: dimension size of dimension (1) of lon_noaa must be equal to the value of ni");
		return(NhlFATAL);
	}
	if(*nj != (int)lon_noaa_dimsizes[0]) {
		NhlPError(NhlFATAL,NhlEUNKNOWN,"sst_noaa2soda_1hr: dimension size of dimension (0) of lon_noaa must be equal to the value of nj");
		return(NhlFATAL);
	}
	if(*ni != (int)lat_noaa_dimsizes[1]) {
		NhlPError(NhlFATAL,NhlEUNKNOWN,"sst_noaa2soda_1hr: dimension size of dimension (1) of lat_noaa must be equal to the value of ni");
		return(NhlFATAL);
	}
	if(*nj != (int)lat_noaa_dimsizes[0]) {
		NhlPError(NhlFATAL,NhlEUNKNOWN,"sst_noaa2soda_1hr: dimension size of dimension (0) of lat_noaa must be equal to the value of nj");
		return(NhlFATAL);
	}
	if(*ni != (int)sst_noaa_dimsizes[1]) {
		NhlPError(NhlFATAL,NhlEUNKNOWN,"sst_noaa2soda_1hr: dimension size of dimension (1) of sst_noaa must be equal to the value of ni");
		return(NhlFATAL);
	}
	if(*nj != (int)sst_noaa_dimsizes[0]) {
		NhlPError(NhlFATAL,NhlEUNKNOWN,"sst_noaa2soda_1hr: dimension size of dimension (0) of sst_noaa must be equal to the value of nj");
		return(NhlFATAL);
	}
	if(*ni != (int)sses_bias_dimsizes[1]) {
		NhlPError(NhlFATAL,NhlEUNKNOWN,"sst_noaa2soda_1hr: dimension size of dimension (1) of sses_bias must be equal to the value of ni");
		return(NhlFATAL);
	}
	if(*nj != (int)sses_bias_dimsizes[0]) {
		NhlPError(NhlFATAL,NhlEUNKNOWN,"sst_noaa2soda_1hr: dimension size of dimension (0) of sses_bias must be equal to the value of nj");
		return(NhlFATAL);
	}
	if(*ni != (int)l2p_flags_dimsizes[1]) {
		NhlPError(NhlFATAL,NhlEUNKNOWN,"sst_noaa2soda_1hr: dimension size of dimension (1) of l2p_flags must be equal to the value of ni");
		return(NhlFATAL);
	}
	if(*nj != (int)l2p_flags_dimsizes[0]) {
		NhlPError(NhlFATAL,NhlEUNKNOWN,"sst_noaa2soda_1hr: dimension size of dimension (0) of l2p_flags must be equal to the value of nj");
		return(NhlFATAL);
	}
	if(*ni != (int)quality_level_dimsizes[1]) {
		NhlPError(NhlFATAL,NhlEUNKNOWN,"sst_noaa2soda_1hr: dimension size of dimension (1) of quality_level must be equal to the value of ni");
		return(NhlFATAL);
	}
	if(*nj != (int)quality_level_dimsizes[0]) {
		NhlPError(NhlFATAL,NhlEUNKNOWN,"sst_noaa2soda_1hr: dimension size of dimension (0) of quality_level must be equal to the value of nj");
		return(NhlFATAL);
	}
	NGCALLF(sst_noaa2soda_1hr,SST_NOAA2SODA_1HR)(ni,nj,time,lon_noaa,lat_noaa,sst_noaa,sses_bias,l2p_flags,quality_level,sst_soda,n_sst_soda);

	return(NhlNOERROR);
}
extern void NGCALLF(sst_noaa2soda_10min,SST_NOAA2SODA_10MIN)();


NhlErrorTypes sst_noaa2soda_10min_W( void ) {
	int i;
	int *lon;
	int *lat;
	int *time;
	float *lon_noaa;
	float *lat_noaa;
	float *sst_noaa;
	float *sses_bias;
	short *l2p_flags;
	byte *quality_level;
	float *sst_soda;
	int *n_sst_soda;
	long lon_noaa_dimsizes[NCL_MAX_DIMENSIONS];
	int lon_noaa_ndims;
	long lat_noaa_dimsizes[NCL_MAX_DIMENSIONS];
	int lat_noaa_ndims;
	long sst_noaa_dimsizes[NCL_MAX_DIMENSIONS];
	int sst_noaa_ndims;
	long sses_bias_dimsizes[NCL_MAX_DIMENSIONS];
	int sses_bias_ndims;
	long l2p_flags_dimsizes[NCL_MAX_DIMENSIONS];
	int l2p_flags_ndims;
	long quality_level_dimsizes[NCL_MAX_DIMENSIONS];
	int quality_level_ndims;

	n_sst_soda = (int*) NclGetArgValue(
		10,
		11,
		NULL,
		NULL,
		NULL,
		NULL,
		NULL,
		1);


	sst_soda = (float*) NclGetArgValue(
		9,
		11,
		NULL,
		NULL,
		NULL,
		NULL,
		NULL,
		1);


	quality_level = (byte*) NclGetArgValue(
		8,
		11,
		&quality_level_ndims,
		quality_level_dimsizes,
		NULL,
		NULL,
		NULL,
		1);


	l2p_flags = (short*) NclGetArgValue(
		7,
		11,
		&l2p_flags_ndims,
		l2p_flags_dimsizes,
		NULL,
		NULL,
		NULL,
		1);


	sses_bias = (float*) NclGetArgValue(
		6,
		11,
		&sses_bias_ndims,
		sses_bias_dimsizes,
		NULL,
		NULL,
		NULL,
		1);


	sst_noaa = (float*) NclGetArgValue(
		5,
		11,
		&sst_noaa_ndims,
		sst_noaa_dimsizes,
		NULL,
		NULL,
		NULL,
		1);


	lat_noaa = (float*) NclGetArgValue(
		4,
		11,
		&lat_noaa_ndims,
		lat_noaa_dimsizes,
		NULL,
		NULL,
		NULL,
		1);


	lon_noaa = (float*) NclGetArgValue(
		3,
		11,
		&lon_noaa_ndims,
		lon_noaa_dimsizes,
		NULL,
		NULL,
		NULL,
		1);


	time = (int*) NclGetArgValue(
		2,
		11,
		NULL,
		NULL,
		NULL,
		NULL,
		NULL,
		1);


	lat = (int*) NclGetArgValue(
		1,
		11,
		NULL,
		NULL,
		NULL,
		NULL,
		NULL,
		1);


	lon = (int*) NclGetArgValue(
		0,
		11,
		NULL,
		NULL,
		NULL,
		NULL,
		NULL,
		1);

	if(*lon != (int)lon_noaa_dimsizes[0]) {
		NhlPError(NhlFATAL,NhlEUNKNOWN,"sst_noaa2soda_10min: dimension size of dimension (0) of lon_noaa must be equal to the value of lon");
		return(NhlFATAL);
	}
	if(*lat != (int)lat_noaa_dimsizes[0]) {
		NhlPError(NhlFATAL,NhlEUNKNOWN,"sst_noaa2soda_10min: dimension size of dimension (0) of lat_noaa must be equal to the value of lat");
		return(NhlFATAL);
	}
	if(*lon != (int)sst_noaa_dimsizes[1]) {
		NhlPError(NhlFATAL,NhlEUNKNOWN,"sst_noaa2soda_10min: dimension size of dimension (1) of sst_noaa must be equal to the value of lon");
		return(NhlFATAL);
	}
	if(*lat != (int)sst_noaa_dimsizes[0]) {
		NhlPError(NhlFATAL,NhlEUNKNOWN,"sst_noaa2soda_10min: dimension size of dimension (0) of sst_noaa must be equal to the value of lat");
		return(NhlFATAL);
	}
	if(*lon != (int)sses_bias_dimsizes[1]) {
		NhlPError(NhlFATAL,NhlEUNKNOWN,"sst_noaa2soda_10min: dimension size of dimension (1) of sses_bias must be equal to the value of lon");
		return(NhlFATAL);
	}
	if(*lat != (int)sses_bias_dimsizes[0]) {
		NhlPError(NhlFATAL,NhlEUNKNOWN,"sst_noaa2soda_10min: dimension size of dimension (0) of sses_bias must be equal to the value of lat");
		return(NhlFATAL);
	}
	if(*lon != (int)l2p_flags_dimsizes[1]) {
		NhlPError(NhlFATAL,NhlEUNKNOWN,"sst_noaa2soda_10min: dimension size of dimension (1) of l2p_flags must be equal to the value of lon");
		return(NhlFATAL);
	}
	if(*lat != (int)l2p_flags_dimsizes[0]) {
		NhlPError(NhlFATAL,NhlEUNKNOWN,"sst_noaa2soda_10min: dimension size of dimension (0) of l2p_flags must be equal to the value of lat");
		return(NhlFATAL);
	}
	if(*lon != (int)quality_level_dimsizes[1]) {
		NhlPError(NhlFATAL,NhlEUNKNOWN,"sst_noaa2soda_10min: dimension size of dimension (1) of quality_level must be equal to the value of lon");
		return(NhlFATAL);
	}
	if(*lat != (int)quality_level_dimsizes[0]) {
		NhlPError(NhlFATAL,NhlEUNKNOWN,"sst_noaa2soda_10min: dimension size of dimension (0) of quality_level must be equal to the value of lat");
		return(NhlFATAL);
	}
	NGCALLF(sst_noaa2soda_10min,SST_NOAA2SODA_10MIN)(lon,lat,time,lon_noaa,lat_noaa,sst_noaa,sses_bias,l2p_flags,quality_level,sst_soda,n_sst_soda);

	return(NhlNOERROR);
}
extern void NGCALLF(write_sst_noaa2soda_5day_avg,WRITE_SST_NOAA2SODA_5DAY_AVG)();


NhlErrorTypes write_sst_noaa2soda_5day_avg_W( void ) {
	int i;
	int *ndat;
	float *sst_soda;

	sst_soda = (float*) NclGetArgValue(
		1,
		2,
		NULL,
		NULL,
		NULL,
		NULL,
		NULL,
		1);


	ndat = (int*) NclGetArgValue(
		0,
		2,
		NULL,
		NULL,
		NULL,
		NULL,
		NULL,
		1);

	NGCALLF(write_sst_noaa2soda_5day_avg,WRITE_SST_NOAA2SODA_5DAY_AVG)(ndat,sst_soda);

	return(NhlNOERROR);
}
extern void NGCALLF(write_out_monthlyavg_nc_soda,WRITE_OUT_MONTHLYAVG_NC_SODA)();


NhlErrorTypes write_out_monthlyAvg_nc_SODA_W( void ) {
	int i;
	int *curr_yr;
	int *curr_mon;
	int *i_var;
	NrmQuark *str_var;
	float *var;
	int str_var_len;
	char *str_var_str;

	var = (float*) NclGetArgValue(
		4,
		5,
		NULL,
		NULL,
		NULL,
		NULL,
		NULL,
		1);


	str_var = (NrmQuark*) NclGetArgValue(
		3,
		5,
		NULL,
		NULL,
		NULL,
		NULL,
		NULL,
		1);


	i_var = (int*) NclGetArgValue(
		2,
		5,
		NULL,
		NULL,
		NULL,
		NULL,
		NULL,
		1);


	curr_mon = (int*) NclGetArgValue(
		1,
		5,
		NULL,
		NULL,
		NULL,
		NULL,
		NULL,
		1);


	curr_yr = (int*) NclGetArgValue(
		0,
		5,
		NULL,
		NULL,
		NULL,
		NULL,
		NULL,
		1);

	str_var_len = strlen(NrmQuarkToString(*str_var));
	str_var_str = malloc(str_var_len + 1);
	strcpy(str_var_str,NrmQuarkToString(*str_var));
	NGCALLF(write_out_monthlyavg_nc_soda,WRITE_OUT_MONTHLYAVG_NC_SODA)(curr_yr,curr_mon,i_var,str_var_str,var,str_var_len);

	str_var_str[str_var_len] = '\0';
	*str_var = NrmStringToQuark(str_var_str);
	free(str_var_str);
	return(NhlNOERROR);
}
extern void NGCALLF(write_out_monthlyavg_nc_soda_339,WRITE_OUT_MONTHLYAVG_NC_SODA_339)();


NhlErrorTypes write_out_monthlyAvg_nc_SODA_339_W( void ) {
	int i;
	int *curr_yr;
	int *curr_mon;
	int *i_var;
	NrmQuark *str_var;
	float *var;
	int str_var_len;
	char *str_var_str;

	var = (float*) NclGetArgValue(
		4,
		5,
		NULL,
		NULL,
		NULL,
		NULL,
		NULL,
		1);


	str_var = (NrmQuark*) NclGetArgValue(
		3,
		5,
		NULL,
		NULL,
		NULL,
		NULL,
		NULL,
		1);


	i_var = (int*) NclGetArgValue(
		2,
		5,
		NULL,
		NULL,
		NULL,
		NULL,
		NULL,
		1);


	curr_mon = (int*) NclGetArgValue(
		1,
		5,
		NULL,
		NULL,
		NULL,
		NULL,
		NULL,
		1);


	curr_yr = (int*) NclGetArgValue(
		0,
		5,
		NULL,
		NULL,
		NULL,
		NULL,
		NULL,
		1);

	str_var_len = strlen(NrmQuarkToString(*str_var));
	str_var_str = malloc(str_var_len + 1);
	strcpy(str_var_str,NrmQuarkToString(*str_var));
	NGCALLF(write_out_monthlyavg_nc_soda_339,WRITE_OUT_MONTHLYAVG_NC_SODA_339)(curr_yr,curr_mon,i_var,str_var_str,var,str_var_len);

	str_var_str[str_var_len] = '\0';
	*str_var = NrmStringToQuark(str_var_str);
	free(str_var_str);
	return(NhlNOERROR);
}
extern void NGCALLF(write_out_monthlyavg_isopycn_nc_soda_339,WRITE_OUT_MONTHLYAVG_ISOPYCN_NC_SODA_339)();


NhlErrorTypes write_out_monthlyAvg_isopycn_nc_SODA_339_W( void ) {
	int i;
	int *curr_yr;
	int *curr_mon;
	NrmQuark *str_var;
	float *var;
	int str_var_len;
	char *str_var_str;

	var = (float*) NclGetArgValue(
		3,
		4,
		NULL,
		NULL,
		NULL,
		NULL,
		NULL,
		1);


	str_var = (NrmQuark*) NclGetArgValue(
		2,
		4,
		NULL,
		NULL,
		NULL,
		NULL,
		NULL,
		1);


	curr_mon = (int*) NclGetArgValue(
		1,
		4,
		NULL,
		NULL,
		NULL,
		NULL,
		NULL,
		1);


	curr_yr = (int*) NclGetArgValue(
		0,
		4,
		NULL,
		NULL,
		NULL,
		NULL,
		NULL,
		1);

	str_var_len = strlen(NrmQuarkToString(*str_var));
	str_var_str = malloc(str_var_len + 1);
	strcpy(str_var_str,NrmQuarkToString(*str_var));
	NGCALLF(write_out_monthlyavg_isopycn_nc_soda_339,WRITE_OUT_MONTHLYAVG_ISOPYCN_NC_SODA_339)(curr_yr,curr_mon,str_var_str,var,str_var_len);

	str_var_str[str_var_len] = '\0';
	*str_var = NrmStringToQuark(str_var_str);
	free(str_var_str);
	return(NhlNOERROR);
}
extern void NGCALLF(write_out_monthlyavg_ice_nc_soda_339,WRITE_OUT_MONTHLYAVG_ICE_NC_SODA_339)();


NhlErrorTypes write_out_monthlyAvg_ice_nc_SODA_339_W( void ) {
	int i;
	int *curr_yr;
	int *curr_mon;
	int *i_var;
	NrmQuark *str_var;
	float *var;
	int str_var_len;
	char *str_var_str;

	var = (float*) NclGetArgValue(
		4,
		5,
		NULL,
		NULL,
		NULL,
		NULL,
		NULL,
		1);


	str_var = (NrmQuark*) NclGetArgValue(
		3,
		5,
		NULL,
		NULL,
		NULL,
		NULL,
		NULL,
		1);


	i_var = (int*) NclGetArgValue(
		2,
		5,
		NULL,
		NULL,
		NULL,
		NULL,
		NULL,
		1);


	curr_mon = (int*) NclGetArgValue(
		1,
		5,
		NULL,
		NULL,
		NULL,
		NULL,
		NULL,
		1);


	curr_yr = (int*) NclGetArgValue(
		0,
		5,
		NULL,
		NULL,
		NULL,
		NULL,
		NULL,
		1);

	str_var_len = strlen(NrmQuarkToString(*str_var));
	str_var_str = malloc(str_var_len + 1);
	strcpy(str_var_str,NrmQuarkToString(*str_var));
	NGCALLF(write_out_monthlyavg_ice_nc_soda_339,WRITE_OUT_MONTHLYAVG_ICE_NC_SODA_339)(curr_yr,curr_mon,i_var,str_var_str,var,str_var_len);

	str_var_str[str_var_len] = '\0';
	*str_var = NrmStringToQuark(str_var_str);
	free(str_var_str);
	return(NhlNOERROR);
}
extern void NGCALLF(write_out_monthlyavg_nc_soda_330,WRITE_OUT_MONTHLYAVG_NC_SODA_330)();


NhlErrorTypes write_out_monthlyAvg_nc_SODA_330_W( void ) {
	int i;
	int *curr_yr;
	int *curr_mon;
	int *i_var;
	NrmQuark *str_var;
	float *var;
	int str_var_len;
	char *str_var_str;

	var = (float*) NclGetArgValue(
		4,
		5,
		NULL,
		NULL,
		NULL,
		NULL,
		NULL,
		1);


	str_var = (NrmQuark*) NclGetArgValue(
		3,
		5,
		NULL,
		NULL,
		NULL,
		NULL,
		NULL,
		1);


	i_var = (int*) NclGetArgValue(
		2,
		5,
		NULL,
		NULL,
		NULL,
		NULL,
		NULL,
		1);


	curr_mon = (int*) NclGetArgValue(
		1,
		5,
		NULL,
		NULL,
		NULL,
		NULL,
		NULL,
		1);


	curr_yr = (int*) NclGetArgValue(
		0,
		5,
		NULL,
		NULL,
		NULL,
		NULL,
		NULL,
		1);

	str_var_len = strlen(NrmQuarkToString(*str_var));
	str_var_str = malloc(str_var_len + 1);
	strcpy(str_var_str,NrmQuarkToString(*str_var));
	NGCALLF(write_out_monthlyavg_nc_soda_330,WRITE_OUT_MONTHLYAVG_NC_SODA_330)(curr_yr,curr_mon,i_var,str_var_str,var,str_var_len);

	str_var_str[str_var_len] = '\0';
	*str_var = NrmStringToQuark(str_var_str);
	free(str_var_str);
	return(NhlNOERROR);
}
extern void NGCALLF(write_out_monthlyavg_nc_soda_341,WRITE_OUT_MONTHLYAVG_NC_SODA_341)();


NhlErrorTypes write_out_monthlyAvg_nc_SODA_341_W( void ) {
	int i;
	int *curr_yr;
	int *curr_mon;
	int *i_var;
	NrmQuark *str_var;
	float *var;
	int str_var_len;
	char *str_var_str;

	var = (float*) NclGetArgValue(
		4,
		5,
		NULL,
		NULL,
		NULL,
		NULL,
		NULL,
		1);


	str_var = (NrmQuark*) NclGetArgValue(
		3,
		5,
		NULL,
		NULL,
		NULL,
		NULL,
		NULL,
		1);


	i_var = (int*) NclGetArgValue(
		2,
		5,
		NULL,
		NULL,
		NULL,
		NULL,
		NULL,
		1);


	curr_mon = (int*) NclGetArgValue(
		1,
		5,
		NULL,
		NULL,
		NULL,
		NULL,
		NULL,
		1);


	curr_yr = (int*) NclGetArgValue(
		0,
		5,
		NULL,
		NULL,
		NULL,
		NULL,
		NULL,
		1);

	str_var_len = strlen(NrmQuarkToString(*str_var));
	str_var_str = malloc(str_var_len + 1);
	strcpy(str_var_str,NrmQuarkToString(*str_var));
	NGCALLF(write_out_monthlyavg_nc_soda_341,WRITE_OUT_MONTHLYAVG_NC_SODA_341)(curr_yr,curr_mon,i_var,str_var_str,var,str_var_len);

	str_var_str[str_var_len] = '\0';
	*str_var = NrmStringToQuark(str_var_str);
	free(str_var_str);
	return(NhlNOERROR);
}



void Init(void){
	void *args;
	long dimsizes[NCL_MAX_DIMENSIONS];
	int nargs;


	nargs = 0;
	args = NewArgs(5);
	dimsizes[0] = 50;
	dimsizes[1] = 330;
	dimsizes[2] = 720;
	SetArgTemplate(args,4,"float",3,dimsizes);nargs++;
	dimsizes[0] = 1;
	SetArgTemplate(args,3,"string",1,dimsizes);nargs++;
	dimsizes[0] = 1;
	SetArgTemplate(args,2,"integer",1,dimsizes);nargs++;
	dimsizes[0] = 1;
	SetArgTemplate(args,1,"integer",1,dimsizes);nargs++;
	dimsizes[0] = 1;
	SetArgTemplate(args,0,"integer",1,dimsizes);nargs++;
	NclRegisterProc(write_out_monthlyAvg_nc_SODA_341_W,args,"WRITE_OUT_MONTHLYAVG_NC_SODA_341",nargs);

	NclRegisterProc(write_out_monthlyAvg_nc_SODA_341_W,args,"write_out_monthlyavg_nc_soda_341",nargs);

	nargs = 0;
	args = NewArgs(5);
	dimsizes[0] = 50;
	dimsizes[1] = 330;
	dimsizes[2] = 720;
	SetArgTemplate(args,4,"float",3,dimsizes);nargs++;
	dimsizes[0] = 1;
	SetArgTemplate(args,3,"string",1,dimsizes);nargs++;
	dimsizes[0] = 1;
	SetArgTemplate(args,2,"integer",1,dimsizes);nargs++;
	dimsizes[0] = 1;
	SetArgTemplate(args,1,"integer",1,dimsizes);nargs++;
	dimsizes[0] = 1;
	SetArgTemplate(args,0,"integer",1,dimsizes);nargs++;
	NclRegisterProc(write_out_monthlyAvg_nc_SODA_330_W,args,"WRITE_OUT_MONTHLYAVG_NC_SODA_330",nargs);

	NclRegisterProc(write_out_monthlyAvg_nc_SODA_330_W,args,"write_out_monthlyavg_nc_soda_330",nargs);

	nargs = 0;
	args = NewArgs(5);
	dimsizes[0] = 5;
	dimsizes[1] = 330;
	dimsizes[2] = 720;
	SetArgTemplate(args,4,"float",3,dimsizes);nargs++;
	dimsizes[0] = 1;
	SetArgTemplate(args,3,"string",1,dimsizes);nargs++;
	dimsizes[0] = 1;
	SetArgTemplate(args,2,"integer",1,dimsizes);nargs++;
	dimsizes[0] = 1;
	SetArgTemplate(args,1,"integer",1,dimsizes);nargs++;
	dimsizes[0] = 1;
	SetArgTemplate(args,0,"integer",1,dimsizes);nargs++;
	NclRegisterProc(write_out_monthlyAvg_ice_nc_SODA_339_W,args,"WRITE_OUT_MONTHLYAVG_ICE_NC_SODA_339",nargs);

	NclRegisterProc(write_out_monthlyAvg_ice_nc_SODA_339_W,args,"write_out_monthlyavg_ice_nc_soda_339",nargs);

	nargs = 0;
	args = NewArgs(4);
	dimsizes[0] = 16;
	dimsizes[1] = 330;
	dimsizes[2] = 720;
	SetArgTemplate(args,3,"float",3,dimsizes);nargs++;
	dimsizes[0] = 1;
	SetArgTemplate(args,2,"string",1,dimsizes);nargs++;
	dimsizes[0] = 1;
	SetArgTemplate(args,1,"integer",1,dimsizes);nargs++;
	dimsizes[0] = 1;
	SetArgTemplate(args,0,"integer",1,dimsizes);nargs++;
	NclRegisterProc(write_out_monthlyAvg_isopycn_nc_SODA_339_W,args,"WRITE_OUT_MONTHLYAVG_ISOPYCN_NC_SODA_339",nargs);

	NclRegisterProc(write_out_monthlyAvg_isopycn_nc_SODA_339_W,args,"write_out_monthlyavg_isopycn_nc_soda_339",nargs);

	nargs = 0;
	args = NewArgs(5);
	dimsizes[0] = 50;
	dimsizes[1] = 330;
	dimsizes[2] = 720;
	SetArgTemplate(args,4,"float",3,dimsizes);nargs++;
	dimsizes[0] = 1;
	SetArgTemplate(args,3,"string",1,dimsizes);nargs++;
	dimsizes[0] = 1;
	SetArgTemplate(args,2,"integer",1,dimsizes);nargs++;
	dimsizes[0] = 1;
	SetArgTemplate(args,1,"integer",1,dimsizes);nargs++;
	dimsizes[0] = 1;
	SetArgTemplate(args,0,"integer",1,dimsizes);nargs++;
	NclRegisterProc(write_out_monthlyAvg_nc_SODA_339_W,args,"WRITE_OUT_MONTHLYAVG_NC_SODA_339",nargs);

	NclRegisterProc(write_out_monthlyAvg_nc_SODA_339_W,args,"write_out_monthlyavg_nc_soda_339",nargs);

	nargs = 0;
	args = NewArgs(5);
	dimsizes[0] = 50;
	dimsizes[1] = 330;
	dimsizes[2] = 720;
	SetArgTemplate(args,4,"float",3,dimsizes);nargs++;
	dimsizes[0] = 1;
	SetArgTemplate(args,3,"string",1,dimsizes);nargs++;
	dimsizes[0] = 1;
	SetArgTemplate(args,2,"integer",1,dimsizes);nargs++;
	dimsizes[0] = 1;
	SetArgTemplate(args,1,"integer",1,dimsizes);nargs++;
	dimsizes[0] = 1;
	SetArgTemplate(args,0,"integer",1,dimsizes);nargs++;
	NclRegisterProc(write_out_monthlyAvg_nc_SODA_W,args,"WRITE_OUT_MONTHLYAVG_NC_SODA",nargs);

	NclRegisterProc(write_out_monthlyAvg_nc_SODA_W,args,"write_out_monthlyavg_nc_soda",nargs);

	nargs = 0;
	args = NewArgs(2);
	dimsizes[0] = 180;
	dimsizes[1] = 360;
	SetArgTemplate(args,1,"float",2,dimsizes);nargs++;
	dimsizes[0] = 1;
	SetArgTemplate(args,0,"integer",1,dimsizes);nargs++;
	NclRegisterProc(write_sst_noaa2soda_5day_avg_W,args,"WRITE_SST_NOAA2SODA_5DAY_AVG",nargs);

	NclRegisterProc(write_sst_noaa2soda_5day_avg_W,args,"write_sst_noaa2soda_5day_avg",nargs);

	nargs = 0;
	args = NewArgs(11);
	dimsizes[0] = 180;
	dimsizes[1] = 360;
	SetArgTemplate(args,10,"integer",2,dimsizes);nargs++;
	dimsizes[0] = 180;
	dimsizes[1] = 360;
	SetArgTemplate(args,9,"float",2,dimsizes);nargs++;
	dimsizes[0] = -1;
	dimsizes[1] = -1;
	SetArgTemplate(args,8,"byte",2,dimsizes);nargs++;
	dimsizes[0] = -1;
	dimsizes[1] = -1;
	SetArgTemplate(args,7,"short",2,dimsizes);nargs++;
	dimsizes[0] = -1;
	dimsizes[1] = -1;
	SetArgTemplate(args,6,"float",2,dimsizes);nargs++;
	dimsizes[0] = -1;
	dimsizes[1] = -1;
	SetArgTemplate(args,5,"float",2,dimsizes);nargs++;
	dimsizes[0] = -1;
	SetArgTemplate(args,4,"float",1,dimsizes);nargs++;
	dimsizes[0] = -1;
	SetArgTemplate(args,3,"float",1,dimsizes);nargs++;
	dimsizes[0] = 1;
	SetArgTemplate(args,2,"integer",1,dimsizes);nargs++;
	dimsizes[0] = 1;
	SetArgTemplate(args,1,"integer",1,dimsizes);nargs++;
	dimsizes[0] = 1;
	SetArgTemplate(args,0,"integer",1,dimsizes);nargs++;
	NclRegisterProc(sst_noaa2soda_10min_W,args,"SST_NOAA2SODA_10MIN",nargs);

	NclRegisterProc(sst_noaa2soda_10min_W,args,"sst_noaa2soda_10min",nargs);

	nargs = 0;
	args = NewArgs(11);
	dimsizes[0] = 180;
	dimsizes[1] = 360;
	SetArgTemplate(args,10,"integer",2,dimsizes);nargs++;
	dimsizes[0] = 180;
	dimsizes[1] = 360;
	SetArgTemplate(args,9,"float",2,dimsizes);nargs++;
	dimsizes[0] = -1;
	dimsizes[1] = -1;
	SetArgTemplate(args,8,"byte",2,dimsizes);nargs++;
	dimsizes[0] = -1;
	dimsizes[1] = -1;
	SetArgTemplate(args,7,"short",2,dimsizes);nargs++;
	dimsizes[0] = -1;
	dimsizes[1] = -1;
	SetArgTemplate(args,6,"float",2,dimsizes);nargs++;
	dimsizes[0] = -1;
	dimsizes[1] = -1;
	SetArgTemplate(args,5,"float",2,dimsizes);nargs++;
	dimsizes[0] = -1;
	dimsizes[1] = -1;
	SetArgTemplate(args,4,"float",2,dimsizes);nargs++;
	dimsizes[0] = -1;
	dimsizes[1] = -1;
	SetArgTemplate(args,3,"float",2,dimsizes);nargs++;
	dimsizes[0] = 1;
	SetArgTemplate(args,2,"integer",1,dimsizes);nargs++;
	dimsizes[0] = 1;
	SetArgTemplate(args,1,"integer",1,dimsizes);nargs++;
	dimsizes[0] = 1;
	SetArgTemplate(args,0,"integer",1,dimsizes);nargs++;
	NclRegisterProc(sst_noaa2soda_1hr_W,args,"SST_NOAA2SODA_1HR",nargs);

	NclRegisterProc(sst_noaa2soda_1hr_W,args,"sst_noaa2soda_1hr",nargs);

	nargs = 0;
	args = NewArgs(5);
	dimsizes[0] = 34;
	dimsizes[1] = 71;
	SetArgTemplate(args,4,"integer",2,dimsizes);nargs++;
	dimsizes[0] = 1;
	SetArgTemplate(args,3,"integer",1,dimsizes);nargs++;
	dimsizes[0] = 1;
	SetArgTemplate(args,2,"integer",1,dimsizes);nargs++;
	dimsizes[0] = 1;
	SetArgTemplate(args,1,"integer",1,dimsizes);nargs++;
	dimsizes[0] = 1;
	SetArgTemplate(args,0,"integer",1,dimsizes);nargs++;
	NclRegisterProc(cal_salt_fcs_obs_diff_dist_W,args,"CAL_SALT_FCS_OBS_DIFF_DIST",nargs);

	NclRegisterProc(cal_salt_fcs_obs_diff_dist_W,args,"cal_salt_fcs_obs_diff_dist",nargs);

	nargs = 0;
	args = NewArgs(5);
	dimsizes[0] = 34;
	dimsizes[1] = 71;
	SetArgTemplate(args,4,"integer",2,dimsizes);nargs++;
	dimsizes[0] = 1;
	SetArgTemplate(args,3,"integer",1,dimsizes);nargs++;
	dimsizes[0] = 1;
	SetArgTemplate(args,2,"integer",1,dimsizes);nargs++;
	dimsizes[0] = 1;
	SetArgTemplate(args,1,"integer",1,dimsizes);nargs++;
	dimsizes[0] = 1;
	SetArgTemplate(args,0,"integer",1,dimsizes);nargs++;
	NclRegisterProc(cal_temp_fcs_obs_diff_dist_W,args,"CAL_TEMP_FCS_OBS_DIFF_DIST",nargs);

	NclRegisterProc(cal_temp_fcs_obs_diff_dist_W,args,"cal_temp_fcs_obs_diff_dist",nargs);

	nargs = 0;
	args = NewArgs(18);
	dimsizes[0] = 1;
	SetArgTemplate(args,17,"integer",1,dimsizes);nargs++;
	dimsizes[0] = 1;
	SetArgTemplate(args,16,"integer",1,dimsizes);nargs++;
	dimsizes[0] = 1;
	SetArgTemplate(args,15,"integer",1,dimsizes);nargs++;
	dimsizes[0] = 1;
	dimsizes[1] = -1;
	dimsizes[2] = -1;
	SetArgTemplate(args,14,"float",3,dimsizes);nargs++;
	dimsizes[0] = 1;
	dimsizes[1] = -1;
	dimsizes[2] = -1;
	SetArgTemplate(args,13,"float",3,dimsizes);nargs++;
	dimsizes[0] = 1;
	dimsizes[1] = -1;
	dimsizes[2] = -1;
	dimsizes[3] = -1;
	SetArgTemplate(args,12,"float",4,dimsizes);nargs++;
	dimsizes[0] = 1;
	dimsizes[1] = -1;
	dimsizes[2] = -1;
	dimsizes[3] = -1;
	SetArgTemplate(args,11,"float",4,dimsizes);nargs++;
	dimsizes[0] = 1;
	dimsizes[1] = -1;
	dimsizes[2] = -1;
	SetArgTemplate(args,10,"float",3,dimsizes);nargs++;
	dimsizes[0] = 1;
	dimsizes[1] = -1;
	dimsizes[2] = -1;
	SetArgTemplate(args,9,"float",3,dimsizes);nargs++;
	dimsizes[0] = 1;
	dimsizes[1] = -1;
	dimsizes[2] = -1;
	SetArgTemplate(args,8,"float",3,dimsizes);nargs++;
	dimsizes[0] = 1;
	dimsizes[1] = -1;
	dimsizes[2] = -1;
	SetArgTemplate(args,7,"float",3,dimsizes);nargs++;
	dimsizes[0] = 1;
	dimsizes[1] = -1;
	dimsizes[2] = -1;
	SetArgTemplate(args,6,"float",3,dimsizes);nargs++;
	dimsizes[0] = 1;
	dimsizes[1] = -1;
	dimsizes[2] = -1;
	dimsizes[3] = -1;
	SetArgTemplate(args,5,"float",4,dimsizes);nargs++;
	dimsizes[0] = 1;
	dimsizes[1] = -1;
	dimsizes[2] = -1;
	dimsizes[3] = -1;
	SetArgTemplate(args,4,"float",4,dimsizes);nargs++;
	dimsizes[0] = 1;
	dimsizes[1] = -1;
	dimsizes[2] = -1;
	dimsizes[3] = -1;
	SetArgTemplate(args,3,"float",4,dimsizes);nargs++;
	dimsizes[0] = 1;
	dimsizes[1] = -1;
	dimsizes[2] = -1;
	dimsizes[3] = -1;
	SetArgTemplate(args,2,"float",4,dimsizes);nargs++;
	dimsizes[0] = 1;
	SetArgTemplate(args,1,"integer",1,dimsizes);nargs++;
	dimsizes[0] = 1;
	SetArgTemplate(args,0,"integer",1,dimsizes);nargs++;
	NclRegisterProc(write_out_5dayAvg2monthlyAvg_SODA_MOM_1979_1981_W,args,"WRITE_OUT_5DAYAVG2MONTHLYAVG_SODA_MOM_1979_1981",nargs);

	NclRegisterProc(write_out_5dayAvg2monthlyAvg_SODA_MOM_1979_1981_W,args,"write_out_5dayavg2monthlyavg_soda_mom_1979_1981",nargs);

	nargs = 0;
	args = NewArgs(7);
	dimsizes[0] = 1;
	SetArgTemplate(args,6,"integer",1,dimsizes);nargs++;
	dimsizes[0] = 1;
	SetArgTemplate(args,5,"integer",1,dimsizes);nargs++;
	dimsizes[0] = 1;
	SetArgTemplate(args,4,"integer",1,dimsizes);nargs++;
	dimsizes[0] = 1;
	dimsizes[1] = -1;
	dimsizes[2] = -1;
	SetArgTemplate(args,3,"float",3,dimsizes);nargs++;
	dimsizes[0] = 1;
	dimsizes[1] = -1;
	dimsizes[2] = -1;
	dimsizes[3] = -1;
	SetArgTemplate(args,2,"float",4,dimsizes);nargs++;
	dimsizes[0] = 1;
	SetArgTemplate(args,1,"integer",1,dimsizes);nargs++;
	dimsizes[0] = 1;
	SetArgTemplate(args,0,"integer",1,dimsizes);nargs++;
	NclRegisterProc(write_out_5dayAvg2monthlyAvg_SODA_MOM_ice_W,args,"WRITE_OUT_5DAYAVG2MONTHLYAVG_SODA_MOM_ICE",nargs);

	NclRegisterProc(write_out_5dayAvg2monthlyAvg_SODA_MOM_ice_W,args,"write_out_5dayavg2monthlyavg_soda_mom_ice",nargs);

	nargs = 0;
	args = NewArgs(13);
	dimsizes[0] = 1;
	SetArgTemplate(args,12,"integer",1,dimsizes);nargs++;
	dimsizes[0] = 1;
	SetArgTemplate(args,11,"integer",1,dimsizes);nargs++;
	dimsizes[0] = 1;
	SetArgTemplate(args,10,"integer",1,dimsizes);nargs++;
	dimsizes[0] = 1;
	dimsizes[1] = -1;
	dimsizes[2] = -1;
	SetArgTemplate(args,9,"float",3,dimsizes);nargs++;
	dimsizes[0] = 1;
	dimsizes[1] = -1;
	dimsizes[2] = -1;
	SetArgTemplate(args,8,"float",3,dimsizes);nargs++;
	dimsizes[0] = 1;
	dimsizes[1] = -1;
	dimsizes[2] = -1;
	dimsizes[3] = -1;
	SetArgTemplate(args,7,"float",4,dimsizes);nargs++;
	dimsizes[0] = 1;
	dimsizes[1] = -1;
	dimsizes[2] = -1;
	SetArgTemplate(args,6,"float",3,dimsizes);nargs++;
	dimsizes[0] = 1;
	dimsizes[1] = -1;
	dimsizes[2] = -1;
	dimsizes[3] = -1;
	SetArgTemplate(args,5,"float",4,dimsizes);nargs++;
	dimsizes[0] = 1;
	dimsizes[1] = -1;
	dimsizes[2] = -1;
	dimsizes[3] = -1;
	SetArgTemplate(args,4,"float",4,dimsizes);nargs++;
	dimsizes[0] = 1;
	dimsizes[1] = -1;
	dimsizes[2] = -1;
	dimsizes[3] = -1;
	SetArgTemplate(args,3,"float",4,dimsizes);nargs++;
	dimsizes[0] = 1;
	dimsizes[1] = -1;
	dimsizes[2] = -1;
	dimsizes[3] = -1;
	SetArgTemplate(args,2,"float",4,dimsizes);nargs++;
	dimsizes[0] = 1;
	SetArgTemplate(args,1,"integer",1,dimsizes);nargs++;
	dimsizes[0] = 1;
	SetArgTemplate(args,0,"integer",1,dimsizes);nargs++;
	NclRegisterProc(write_out_5dayAvg2monthlyAvg_SODA_MOM_W,args,"WRITE_OUT_5DAYAVG2MONTHLYAVG_SODA_MOM",nargs);

	NclRegisterProc(write_out_5dayAvg2monthlyAvg_SODA_MOM_W,args,"write_out_5dayavg2monthlyavg_soda_mom",nargs);

	nargs = 0;
	args = NewArgs(11);
	dimsizes[0] = -1;
	SetArgTemplate(args,10,"float",1,dimsizes);nargs++;
	dimsizes[0] = 1;
	SetArgTemplate(args,9,"integer",1,dimsizes);nargs++;
	dimsizes[0] = 1;
	SetArgTemplate(args,8,"integer",1,dimsizes);nargs++;
	dimsizes[0] = 1;
	SetArgTemplate(args,7,"logical",1,dimsizes);nargs++;
	dimsizes[0] = 1;
	SetArgTemplate(args,6,"logical",1,dimsizes);nargs++;
	dimsizes[0] = 1;
	SetArgTemplate(args,5,"logical",1,dimsizes);nargs++;
	dimsizes[0] = -1;
	dimsizes[1] = -1;
	SetArgTemplate(args,4,"float",2,dimsizes);nargs++;
	dimsizes[0] = -1;
	dimsizes[1] = -1;
	SetArgTemplate(args,3,"float",2,dimsizes);nargs++;
	dimsizes[0] = 1;
	SetArgTemplate(args,2,"integer",1,dimsizes);nargs++;
	dimsizes[0] = -1;
	dimsizes[1] = -1;
	SetArgTemplate(args,1,"float",2,dimsizes);nargs++;
	dimsizes[0] = -1;
	dimsizes[1] = -1;
	SetArgTemplate(args,0,"float",2,dimsizes);nargs++;
	NclRegisterProc(compute_pr_ets_us_cwrf_W,args,"COMPUTE_PR_ETS_US_CWRF",nargs);

	NclRegisterProc(compute_pr_ets_us_cwrf_W,args,"compute_pr_ets_us_cwrf",nargs);

	nargs = 0;
	args = NewArgs(4);
	dimsizes[0] = 1;
	SetArgTemplate(args,3,"integer",1,dimsizes);nargs++;
	dimsizes[0] = 1;
	SetArgTemplate(args,2,"float",1,dimsizes);nargs++;
	dimsizes[0] = 1;
	SetArgTemplate(args,1,"integer",1,dimsizes);nargs++;
	dimsizes[0] = -1;
	SetArgTemplate(args,0,"float",1,dimsizes);nargs++;
	NclRegisterProc(search_index_binary_W,args,"SEARCH_INDEX_BINARY",nargs);

	NclRegisterProc(search_index_binary_W,args,"search_index_binary",nargs);

	nargs = 0;
	args = NewArgs(2);
	dimsizes[0] = 1;
	SetArgTemplate(args,1,"integer",1,dimsizes);nargs++;
	dimsizes[0] = -1;
	SetArgTemplate(args,0,"float",1,dimsizes);nargs++;
	NclRegisterProc(qsort_W,args,"QSORT",nargs);

	NclRegisterProc(qsort_W,args,"qsort",nargs);

}
