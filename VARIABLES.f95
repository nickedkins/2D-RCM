MODULE VARIABLES

    integer,parameter				::MAXLAYM=203,MXMOLM=38,NBANDSM=16,MAXNCOLS=99,MAXNCLOUDCOLS=10,maxextraclds=2
    real,dimension(MAXLAYM)			::pavelm,tavelm,wbrodlm
    real,dimension(0:MAXLAYM)		::altzm,tzm,pzm,totuflum,totdflum,fnetm,htrm=0.,logpzm
    real,dimension(MXMOLM,MAXLAYM)	::wklm
    real,dimension(NBANDSM)			::semism
    real,parameter					::gravity=9.81 !m s-2 on Earth !NJE
    real,parameter					::mco2=44.01,mo2=32.0,mn2=28.0,mar=40.0 !molecular weight of co2
    real,parameter					::mmwn2 = 28.0134*1e-3,mmwo2 = 31.9988*1e-3,mmwco2 = 44.01*1e-3,mmwh2o=18.01528*1e-3
    real,parameter					::cpn2 = 1040.0,cpo2 = 919.0,cvn2 = 743.0,cvo2 = 659.0
    real,parameter					::avogadro=6.02214e23
    real,parameter					::boltzmann=1.3806488e-23 !J/K
    real,parameter					::universal=avogadro*boltzmann
    real,parameter					::solar_constant=1362.0
    real,parameter					::Liv=2836631.0,eps=0.622,cplh=1003.0,hv = 2501000.0
    real,parameter :: r_earth = 6371e3

    integer	::nlayersm,ioutm,iatmm,ixsectm,iscatm,numangsm,icldm,iemissm,nmolm,ireflecm
    real	::tboundm
    real	::scalep,surfacep
    integer :: cloudcol
    real :: cloudcolfrac=0.,cloudcolalt=0.,cloudcoltau=0.
!    real,dimension(maxextraclds) :: extra_cld_taus,extra_cld_alts,extra_cld_fracs,extra_cld_latcols,extra_cld_cldcols
    REAL :: extra_cld_tau,extra_cld_alt,extra_cld_frac,extra_cld_latcol,extra_cld_cldcol
    integer :: extra_cloudindex
    
    real,dimension(MAXNCLOUDCOLS,MAXNCOLS ) :: latcol_cloudcol_olrs
    
    real :: min_press, inversion_strength
    integer :: inversion_col
    real :: ks,kl,eta,phim,planet_radius,planet_rotation,twarm,tcold,dqstar_by_dt,Lv
    
    real,dimension(MAXNCOLS) :: ddry,lambda,d_vl,x_lats,tot_transp
    
    real :: avl,bvl,cvl,dvl,evl,fvl,gvl

    real :: tcels,t1_vl=288.0,t2_vl=268.0,delta_pv_star,delta_x_lats,delta_temp
    real, dimension(0:maxncols) :: x_edge,delta_x_edge,delta_T_edge,meridtransp_edge,delta_y_edge
    real,dimension(MAXNCOLS) :: delta_meridtransp_edge
    integer :: adj1,adj2,adj3,edge
    

    ! ----------------------------------- NJE trying something here

    real :: massatmo,totmolec,undrelax,lapsecrit
    real,dimension(MAXNCOLS) :: lapsecritcols,height_lccols,od_lccols,frac_lccols,surf_rhcols,R_gcols
    real    :: mixo2
    integer ::i,j,k,l,m,n,x,timesteps,lev,ipe,irdcld,pb_new,loopvar
    real,dimension(maxlaym) :: mixo3,mixh2o,altlaym,mixco2
    real,dimension(maxlaym) :: cmo3,totcmo3,solabso3vis,solabso3uv,solabso3tot,laysolabso3tot,htro3
    real,dimension(maxlaym) :: mperlayr
    real, dimension(maxlaym) :: totcmh2o,solabsh2o,laysolabsh2o,htrh2o,cmh2o
    real :: tot_albedo,sol_inc
    real,dimension(maxlaym) :: rel_hum,es
    real :: rmin
    real :: co2_inv, air_inv,n2_inv,o2_inv
    real,dimension(maxlaym) :: co2_molecs, cptot,cpco2
    real :: massatmo_air,totmolec_air,totmolec_co2
    real, dimension(maxlaym) :: mperlayr_air,mperlayr_co2
    real :: toa_precision,adj_speed,swo3,swh2o,maxhtr,adjspeedfactor, adj_freq
    real,dimension(maxncols) :: fixed_trop
    integer :: OLR_layer,pb,fixed_sw_on,fp,convecttype,stepssinceadj
    integer,dimension(MAXNCOLS) :: conv_trop_ind
    real :: fixed_sw,surf_rh
    real,dimension(MAXNCOLS) :: abs_sw
    real,dimension(maxlaym,MAXNCOLS) :: conv
    real,dimension(maxlaym) :: u_lw
    real :: H
    real,dimension(17,maxlaym) :: band_flux_up, band_flux_down
    ! real, parameter :: PI = 3.14159265359
    real :: startsecs, startmins, endsecs, endmins, starthours, endhours, totsecs, secsperloop
    integer :: fixed_trop_ind
    real :: pico2, pin2, pio2, massatmo_n2, massatmo_co2, massatmo_o2
    real :: molec_co2, molec_n2, molec_o2, vol_mixn2, vol_mixo2, mass_mixn2, mass_mixo2, mmwtot, mass_mixco2, vol_mixco2
    real, dimension(maxlaym) :: rspecific_co2, cvco2, cvtot, rsp_tot, malr
    real :: htransp, currentmaxhtr
    real, dimension(maxlaym) :: newur
    integer :: detailprint
    real,dimension(0:maxlaym) :: kappa,theta,montgomery,exner,geopotential
    real,dimension(0:maxlaym,MAXNCOLS) :: tzmcols,altzmcols,totuflumcols,totdflumcols,htrmcols,pzmcols
    real,dimension(maxlaym,MAXNCOLS) :: tavelmcols,htrh2ocols,htro3cols,wklm1cols,wklm2cols,wklm3cols,wbrodlmcols,abspncols,&
    a_oz_lcols,altlaymcols,pavelmcols,tboundmcols,tau_cldcols,fracscols,rel_hum_cols
    real,dimension(maxncols) :: abs_surf_lhcols
    integer :: col=0, transpcalled, stepssinceboxadj
    real,dimension(MAXNCOLS):: olrcols,insolcols,boxnetradflux,boxnettotflux,meridtransp,currentmaxhtrcols,tempchanges,&
    solar_constants,zencols,boxnetradflux_prev,boxnettotflux_prev
    real,dimension(0:MAXNCOLS) :: tair_lowest_edges
    real :: mtranspfac, boxnetfluxfac
    real :: tglobsum, tglobmean
    real,dimension(MAXNCOLS) :: solar_constantcols
    real :: startlat,endlat
    real,dimension(0:MAXNCOLS) :: latbounds
    real,dimension(MAXNCOLS) :: boxlats
    integer :: day,hour
    real :: Hrad,Xrad,Yrad,hourang,declin,cossums
    real,dimension(MAXNCOLS,365,24) :: insol,zen
    character(len=1024) :: qfn,o3fn,ccfn,clwcfn,ciwcfn,ccfracsfn,cctausfn,ccaltsfn,t_fn !filenames for input distbns
    real,dimension(maxlaym) :: clwc,lwp,ciwc,iwp
    integer :: pertlay, pertcol
    character(len=100) :: ttsfile
    character(len=4):: pertvar
    ! Lacis and Hansen SW variables
    real*4 :: R_g, mu_0,mag,Ra,Ta,Rb,Tb,Rb_star,Tb_star,Ra_star,Ta_star
    real*4,dimension(8) :: k_n,dk,Nlh
    real,dimension(maxlaym) :: wl,delta_z,zlh,rho_w,plh,tlh,tau_cld,abspn,delta_p,htrlh,temp_plh,zlh1,tau_cld_sw
    real*8,dimension(maxlaym,8) :: tau,ssa,ulh,temp,trlh,Translh,Rlh,Rl_star,Tl_star,R_1l,T_1l,&
    &R_1l_star,T_1l_star
    real,dimension(0:maxlaym,8) :: R_lg,Ul,Dl,A_1l,Al
    real :: p_0,t_0,asym,refl_lay_z
    integer :: cld,refl_lay_ind,Llh
    real :: height_hc,frac_hc,od_hc,height_mc,frac_mc,od_mc,height_lc,frac_lc,od_lc,days,mixco2_in,cldfracnje,normalizer
    integer :: ncloudcols,cctemp,ncols
    real,dimension(maxncloudcols) :: ccfracs,cctaus,ccalts
    ! real :: cloudcoltau,cloudcolp
    real,dimension(0:maxlaym) :: htrmwghtd,totuflumwghtd,totdflumwghtd
    real,dimension(maxlaym) :: htrlhwghtd,htro3_lhwghtd,abspnwghtd,A_oz_lwghtd
    real :: tot_sol_abs_lhwghtd
    real :: currentmaxhtr1, currentmaxhtr2
    integer :: htrmaxloc
    real,dimension(maxlaym,2) :: htrm_store,undrelax_store,tavelm_store
    real,dimension(maxlaym) :: temp_tendency
    real,dimension(0:maxlaym) :: sigma
    real :: t_min
    real,dimension(0:maxlaym) :: htrm_over_newur
    real :: ur_htr,ur_seb
    real,dimension(maxncols) :: ur_toafnet
    integer :: couple_tgta
    real,dimension(maxncols) :: abs_h2o_cols,abs_o3_cols,abs_surf_cols
    real :: xglobsum,xglobmean,netradflux_globmean,meridtransp_globmean,nettotflux_globmean
    real :: olr_globmean,abs_sw_globmean,abs_h2o_globmean,abs_o3_globmean,abs_surf_globmean,seb_globmean
    real :: gas_amt_fac_h2o,gas_amt_fac_co2,gas_amt_fac_o3,gas_amt_p_high_h2o,gas_amt_p_low_h2o,gas_amt_p_high_co2,&
    gas_amt_p_low_co2,gas_amt_p_high_o3,gas_amt_p_low_o3
    integer :: gas_amt_pert_h2o,gas_amt_pert_co2,gas_amt_pert_o3,mixco2_prescribed_on,steps_before_toa_adj,cloudloctype
    real :: psurf_override,mixco2_prescribed,a_green,b_green,c_green,H_green
    real :: h_scale,f_cor,beta,gamma_d
    real, dimension(maxncols) :: d_mid,d_trop
    integer :: lapse_type,h2o_sb,h2o_for,h2o_source,mtransp_type,steps_before_first_eqbcheck
    real :: ur_mt,gas_addmolec_h2o,gas_addmolec_co2,gas_addmolec_o3,max_rh,omega_rh


    !Lacis and Hansen Ozone variables
    real :: a,b,c,tau_c,rbar_a,rbarbarstar_a,rbar,rbarbar_a
    integer :: hc
    real,dimension(maxlaym) ::x_o3,cld_fracs,u,A_oz_vis_x,xstar,A_oz_vis_xstar,a_oz_uv_x,a_oz_uv_xstar,a_oz_x,a_oz_xstar,a_oz_l,&
    &htro3_lh
    real,dimension(maxlaym,5) :: htro3_lh_tot,htrlh_tot,abspn_tot
    real, dimension(4) :: fracweight,abs_surf_tot
    real, dimension(8) :: pk

    !Surface
    real :: w_t,ylh,A_wv,Ag1,Ag2,Rbar_r,abs_h2o,abs_o3,abs_surf,tot_sol_abs_lh,abs_surf_lh
    real :: abs_surf_lhwghtd,abs_h2owghtd,abs_o3wghtd
    real :: sebfac,seb,lhf,shf,bowen,c_drag,meanwind,density
    real, dimension(maxncols) :: sebcols
    integer :: sfc_heating,playtype,mtranspon,surf_emiss_on

    !To keep the terminal open
    character*1 KeyBuf

    !Declaration of timing variables
    character(8)  :: date
    character(10) :: time
    character(5)  :: zone
    integer,dimension(8) :: values=0.
    character*20 my_output_file
    character*89 outputfileloc
    real :: time1,time2

    ! ---------------------------------------- NJE end trying

    ! NJE printtime variables

    ! character(8)  :: date
    ! character(10) :: time
    ! character(5)  :: zone
    ! integer,dimension(8) :: values
    character(len=6) :: whencalled

    ! ! NJE calcswhtr variables

    ! real,dimension(0:MAXLAYM) :: altzm,tzm,pzm
    ! real,dimension(MXMOLM,MAXLAYM)  ::wklm
    ! real :: a,b,c,tau_c,rbar_a,rbarbarstar_a,rbar,rbarbar_a
    ! integer :: hc
    ! real,dimension(maxlaym) ::x_o3,u,A_oz_vis_x,xstar,A_oz_vis_xstar,a_oz_uv_x,a_oz_uv_xstar,a_oz_x,&
    ! a_oz_xstar,a_oz_l,htro3_lh
    ! real,dimension(maxlaym) :: htro3_lh_tot,htrlh_tot,abspn_tot
    ! real :: abs_surf_tot
    ! real, dimension(8) :: pk
    ! !Surface
    ! real :: w_t,ylh,A_wv,Ag1,Ag2,Rbar_r,abs_h2o,abs_o3,abs_surf,tot_sol_abs_lh,abs_surf_lh
    ! real,dimension(maxlaym) :: mixh2o
    ! real :: mmwtot
    ! real*4 :: R_g, mu_0,mag,Ra,Ta,Rb,Tb,Rb_star,Tb_star,Ra_star,Ta_star
    ! real*4,dimension(8) :: k_n,dk,Nlh
    ! real,dimension(maxlaym) :: wl,delta_z,zlh,rho_w,plh,tlh,tau_cld,abspn,delta_p,htrlh,temp_plh,zlh1
    ! real*8,dimension(maxlaym,8) :: tau,ssa,ulh,temp,trlh,Translh,Rlh,Rl_star,Tl_star,R_1l,T_1l,&
    ! &R_1l_star,T_1l_star
    ! real,dimension(0:maxlaym,8) :: R_lg,Ul,Dl,A_1l,Al
    ! real :: p_0,t_0,asym,refl_lay_z
    ! integer :: cld,refl_lay_ind,Llh
    ! real :: days,mixco2_in
    ! ! real :: cloudcoltau,cloudp
    ! integer :: cloudindex

    ! ! NJE levconvect variables
    ! real,dimension(0:MAXLAYM) :: altzm,tzm
    ! real,dimension(maxlaym) :: malr,altlaym
    ! real :: lapsecrit,fixed_trop
    ! integer :: convecttype,col,i
    ! real,dimension(maxlaym,MAXNCOLS) :: conv

END MODULE VARIABLES
