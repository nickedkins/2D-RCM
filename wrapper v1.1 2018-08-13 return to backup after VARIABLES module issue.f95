! this program written by Roger Davies 27/8/2012ff
! designed initially to provide a wrapper round prrtm, to preserve ID of prttm
! This version tidied up and using subroutines where possible 2018-05-09 NJE
subroutine wrapper
    use VARIABLES
    use MYSUBS
    implicit none

    call printtime

    ! Initialise some variables to zero (unsure why but this is sometimes needed to avoid errors)
    adj1 = 0
    adj2 = 0
    adj3 = 0
    mixh2o = 0.
    rel_hum = 0.
    rho_w = 0.
    wl = 0.
    abs_h2o = 0.
    abs_surf = 0.
    R_g = 0.
    mu_0 = 0.
    mag = 0.
    Ra = 0.
    Ta = 0.
    Rb = 0.
    Tb = 0.
    Rb_star = 0.
    Tb_star = 0.
    Ra_star = 0.
    Ta_star = 0.
    k_n = 0.
    dk = 0.
    Nlh = 0.
    wl = 0.
    delta_z = 0.
    zlh = 0.
    rho_w = 0.
    plh = 0.
    tlh = 0.
    tau_cld = 0.
    abspn = 0.
    delta_p = 0.
    htrlh = 0.
    htro3_lh = 0.
    temp_plh = 0.
    tau = 0.
    ssa = 0.
    ulh = 0.
    temp = 0.
    trlh = 0.
    Translh = 0.
    Rlh = 0.
    Rl_star = 0.
    Tl_star = 0.
    R_1l = 0.
    T_1l = 0. 
    R_1l_star = 0.
    T_1l_star = 0.
    R_lg = 0.
    Ul = 0.
    Dl = 0.
    A_1l = 0.
    Al = 0.
    p_0=0
    t_0 = 0.
    asym = 0.
    refl_lay_z = 0.
    cld=0
    refl_lay_ind=0
    Llh=0

    Lv = 2.25e6 !latent heat of vaporisation water [J/kg]

    !Set the output file name to the current date and time at the start of the run
    call date_and_time(date,time,zone,values)
    call date_and_time(DATE=date,ZONE=zone)
    call date_and_time(TIME=time)
    call date_and_time(VALUES=values)

    write(my_output_file,1100) values(1),values(2),values(3),&
    &values(5),values(6),values(7)

    outputfileloc = trim('_Raw Output Data/')&
    &//trim(my_output_file)    

    OPEN (50,FILE=outputfileloc,FORM='FORMATTED')

    ! Read input parameters from the file 'Earth RCM Parameters' with is created by the Python script 'Call Fortran from Python.py'
    open(73,file=trim('Earth RCM Parameters'),form='formatted')

    read(73,*) ncols !number of latitude columns
    read(73,*) ncloudcols !number of independent cloud columns within one latitude column
    read(73,*) !tot_albedo !NJE delete this
    read(73,*) !solar_constant    
    read(73,*) tboundmcols(:ncols) !temperature of lower boundary (surface)
    read(73,*) lapsecritcols(:ncols) !critical lapse rate (may be overridden by convecttype)
    read(73,*) days !number of model days
    read(73,*) mixco2_in !CO2 mixing ratio (volumetric)
    read(73,*) undrelax !under-relaxation constant
    read(73,*) icldm !clouds on/off
    read(73,*) rmin !minimum water vapour mixing ratio
    read(73,*) height_hc !temperature of the highest cloud
    read(73,*) frac_hc !fraction of the highest cloud
    read(73,*) od_hc !optical depth of the highest cloud
    read(73,*) height_mc !temperature of the middle cloud
    read(73,*) frac_mc !fraction of the middle cloud
    read(73,*) od_mc !optical depth of the middle cloud
    read(73,*) height_lccols(:ncols) !temperature of the lowest cloud
    read(73,*) frac_lccols(:ncols) !fraction of the lowest cloud
    read(73,*) od_lccols(:ncols) !optical depth of the lowest cloud
    read(73,*) toa_precision !eqb is reached when |OLR - abs SW| < toa_precision
    read(73,*) R_gcols(:ncols) !surface albedo
    read(73,*) fixed_trop(:ncols) !height to which convection is forced to go
    read(73,*) OLR_layer !layer regarded as top (NJE remove this)
    read(73,*) adj_speed !NJE improve this (make this the initial guess, then use Newtonian)
    read(73,*) !cs !cloudsurface on/off
    read(73,*) pb !NJE remove this
    read(73,*) fixed_sw_on !use a fixed value for total abs SW instead of that calculated by Lacis and Hansen
    read(73,*) fixed_sw !the value of the total abs SW when fixed
    read(73,*) fp !use a fixed temperature profile
    read(73,*) surf_rhcols(:ncols) !surface relative humidity for Manabe and Wetherald
    read(73,*) !ps1 !set the surface pressure to 1 bar
    read(73,*) adj_freq !NJE remove this
    read(73,*) convecttype !dry adiabatic lapse rate on/off
    read(73,*) pb_new !new method of turning pressure broadening off    
    read(73,*) swo3 !include O3 SW heating rates
    read(73,*) swh2o !include H2O SW heating rates
    read(73,*) nlayersm !number of model layers
    read(73,*) maxhtr !maximum acceptable heating rate in stratosphere for RCE
    read(73,*) adjspeedfactor !factor applied to adj_speed to decrease time to eqb
    read(73,*) !tuf !factor applied to undrelax to speed things up when near eqb
    read(73,*) pico2 !co2 inventory in bar
    read(73,*) pin2 !n2 inventory in bar
    read(73,*) pio2 !o2 inventory in bar
    read(73,*) htransp !factor reducing lapse rate to account for horizontal transport
    read(73,*) ipe !inert pressure effect (1 means psurf = piair + pico2, 0 mean s psurf = piair)
    read(73,*) detailprint !1: print heating rates, 0: print only start and end times
    read(73,*) mtranspfac !for simple diffusion meridional heat transport, multiply difference between box T and global mean T by this factor
    read(73,*) boxnetfluxfac !obsolete
    read(73,*) pertlay !obsolete
    read(73,*) pertcol !obsolete
    read(73,*) boxlats(:ncols) !latitude at centre of each latitudinal box
    read(73,*) inversion_strength !for manual inversion
    read(73,*) inversion_col !lat column of manual inversion
    read(73,*) twarm !temperature of warm edge of heat engine in Vladilo's heat transport parameterisation
    read(73,*) tcold !temperature of warm edge of heat engine in Vladilo's heat transport parameterisation
    read(73,*) phim !latitude of maximum kinetic energy in Vladilo's heat transport parameterisation
    read(73,*) ks !correlation coefficient in Vladilo's heat transport parameterisation
    read(73,*) kl !correlation coefficient in Vladilo's heat transport parameterisation
    read(73,*) eta !fraction of generated kinetic energy used by heat transporting eddies
    read(73,*) planet_radius !planetary radius (km)
    read(73,*) planet_rotation !planetary rotation rate (rad/s)
    read(73,*) latbounds(0:ncols) !latitudes at the edges of the latitude boxes
    read(73,*) t_min !minimum atmospheric temperature allowed (K)
    read(73,*) sebfac !obsolete
    read(73,*) sfc_heating !apply heating rate from surface energy budget to surface temperature? 0=no, 1=yes
    read(73,*) playtype !pressure layering type. 0=equally spaced in p, 1=sigma formulation
    read(73,*) ur_htr !under-relaxation constant applied to atmospheric heating rates
    read(73,*) ur_toafnet(:ncols) !under-relaxation constant applied to TOA net flux (when changing surface temperature to balance TOA net flux)
    read(73,*) ur_seb !under-relaxation constant applied to surface energy budget if sfc_heating=1
    read(73,*) couple_tgta !0:lowest atmospheric temperature decoupled from surface temperature, 1:lowest atmospheric temperature equals surface temperature
    read(73,*) mtranspon !1:meridional heat transport on
    read(73,*) min_press !lowest pressure in model, at the top of the highest layer (hPa)
    read(73,*) gas_amt_fac_h2o !factor applied to h2o amount within a specified pressure range, for perturbation experiments (ditto below)
    read(73,*) gas_amt_fac_co2
    read(73,*) gas_amt_fac_o3
    read(73,*) gas_amt_p_high_h2o !high end of pressure range over which factor is applied for h2o
    read(73,*) gas_amt_p_low_h2o !low end of pressure range over which factor is applied for h2o
    read(73,*) gas_amt_p_high_co2
    read(73,*) gas_amt_p_low_co2
    read(73,*) gas_amt_p_high_o3
    read(73,*) gas_amt_p_low_o3
    read(73,*) gas_amt_pert_h2o !obsolete
    read(73,*) gas_amt_pert_co2 !obsolete
    read(73,*) gas_amt_pert_o3 !obsolete
    read(73,*) psurf_override !manual override of surface pressure calculated from sum of gas inventories (set to 0 for no override)
    read(73,*) mixco2_prescribed_on !1:specify CO2 mixing ratio instead of inventory
    read(73,*) mixco2_prescribed !prescribed CO2 mixing ratio
    read(73,*) steps_before_toa_adj !number of timesteps before the whole temperature profile is adjusted to approach TOA equilibrium
    read(73,*) a_green !Green's analytic ozone expression: total ozone amount
    read(73,*) b_green !Green's analytic ozone expression: altitude of peak ozone
    read(73,*) c_green !Green's analytic ozone expression: ozone scale height
    read(73,*) H_green !Green's analytic ozone expression: atmosphere scale height
    read(73,*) cloudloctype !vertical location of cloud is in 1:altitude (km) 2:pressure (hPa) 3:temperature (K) 
    read(73,*) surf_emiss_on !emission from the surface 0:no emission
    read(73,*) lapse_type !type of lapse rate 0: critical lapse rate input 1: Held's parameterisation
    read(73,*) h2o_sb !H2O self-broadening 0:off 1:on
    read(73,*) h2o_for !H2O foreign broadening 0:off 1:on
    read(73,*) h2o_source ! 0=ERA-Interim mixh2o, 1=MW67 RH, 2=Cess, 3=Kasting, 4=Ramirez
    read(73,*) ur_mt !under-relaxation constant for meridional transport
    read(73,*) mtransp_type !source of meridional transport 1:simple diffusion 2:Vladilo
    read(73,*) steps_before_first_eqbcheck !number of steps before first adjustment to approach TOA equilibrium
    read(73,*) gas_addmolec_h2o
    read(73,*) gas_addmolec_co2
    read(73,*) gas_addmolec_o3
    read(73,*) max_rh
    read(73,*) snapshot
    read(73,*) piar
    read(73,*) pich4
    read(73,*) gas_amt_p_high_ch4
    read(73,*) gas_amt_p_low_ch4
    read(73,*) gas_amt_pert_ch4
    read(73,*) gas_amt_fac_ch4
    read(73,*) gas_addmolec_ch4
    read(73,*) forcing_expt
    read(73,*) o3_source
    read(73,*) t_intp_type

    close(73)


    ! open(81,file=('Input Distributions/fal lats'),form='formatted')

    ! !Read in lat profiles of 1D variables
    ! do col=1,ncols
    !     read(81,*) R_gcols(col)
    !     R_gcols(col) = R_gcols(col)
    ! enddo

    ! close(81)

    write(50,*) ncols,nlayersm,tot_albedo,solar_constant,tboundm,lapsecrit,days,mixco2_in,undrelax,cld,rmin,height_hc,&
    frac_hc,od_hc,height_mc,&
    frac_mc,od_mc,height_lc,frac_lc,od_lc,toa_precision,R_g,fixed_trop,OLR_layer,adj_speed,pb,fixed_sw_on,fixed_sw,fp,surf_rh,&
    adj_freq,convecttype,pb_new,swo3,swh2o,maxhtr,adjspeedfactor,pico2,pin2,pio2,htransp,ipe,detailprint,mtranspfac,&
    boxnetfluxfac,pertlay,pertcol,ncloudcols

    ! Calculate the average insolation in each latitude box
    do col = 1,ncols
        do day = 1,365
            do hour=1,24
                hourang = 15.0 * (float(hour)-12.0)
                declin = -23.5 * cosd(360.0/365.0 + (float(day) + 10.0))
                Xrad = boxlats(col) * 2.0*3.14/360.0
                Yrad = declin * 2.0*3.14/360.0
                Hrad = hourang * 2.0*3.14/360.0
                if (sin(Xrad)*sin(Yrad) + cos(Xrad)*cos(Yrad)*cos(Hrad) > 0) then
                    zen(col,day,hour) = acos(sin(Xrad)*sin(Yrad) + cos(Xrad)*cos(Yrad)*cos(Hrad))
                    insol(col,day,hour) = solar_constant * cos(zen(col,day,hour))
                else
                    insol(col,day,hour) = 0.0
                endif
            enddo
        enddo
    enddo   

    do col=1,ncols
        x_lats(col) = sin(boxlats(col)*3.14/180.0)
    enddo

    do i=0,ncols
        x_edge(i) = sin(latbounds(i)*3.14/180.0)
    enddo

    do col=1,ncols
        insolcols(col) = sum(insol(col,:,:)) / size(insol(col,:,:))
        zencols(col) = sum(zen(col,:,:)) / size(zen(col,:,:))
    enddo

    timesteps = int(days*ur_htr) ! number of time steps that the model runs for 
    tboundm = tboundmcols(1)

    ! undrelax is a single float, but now we have newur, which is an under-relaxation constant for each model layer
    do i=1,nlayersm
        newur(i) = ur_htr
    enddo

    totuflum = 0.
    totdflum = 0.
    fnetm = 0.
    htrm = 0.

    iatmm=0 !0 for layer values, 1 for level values
    ixsectm=0 !could be 1, but why?
    iscatm=0 !just absorption and emission
    !numangsm=4 !can be 0-4 for higher precision
    numangsm=0
    !ioutm=99 !for info in all spectral bands
    ioutm=0 !for broadband only
    !ioutm=-1 !for broadband, no printings
    !icldm=0 !for clear sky
    ! icldm=1  !for grey clouds
    iemissm=1 !surface emissivity. Keep this fixed for now.
    ireflecm=0 !for Lambert reflection
    semism=1. !all spectral bands the same as iemissm

    if (iatmm.ne.0) stop 'not ready for this'

    mixo2 = 0.21

    ! Gas inventories
    pin2 = pin2 * 1e5 !convert the input in bar to Pa
    pico2 = pico2 * 1e5 !convert the input in bar to Pa
    pio2 = pio2 * 1e5
    piar = piar * 1e5 !convert the input in bar to Pa
    pich4 = pich4 * 1e5 !convert the input in bar to Pa

    piair = pin2 + pio2 + piar

    massatmo_co2 = pico2 / gravity ! [kg]
    massatmo_n2 = pin2 / gravity ! [kg]
    massatmo_o2 = pio2 / gravity ! [kg]
    massatmo_ar = piar / gravity ! [kg]
    massatmo_ch4 = pich4 / gravity ! [kg]

    massatmo = massatmo_co2 + massatmo_n2 + massatmo_o2 + massatmo_ar + massatmo_ch4

    ! Gas mass mixing ratios 
    mass_mixco2 = massatmo_co2 / massatmo
    mass_mixn2 = massatmo_n2 / massatmo
    mass_mixo2 = massatmo_o2 / massatmo
    mass_mixar = massatmo_ar / massatmo
    mass_mixch4 = massatmo_ch4 / massatmo

    ! Number of molecules of each gas
    molec_co2 = massatmo_co2 / mmwco2 * avogadro
    molec_n2 = massatmo_n2 / mmwn2 * avogadro
    molec_o2 = massatmo_o2 / mmwo2 * avogadro
    molec_ar = massatmo_ar / mmwar * avogadro
    molec_ch4 = massatmo_ch4 / mmwch4 * avogadro

    totmolec = molec_co2 + molec_n2 + molec_o2 + molec_ar + molec_ch4

    ! Gas volume mixing ratios
    vol_mixco2 = molec_co2 / totmolec
    vol_mixn2 = molec_n2 / totmolec
    vol_mixo2 = molec_o2 / totmolec
    vol_mixar = molec_ar / totmolec
    vol_mixch4 = molec_ch4 / totmolec

    ! Mean molecular weight of the atmosphere
    mmwtot = mmwco2 * vol_mixco2 + mmwn2 * vol_mixn2 + mmwo2 * vol_mixo2 + mmwar*vol_mixar + mmwch4*vol_mixch4

    do i=1,nlayersm
        cpco2(i) = 0.000879 * tavelm(i) + 0.589073
        cvco2(i) = 0.000974 * tavelm(i) + 0.36620 ! Linear fits from Engineering Toolbox data (https://www.engineeringtoolbox.com/carbon-dioxide-d_974.html)
        rspecific_co2(i) = cpco2(i) - cvco2(i)
        cptot(i) = vol_mixco2 * cpco2(i) + vol_mixn2 * cpn2 + vol_mixo2 * cpo2
        cvtot(i) = vol_mixco2 * cvco2(i) + vol_mixn2 * cvn2 + vol_mixo2 * cvo2
        !        rsp_tot(i) = cptot(i) - cvtot(i)
        ! cptot(i) = 1.003 !default
        rsp_tot(i) = 0.287*1000. !default
    enddo

    nmolm=7 !number of molecular species, le 7 for now.

    if (psurf_override > 0.0) then 
        surfacep = psurf_override * 100.
    else
        surfacep = air_inv+co2_inv !Pascals
    end if

    if (ipe == 0) surfacep = air_inv !If you want to neglect the effect of CO2 on the surface pressure

    wklm = 0.; wbrodlm = 0.; tzm = 0.; pzm = 0.; altzm = 0.; tavelm = 0.

    pzm(0)=surfacep/100.
    tzm(0)=tboundm

    ! Initialise with an isothermal profile
    do i=1,nlayersm
        tavelm(i) = tboundm
        tzm(i) = tboundm
    enddo

    ! Surface pressure is just the sum of the gas inventories
    if (psurf_override > 0.0) then 
        surfacep = psurf_override * 100.
    else
        surfacep = air_inv+co2_inv !Pascals !NJE changed
    end if
    if (ipe == 0) surfacep = air_inv !If you want to neglect the effect of CO2 on the surface pressure
    pzm(0) = surfacep/100

    logpzm(0) = log(surfacep/100.0)
    sigma(0) = 1.0

    if (fp == 0) then
        do i =1,nlayersm
            select case(playtype)
            case(0)
                pzm(i)=pzm(i-1)-(surfacep-min_press*100.)/nlayersm/100. !in mb !Split the atmosphere into layers of equal pressure-thickness
            case(1)
                sigma(i)=sigma(i-1) - 0.9 / nlayersm
                pzm(i) = pzm(0) * sigma(i)**4.0 * (3.0 - 2.0*sigma(i))
            end select            
            if (pzm(i) < min_press) then
                pzm(i) = min_press
            endif            
            pavelm(i)=(pzm(i)+pzm(i-1))/2 !Set the average pressure for each layer
            altzm(i)=altzm(i-1)+(pzm(i-1)-pzm(i))*rsp_tot(i)*tavelm(i)/pavelm(i)/gravity !Rearrange the hydrostatic equation
            altlaym(i)=(altzm(i)+altzm(i-1))/2 !Set the altitude of the middle of the layer
        enddo
    endif

    do i=1,nlayersm
        mperlayr(i) = totmolec/nlayersm !Divide the molecules equally between layers
        mperlayr_air(i) = (molec_n2 + molec_o2)/(nlayersm)
        mperlayr_co2(i) = (molec_co2)/(nlayersm)
        mperlayr_co2(i) = (molec_co2)/(nlayersm)
    enddo

    !Set up mixing ratio of broadening molecules (N2 and O2 mostly)
    do i =1,nlayersm
        wbrodlm(i) = mperlayr_air(i) * 1.0e-4
        wklm(1,i) = mperlayr(i) * 1.0e-4 * mixh2o(i)
        if (mixco2_prescribed_on ==1) then
            wklm(2,i) = mperlayr(i) * 1.0e-4 * mixco2_prescribed
        else
            ! wklm(2,i) = mperlayr_co2(i) * 1.0e-4
            wklm(2,i) = mperlayr(i) * vol_mixco2
        endif
        wklm(3,i) = mperlayr(i) * 1.0e-4 * mixo3(i)
        wklm(6,i) = mperlayr(i) * 1.0e-4 * vol_mixch4
        wklm(7,i) = mperlayr(i) * 1.0e-4 * vol_mixo2
        ! wklm(3,i) = (2.69e19) * (1000. / pzm(0) ) * (u_lw(i-1) - u_lw(i)) 
    enddo

    ! Use Green's analytic ozone formula to calculate u_lw(i), which is the amount of ozone (in cm NTP) above layer i in the atmosphere. (The subscript lw means that we use this ozone amount for the longwave radiative transfer)
    do i=1,nlayersm
        u_lw(i) = ( a_green + a_green*exp(-b_green/c_green) ) / ( 1.0 + exp( (-H_green*log( pzm(i)/1000.0 ) -b_green ) /c_green ))
        if (pzm(i) < 0.0) u_lw(i) = 1.0e-4
    enddo

    do col=1,ncols
        totuflum = 0.
        totdflum = 0.
        fnetm = 0.
        htrm = 0.
        currentmaxhtr = 0.
        if (col ==1) then
            ! write(*,1109,advance='no')
        endif
        if (col < ncols) then
            ! write(*,1104,advance='no') boxlats(col)
        else
            ! write(*,1104) boxlats(col)
        endif
    enddo

    do col=1,ncols
        if (col < ncols) then
            ! write(*,1105,advance='no') '-------------'
        else
            ! write(*,1105) '-------------'
        endif
    enddo

    if(fp==1) then
        open(91,file='fixed_profile')
        do col=1,ncols
            do i=0,nlayersm
                read(91,1110) tzmcols(i,col)
            enddo
        enddo
    endif

    close(91)

    ! Main do loop - calculate the heating rates, adjust the temperature profile, and iterate until radiative (convective) equilibrium
    ! MAIN LOOP
    ! MAIN LOOP

    if(snapshot==1) timesteps=1

    do j=1,timesteps

        do col=1,ncols

            tboundm = tboundmcols(col)
            sol_inc = insolcols(col)
            lapsecrit = lapsecritcols(col)
            height_lc = height_lccols(col)
            frac_lc = frac_lccols(col)
            od_lc = od_lccols(col)
            surf_rh = surf_rhcols(col)
            R_g = R_gcols(col)
            mu_0 = zencols(col)

            tzm(0) = tboundm

            do i=1,nlayersm
                tavelm(i) = tboundm
                tzm(i) = tboundm
            enddo

            write(qfn,"(A83,I2)") 'Input Distributions/q vert col ', col-1
            write(o3fn,"(A84,I2)") 'Input Distributions/o3 vert col ', col-1
            write(t_fn,"(A83,I2)") 'Input Distributions/t vert col ', col-1

            write(ccfracsfn,"(A84,I2)") 'Input Distributions/ccfracs col '&
            &, col-1      
            write(cctausfn,"(A84,I2)") 'Input Distributions/cctaus col ', col-1   
            write(ccaltsfn,"(A84,I2)") 'Input Distributions/ccalts col ', col-1   

            open(82,file=(trim(qfn)),form='formatted')
            open(83,file=(trim(o3fn)),form='formatted')
            open(90,file=(trim(t_fn)),form='formatted')
            !            open(84,file=(trim(ccfn)),form='formatted')
            !            open(85,file=(trim(clwcfn)),form='formatted')
            !            open(86,file=(trim(ciwcfn)),form='formatted')

            open(87,file=trim(ccfracsfn),form='formatted')
            open(88,file=trim(cctausfn),form='formatted')
            open(89,file=trim(ccaltsfn),form='formatted')

            open(92,file=('extra_clds'),form='formatted')

            read(92,*) extra_cld_tau
            read(92,*) extra_cld_frac
            read(92,*) extra_cld_alt
            read(92,*) extra_cld_latcol
            read(92,*) extra_cld_cldcol

            close(92)


            do i=1,nlayersm
                select case(o3_source)
                case(1)
                    read(83,*) mixo3(i)
                    mixo3(i) = mixo3(i) * mmwtot / (48.0*1e-3)
                case(2)
                    mixo3(i) = (3.6478*(pzm(i)**0.83209))*exp(-pzm(i)/11.3515)*1e-6
                end select
                wklm(3,i) = mperlayr(i) * 1.0e-4 * mixo3(i)
                read(90,*) tzm(i)
                !                read(84,*) fracs(i)
                !                read(85,*) clwc(i)
                !                read(86,*) ciwc(i)
            enddo

            tzm(0) = tzm(1)
            ! tboundm = tzm(0) !don't need to do this since I change tboundm directly


            close(83)

            do cloudcol = 1,ncloudcols
                read(87,*) ccfracs(cloudcol)
                read(88,*) cctaus(cloudcol)
                read(89,*) ccalts(cloudcol)
            end do

            close(87)
            close(88)
            close(89)

            !Calculate liquid and ice water paths and cloud optical thicknesses from these following Claire Radley's thesis figs
            ! do i=1,nlayersm
            !     lwp(i) = clwc(i) * (pzm(i-1) - pzm(i)) / gravity * 1000.0 !liquid water path in gram metre^-2
            !     iwp(i) = ciwc(i) * (pzm(i-1) - pzm(i)) / gravity * 1000.0 !liquid water path in gram metre^-2
            !     tau_cld(i) = (lwp(i) / 15.0 + iwp(i) / 75.0) * 1500.0 !nje tcwp
            !     tau_cld_sw(i) = (lwp(i) / 15.0 + iwp(i) / 75.0) * 0.0
            ! enddo

            if (transpcalled == 1) then

                tboundm = tboundmcols(col) + tempchanges(col)
                if (couple_tgta == 1) then
                    tzm(0) = tboundm
                    tavelm(1) = tzm(0)
                end if

                do i=0,nlayersm
                    tzm(i) = tzmcols(i,col) + tempchanges(col)
                    if (i>1) tavelm(i) = tavelmcols(i,col) + tempchanges(col)
                enddo

            elseif (transpcalled == 0 .and. j > 1) then
                tboundm = tboundmcols(col)
                if (couple_tgta == 1) then
                    tzm(0) = tboundm
                    tavelm(1) = tzm(0)
                end if

                do i=1,nlayersm
                    tzm(i) = tzmcols(i,col)
                    tavelm(i) = tavelmcols(i,col) 
                enddo

            endif

            ! if (col==inversion_col) then
            !     tzm(0) = tboundm + inversion_strength
            ! endif

            sigma(0) = 1.0
            if (fp==1) then
                do i=0,nlayersm
                    tzm(i) = tzmcols(i,col)
                enddo
                do i =1,nlayersm
                    tavelm(i) = (tzm(i-1) + tzm(i)) / 2.0
                    select case(playtype)
                    case(0)
                        pzm(i)=pzm(i-1)-(surfacep-min_press*100.)/nlayersm/100. !in mb !Split the atmosphere into layers of equal pressure-thickness
                    case(1)
                        sigma(i)=sigma(i-1) - 0.9 / nlayersm
                        pzm(i) = pzm(0) * sigma(i)**4.0 * (3.0 - 2.0*sigma(i))
                    end select            
                    if (pzm(i) < min_press) then
                        pzm(i) = min_press
                    endif            
                    pavelm(i)=(pzm(i)+pzm(i-1))/2 !Set the average pressure for each layer
                    altzm(i)=altzm(i-1)+(pzm(i-1)-pzm(i))*rsp_tot(i)*tavelm(i)/pavelm(i)/gravity !Rearrange the hydrostatic equation
                    altlaym(i)=(altzm(i)+altzm(i-1))/2 !Set the altitude of the middle of the layer
                enddo
            endif 

            do i=1,nlayersm
                u_lw(i) = ( a_green + a_green*exp(-b_green/c_green) ) / ( 1.0 + exp((-H_green*log(pzm(i)/1000.0)-b_green)/c_green))
                if (pzm(i) < 0.0) u_lw(i) = 1.0e-4
            enddo

            coldpoint_trop_ind(col) = minloc(tzm(1:nlayersm),dim=1)

            do i=1,nlayersm
               esat_liq(i) = 6.1094*exp(17.625*(tzm(i)-273.15)/(tzm(i)-273.15+243.04))

                select case(h2o_source)
                case(0) !ERA-Interim
                    read(82,*) mixh2o(i)
                    mixh2o(i) = mixh2o(i) * mmwtot / (18.014*1e-3)
                    rel_hum_cols(i,col) = mixh2o(i) * pavelm(i) / (esat_liq(i) * (0.622 + mixh2o(i)) )
                case(1) !MW67
                    rel_hum(i) = surf_rh*(pzm(i)/1000.0 - 0.02)/(1.0-0.02)
                    mixh2o(i) = 0.622*rel_hum(i)*esat_liq(i)/(pavelm(i)-rel_hum(i)*esat_liq(i))
                case(2) !Cess
                    omega_rh = 1. - 0.03*(tzm(0)-288.)
                    rel_hum(i) = surf_rh*(pzm(i)/1000.0)**omega_rh
                    mixh2o(i) = 0.622*rel_hum(i)*es(i)/(pavelm(i)-rel_hum(i)*esat_liq(i))
                case(3) !Kasting and Ackerman
                    omega_rh = 1. - (esat_liq(1)/pzm(0) - 0.0166 ) / ( 0.1 - 0.0166 )
                    if (omega_rh > 1.) omega_rh=1
                    if (omega_rh < 0.) omega_rh=0
                    rel_hum(i) = surf_rh*(pzm(i)/1000.0 - 0.02)/(1.0-0.02)**omega_rh
                    mixh2o(i) = 0.622*rel_hum(i)*esat_liq(i)/(pavelm(i)-rel_hum(i)*esat_liq(i))
                case(4) !Ramirez 2013
                    omega_rh = 1. - (esat_liq(1)/pzm(0) - 0.0166 ) / ( 0.1 - 0.0166 )
                    if (omega_rh > 1.) omega_rh=1
                    if (omega_rh < 0.) omega_rh=0
                    density = pavelm(1) * 100. / (rsp_tot(1) * tavelm(1))
                    tcels = tzm(0) - 273.15
                    t1_vl = tcels + 0.5
                    t2_vl = tcels - 0.5
                    s_r13 = (6.1094 * exp(17.625*t1_vl/(t1_vl+243.04)) - 6.1094 * exp(17.625*t2_vl/(t2_vl+243.04)))
                    gamma_r13 = s_r13 / (s_r13+(cptot(1)/2.26e6))
                    bowen_r13 = (1.0 - gamma_r13*1.26) / (gamma_r13*1.26)
                    if (bowen_r13 < 0.1) bowen_r13 = 0.1
                    surf_rh = 1. - (1. - 0.8) * ( ( abs_surf_lh - fnetm(0) ) / ( 97. )) * ( 1.225 / density ) *&
                    & ( 0.0166 / (es(1)/pzm(0)) ) * ( ( 1. + 0.59 )/( 1. + bowen_r13 ) )
                    rel_hum(i) = surf_rh*(pzm(i)/1000.0 - 0.02)/(1.0-0.02)**omega_rh
                    mixh2o(i) = 0.622*rel_hum(i)*esat_liq(i)/(pavelm(i)-rel_hum(i)*esat_liq(i))
                case(5) !option to keep h2o constant in all columns
                    h2o_dummy_temp(i) = max(280.+altzm(i)/1000.*lapsecritcols(col), 200.)
                   esat_liq(i) = 6.1094*exp(17.625*(h2o_dummy_temp(i)-273.15)/&
                        &(h2o_dummy_temp(i)-273.15+243.04)) 
                    rel_hum(i) = 0.8*(pzm(i)/1000.0 - 0.02)/(1.0-0.02)
                    mixh2o(i) = 0.622*rel_hum(i)*esat_liq(i)/(pavelm(i)-rel_hum(i)*esat_liq(i))
                case(6)
                    alpha(i) = ( ( tzm(i) - 273.16 + 23 ) / ( 23 ) )**2
                    if (alpha(i) < 0) alpha(i) = 0.
                    if (alpha(i) > 1) alpha(i) = 1.
                    esat_ice(i) = 6.1121 * exp( 22.587 * ( tzm(i) - 273.16) / ( tzm(i) + 0.7 ) ) !note different e0 for ice and liq - should synthesise
                    esat_tot(i) = alpha(i)*esat_liq(i) + (1.-alpha(i)) * esat_ice(i)                    
                    rel_hum(i) = 0.8*(pzm(i)/1000.0 - 0.02)/(1.0-0.02)
                    mixh2o(i) = 0.622*rel_hum(i)*esat_tot(i)/(pavelm(i)-rel_hum(i)*esat_tot(i))
                    if (i>coldpoint_trop_ind(col)) then
                        mixh2o(i) = mixh2o(coldpoint_trop_ind(col))
                    endif
                end select
                if (mixh2o(i) < rmin) mixh2o(i) = rmin
                if (rel_hum(i) < 1e-3) rel_hum(i) = 1e-3
                if (rel_hum_cols(i,col) > max_rh) then
                    rel_hum_cols(i,col) = max_rh
                    mixh2o(i) = 0.622*max_rh*esat_liq(i)/(pavelm(i)-max_rh*esat_liq(i))
                end if
            enddo

            do i =1,nlayersm
                wbrodlm(i) = mperlayr_air(i) * 1.0e-4
                wklm(1,i) = mperlayr(i) * 1.0e-4 * mixh2o(i) 
                if (mixco2_prescribed_on ==1) then
                    wklm(2,i) = mperlayr(i) * 1.0e-4 * mixco2_prescribed
                else
                    wklm(2,i) = mperlayr_co2(i) * 1.0e-4
                endif
                wklm(3,i) = mperlayr(i) * 1.0e-4 * mixo3(i)
                if(pzm(i) < gas_amt_p_high_h2o .and. pzm(i) > gas_amt_p_low_h2o) then 
                    if (gas_amt_pert_h2o == 1) then
                        wklm(1,i) = wklm(1,i) * gas_amt_fac_h2o
                        wklm(1,i) = wklm(1,i) + gas_addmolec_h2o
                    end if
                end if
                if(pzm(i) < gas_amt_p_high_co2 .and. pzm(i) > gas_amt_p_low_co2) then 
                    if (gas_amt_pert_co2 == 1) then
                        wklm(2,i) = wklm(2,i) * gas_amt_fac_co2
                        wklm(2,i) = wklm(2,i) + gas_addmolec_co2
                    end if
                end if
                if(pzm(i) < gas_amt_p_high_o3 .and. pzm(i) > gas_amt_p_low_o3) then 
                    if (gas_amt_pert_o3 == 1) then
                        wklm(3,i) = wklm(3,i) * gas_amt_fac_o3
                        wklm(3,i) = wklm(3,i) + gas_addmolec_o3
                    end if
                end if
                if(pzm(i) < gas_amt_p_high_ch4 .and. pzm(i) > gas_amt_p_low_ch4) then 
                    if (gas_amt_pert_ch4 == 1) then
                        wklm(6,i) = wklm(6,i) * gas_amt_fac_ch4
                        wklm(6,i) = wklm(6,i) + gas_addmolec_ch4
                    end if
                end if
            enddo

            do i=1,nlayersm
                cpco2(i) = 0.000879 * tavelm(i) + 0.589073
                cvco2(i) = 0.000974 * tavelm(i) + 0.36620
                rspecific_co2(i) = cpco2(i) - cvco2(i)
                cptot(i) = vol_mixco2 * cpco2(i) + vol_mixn2 * cpn2 + vol_mixo2 * cpo2
                cvtot(i) = vol_mixco2 * cvco2(i) + vol_mixn2 * cvn2 + vol_mixo2 * cvo2
                rsp_tot(i) = cptot(i) - cvtot(i)
                ! cptot(i) = 1.003
                rsp_tot(i) = 0.287 * 1000.
            enddo

            ! If pressure broadening is off, feed a pressure profile with surface pressure = 1 bar to the subroutine that picks the correlated-k distribution
            if (pb_new == 0) then
                pzm(0) = 1000.0
                sigma(0) = 1.0
                do i =1,nlayersm
                    select case(playtype)
                    case(0)
                        pzm(i)=pzm(i-1)-(100000.-min_press)/nlayersm/100. !in mb !Split the atmosphere into layers of equal pressure-thickness
                    case(1)
                        sigma(i)=sigma(i-1) - 0.9 / nlayersm
                        pzm(i) = 1000. * sigma(i)**4.0 * (3.0 - 2.0*sigma(i))
                    end select            
                    if (pzm(i) < min_press) then
                        pzm(i) = min_press
                    endif            
                    pavelm(i)=(pzm(i)+pzm(i-1))/2 !Set the average pressure for each layer
                enddo
            endif 

            htrmwghtd = 0.
            totuflumwghtd = 0.
            totdflumwghtd = 0.
            htrlhwghtd = 0.
            htro3_lhwghtd = 0.

            abspnwghtd = 0.
            A_oz_lwghtd = 0.
            abs_surf_lhwghtd = 0.
            abs_h2owghtd = 0.
            abs_o3wghtd = 0.

            tot_sol_abs_lhwghtd = 0.

            do cloudcol = 1,ncloudcols

                cloudcolfrac = ccfracs(cloudcol)
                cloudcoltau = cctaus(cloudcol)
                cloudcolalt = ccalts(cloudcol)

                call createcloudfile(cloudcolalt,cloudcoltau)

                ! Turn surface emission off
                if(surf_emiss_on == 0) then
                    tboundm = tboundm/3.0
                end if

                ! Call the subroutine rrtm, which calculates the upward and downward fluxes and the heating rates
                call rrtm

                if(surf_emiss_on == 0) then
                    tboundm = tboundm*3.0
                end if

                call calcswhtr

                htrmwghtd = htrmwghtd + htrm * cloudcolfrac
                totuflumwghtd = totuflumwghtd + totuflum * cloudcolfrac
                totdflumwghtd = totdflumwghtd + totdflum * cloudcolfrac
                htrlhwghtd  = htrlhwghtd + htrlh * cloudcolfrac
                htro3_lhwghtd  = htro3_lhwghtd + htro3_lh * cloudcolfrac

                abspnwghtd = abspnwghtd + abspn * cloudcolfrac
                A_oz_lwghtd  = A_oz_lwghtd + A_oz_l * cloudcolfrac

                abs_h2owghtd = abs_h2owghtd + abs_h2o * cloudcolfrac
                abs_o3wghtd = abs_o3wghtd + abs_o3 * cloudcolfrac
                abs_surf_lhwghtd = abs_surf_lhwghtd + abs_surf_lh * cloudcolfrac

                tot_sol_abs_lhwghtd = tot_sol_abs_lhwghtd + tot_sol_abs_lh * cloudcolfrac

            enddo !cloudcol

            htrm = htrmwghtd
            totuflum = totuflumwghtd
            totdflum = totdflumwghtd
            htrlh = htrlhwghtd
            htro3_lh = htro3_lhwghtd
            abspn = abspnwghtd
            A_oz_l = A_oz_lwghtd
            abs_surf_lh = abs_surf_lhwghtd
            tot_sol_abs_lh = tot_sol_abs_lhwghtd
            abs_h2o = abs_h2owghtd
            abs_o3 = abs_o3wghtd
            abs_surf_lh = abs_surf_lhwghtd

            ! call add_seb_to_tboundm
            if (couple_tgta == 1) then
                tzm(0) = tboundm
                tavelm(1) = tzm(0)
            end if

            !Indexing to match indexing scheme of Lacis and Hansen (which is different, not sure why)
            do i=1,nlayersm
                if(swh2o == 1) htrm(i-1) = htrm(i-1) + htrlh(i)
                if(swo3 == 1 .and. i > 1)  htrm(i) = htrm(i) + htro3_lh(i-1)
            enddo


            ! Undo the if block from before the call to rrtm (i.e., set the surface pressure back to whatever it was before it was set to 1 bar)
            if (pb_new == 0) then
                pzm(0)=surfacep/100.
                sigma(0) = 1.0
                do i=1,nlayersm
                    select case(playtype)
                    case(0)
                        pzm(i)=pzm(i-1)-(surfacep-min_press*100.)/nlayersm/100. !in mb !Split the atmosphere into layers of equal pressure-thickness
                    case(1)
                        sigma(i)=sigma(i-1) - 0.9 / nlayersm
                        pzm(i) = pzm(0) * sigma(i)**4.0 * (3.0 - 2.0*sigma(i))
                    end select            
                    if (pzm(i) < min_press) then
                        pzm(i) = min_press
                    endif            
                    pavelm(i)=(pzm(i)+pzm(i-1))/2 !Set the average pressure for each layer
                enddo
            endif

            if (fixed_sw_on == 1) then
                abs_sw(col) = fixed_sw ! Feed in a fixed value of absorbed SW (for the purposes of calculating Fnet at the TOA)
            else
                abs_sw(col) = tot_sol_abs_lh ! Use the total absorbed SW 
                abs_h2o_cols(col) = abs_h2o
                abs_o3_cols(col) = abs_o3
                abs_surf_cols(col) = abs_surf_lh
            endif


            ! Auto-T

            !            htrmaxloc = maxloc(abs(htrm(conv_trop_ind(col)+4:nlayersm-4)),dim=1)
            !            currentmaxhtr1 = maxval(abs(htrm(conv_trop_ind(col)+4:htrmaxloc)))
            !            currentmaxhtr2 = maxval(abs(htrm(htrmaxloc:nlayersm-4)))
            !            
            !            currentmaxhtr = max(currentmaxhtr1, currentmaxhtr2 )

            ! currentmaxhtr = maxval( abs( htrm( conv_trop_ind(col)+1:nlayersm-5 ) ) ) ! Current maximum heating rate for the whole atmosphere
            currentmaxhtr = maxval( abs( htrm( conv_trop_ind(col)+1:nlayersm-1 ) ) )
            !            if (conv_trop_ind(col)+4 > nlayersm-4) then 
            !              currentmaxhtr = maxval(abs(htrm(conv_trop_ind(col)+1:-1)))
            !            endif


            stepssinceadj = stepssinceadj + 1 

            ! Progressively decrease the under-relaxation constant as currentmaxhtr increases (which means equilibrium is approaching)
            do i=1,nlayersm
                newur(i) = ur_htr / 1.0 * (pzm(0) - pzm(1)) / (pzm(i-1) - pzm(i))
                if (currentmaxhtr < 10.0 .and. stepssinceadj > 5 .and. adj1<1) then
                    newur(i) = ur_htr / 2.0 * (pzm(0) - pzm(1)) / (pzm(i-1) - pzm(i))
                    if (i==nlayersm) adj1 = 1
                endif
                if (currentmaxhtr < 5.0 .and. stepssinceadj > 5 .and. adj1<1) then
                    newur(i) = ur_htr / 4.0 * (pzm(0) - pzm(1)) / (pzm(i-1) - pzm(i))
                    if (i==nlayersm) adj1 = 1
                endif
                if (currentmaxhtr < 2.0 .and. stepssinceadj > 5 .and. adj1<1) then
                    newur(i) = ur_htr / 8.0 * (pzm(0) - pzm(1)) / (pzm(i-1) - pzm(i))
                    if (i==nlayersm) adj1 = 1
                endif
                if (currentmaxhtr < 0.5 .and. stepssinceadj > 5 .and. adj2<1) then
                    newur(i) = ur_htr / 16.0 * (pzm(0) - pzm(1)) / (pzm(i-1) - pzm(i))
                    if (i==nlayersm) adj2 = 1
                endif
                if (currentmaxhtr < 0.2 .and. stepssinceadj > 5 .and. adj3<1) then
                    newur(i) = ur_htr / 32.0 * (pzm(0) - pzm(1)) / (pzm(i-1) - pzm(i))
                    if (i==nlayersm) adj3 = 1
                endif
            enddo


            ! htrm(i-1) is used because rtr.f assigns htr(L) to layer L, when it should actually heat layer L+1 (which rtr.f calls LEV)
            ! applyhtr
            if(snapshot /= 1) then
                do i=1,nlayersm
                    if(htrm(i-1)/(newur(i)) > 5.) then
                        tavelm(i) = tavelm(i)+5.
                    else if(htrm(i-1)/(newur(i)) < -5.) then
                        tavelm(i) = tavelm(i)-5.
                    else
                        tavelm(i) = tavelm(i) + htrm(i-1)/(newur(i))
                        ! tavelm(i) = tavelm(i) + (totuflum(i)-totuflum(i-1)-(totdflum(i)-totdflum(i-1)))/(pzm(i)-pzm(i-1))*8.49
                    end if
                    if (tavelm(i) < t_min) then 
                        tavelm(i) = t_min
                    end if
                end do
            end if

            ! do i=1,nlayersm
            !   tzm(i) = tzm(i) + htrm(i-1)/(newur(i))
            ! ENDDO

            ! tboundm = tboundmcols(col)
            ! if (col==inversion_col) then
            !     tzm(0) = tboundm + inversion_strength
            ! endif


            if (couple_tgta == 1) then
                tzm(0) = tboundm
                tavelm(1) = tzm(0)
            end if


            if(snapshot/=1) then
                select case(t_intp_type)                
                case(0)
                    do i=1,nlayersm-1
                        tzm(i)=(tavelm(i)+tavelm(i+1))/2
                        if (tzm(i) < t_min) tzm(i) = t_min
                    enddo
                    !Set temperature of very top level
                    tzm(nlayersm) = 2.0*tavelm(nlayersm)-tzm(nlayersm-1)
                case(1)
                    ! quadratic interpolation of tzm
                    do i=2,nlayersm-1
                        tzm(i)=&
((pzm(i)-pavelm(i))*(pzm(i)-pavelm(i+1)))/((pavelm(i-1)-pavelm(i))*(pavelm(i-1)-pavelm(i+1)))*tavelm(i-1)+&
((pzm(i)-pavelm(i-1))*(pzm(i)-pavelm(i+1)))/((pavelm(i)-pavelm(i-1))*(pavelm(i)-pavelm(i+1)))*tavelm(i)+&
((pzm(i)-pavelm(i-1))*(pzm(i)-pavelm(i)))/((pavelm(i+1)-pavelm(i-1))*(pavelm(i+1)-pavelm(i)))*tavelm(i+1)
                    end do
                    tzm(nlayersm)=&
((pzm(nlayersm)-pavelm(i-1))*(pzm(nlayersm)-pavelm(i)))/((pavelm(i-2)-pavelm(i-1))*(pavelm(i-2)-pavelm(i)))*tavelm(i-2)+&
((pzm(nlayersm)-pavelm(i-2))*(pzm(nlayersm)-pavelm(i)))/((pavelm(i-1)-pavelm(i-2))*(pavelm(i-1)-pavelm(i)))*tavelm(i-1)+&
((pzm(nlayersm)-pavelm(i-2))*(pzm(nlayersm)-pavelm(i-1)))/((pavelm(i)-pavelm(i-2))*(pavelm(i)-pavelm(i-1)))*tavelm(i)
                end select


                conv(:,col) = 0

                ! NJE Level Convection:
                conv_trop_ind(col) = 0
                do i=1,nlayersm
                    malr(i) = gravity * ( 1.0 + ( hv * mixh2o(i) / ( rsp_tot(i) * tavelm(i) ) ) )&
                    /( cptot(i) + ( (hv**2) * mixh2o(i) / (461.5 * (tavelm(i)**2.0) ) ) ) * (-1.0) * htransp ! Calculate the moist adiabatic lapse rate  
                end do

                !            call levconvect(col,tzm,altzm,altlaym,lapsecrit,fixed_trop,conv,int(convecttype),malr)
                call levconvect !nje t0

                conv_trop_ind(col) = minloc(conv(:,col),dim=1)
                coldpoint_trop_ind(col) = minloc(tzmcols(:,col),dim=1)

                ! do i=1,nlayersm
                !     tavelm(i) = (tzm(i-1) + tzm(i)) / 2.0
                ! end do

                do i=1,conv_trop_ind(col)
                    tavelm(i) = (tzm(i-1) + tzm(i)) / 2.0
                end do

            end if

            sigma(0) = 1.0

            do i =1,nlayersm
                !                logpzm(i) = logpzm(i-1) - (logpzm(0)-log(min_press/100.0))/nlayersm
                !                pzm(i) = exp(logpzm(i))
                select case(playtype)
                case(0)
                    pzm(i)=pzm(i-1)-(surfacep-min_press*100.)/nlayersm/100. !in mb !Split the atmosphere into layers of equal pressure-thickness
                case(1)
                    sigma(i)=sigma(i-1) - 0.9 / nlayersm
                    pzm(i) = pzm(0) * sigma(i)**4.0 * (3.0 - 2.0*sigma(i))
                end select            
                if (pzm(i) < min_press) then
                    pzm(i) = min_press
                endif
                pavelm(i)=(pzm(i)+pzm(i-1))/2.0
                altzm(i)=altzm(i-1)+(pzm(i-1)-pzm(i))*rsp_tot(i)*tavelm(i)/pavelm(i)/gravity
            enddo

            do i=1,nlayersm
                altlaym(i)=(altzm(i)+altzm(i-1))/2.0
            enddo

            ! Calculate quantities relevant to isentropic coordinates

            ! do i=0,nlayersm
            !     kappa(i) = rsp_tot(i)/cptot(i)
            !     theta(i) = tzm(i) * ((pzm(0) / pzm(i)) ** kappa(i)) !potential temperature
            !     exner(i) = cptot(i) * ((pzm(i) / pzm(0)) ** kappa(i))
            !     geopotential(i) = gravity * altzm(i)
            !     montgomery(i) = cptot(i) * tzm(i) + geopotential(i)
            ! enddo

            tzmcols(:,col) = tzm
            altzmcols(:,col) = altzm
            tavelmcols(:,col) = tavelm
            totuflumcols(:,col) = totuflum
            totdflumcols(:,col) = totdflum
            htrmcols(:,col) = htrm
            htrh2ocols(:,col) = htrlh
            htro3cols(:,col) = htro3_lh
            wklm1cols(:,col) = wklm(1,:)
            wklm2cols(:,col) = wklm(2,:)
            wklm3cols(:,col) = wklm(3,:)
            wbrodlmcols(:,col) = wbrodlm
            abspncols(:,col) = abspn
            A_oz_lcols(:,col) = A_oz_l
            abs_surf_lhcols(col) = abs_surf_lh

            altlaymcols(:,col) = altlaym
            pzmcols(:,col) = pzm
            pavelmcols(:,col) = pavelm
            tavelmcols(:,col) = tavelm
            tboundmcols(col) = tboundm
            tau_cldcols(:,col) = tau_cld
            fracscols(:,col) = cld_fracs
            R_gcols(col) = R_g
            fnetmcols(:,col) = fnetm


            olrcols(col) = totuflum(OLR_layer)
            currentmaxhtrcols(col) = currentmaxhtr
            tboundmcols(col) = tboundm

            if (detailprint==1) then

                if (col ==1) then
                    write(*,1107,advance='no') j
                endif


                if (col < ncols) then
                    write(*,1102, advance='no')  conv_trop_ind(col), currentmaxhtrcols(col)
                else
                    write(*,1102) conv_trop_ind(col), currentmaxhtrcols(col)
                endif



            endif

            ! Calculate surface energy budget (SEB) and perturb tboundm if chosen
            lhf = 10.
            shf = 30.
            bowen = 0.2125
            seb = 0.
            c_drag = 0.001
            meanwind = 5.0
            density = pavelm(1) * 100. / (rsp_tot(1) * tavelm(1))
            shf = cptot(1) * density * c_drag * meanwind * (tboundm - tavelm(1))
            lhf = (Lv * c_drag * meanwind) / (461.52 * tzm(0)) * &
            & 6.1094*100. * &
            &( exp( 17.625 * ( tboundm - 273.15 ) / ( tboundm - 273.15 + 243.04 ) ) -&
                & rel_hum(1) * exp( 17.625 * ( tzm(0) - 273.15 ) / ( tzm(0) - 273.15 + 243.04 ) ) ) 
            seb = totdflum(0) + abs_surf_lh - totuflum(0) - lhf - shf
            sebcols(col) = seb
            select case(sfc_heating)
            case(0) ! SEB doesn't affect surface temperature
                tboundm = tboundm 
            case(1) ! SEB warms/cools surface
                tboundm = tboundm + seb / ur_seb   
            end select
            tboundmcols(col) = tboundm

        end do !columns do loop

        tglobsum = 0.
        tglobmean = 0.
        cossums = 0.

        do col=1,ncols
            tglobsum = tglobsum + tboundmcols(col) * cosd(boxlats(col))
            cossums = cossums + cosd(boxlats(col))
        enddo

        tglobmean = tglobsum / cossums


        if (ncols > 1) then
            do i=1,ncols-1
                tair_lowest_edges(i) = (tavelmcols(1,i+1) + tavelmcols(1,i)) / 2.0
            enddo

            tair_lowest_edges(0) = 2.0 * tavelmcols(1,1) -  tair_lowest_edges(1)
            tair_lowest_edges(ncols) = 2.0 * tavelmcols(1,ncols) - tair_lowest_edges(ncols-1)
        end if



        do col=1,ncols
            ddry(col) = ks * cptot(1) * eta**(3.0/5) * cos(phim)**(-4.0/5) * planet_radius**(-6.0/5) * pzm(0)*100.0 / gravity *&
            & planet_rotation**(-4.0/5) * ( (twarm-tcold) / twarm * tot_sol_abs_lh/(pzm(0)*100.0/gravity))**(3.0/5)
            tcels = twarm - 273.15
            t1_vl = tcels + 0.1
            t2_vl = tcels - 0.1
            delta_pv_star = (6.1094 * exp(17.625*t1_vl/(t1_vl+243.04)) - 6.1094 * exp(17.625*t2_vl/(t2_vl+243.04)))
            lambda(col) = (kl * Lv * mmwh2o * rel_hum(1)) / (ks * cptot(1) * mmwtot * pzm(0)*100.0) * delta_pv_star/0.2
            d_vl(col) = ddry(col) * (1.0 + lambda(col))


            ! if(col==ncols) then
            !     ! delta_x_lats = x_lats(col-1) - x_lats(col) 
            !     delta_x_lats = boxlats(col-1) - boxlats(col) 
            !     delta_temp = (tboundmcols(col-1) - tboundmcols(col))
            ! else
            !     ! delta_x_lats = x_lats(col) - x_lats(col+1)
            !     delta_x_lats = boxlats(col) - boxlats(col+1) 
            !     delta_temp = (tboundmcols(col) - tboundmcols(col+1))
            ! endif

            ! meridtransp(col) = - 1.0 * (d_vl(col) * delta_temp / delta_x_lats) * cos(boxlats(col)*3.14/180.0)

            !  tot_transp(col) = -1.0 * (d_vl(col)*(1.0 - (x_lats(col)**2.0) * delta_temp/delta_x_lats))
            ! if(col == ncols) then
            !     meridtransp(col) = (tot_transp(col-1) - tot_transp(col)) / delta_x_lats
            ! else
            !     meridtransp(col) = (tot_transp(col) - tot_transp(col+1)) / delta_x_lats
            ! endif
            ! print*, 'delta_x_lats, delta_temp, meridtransp(col),tot_transp(col)'
            ! print*, delta_x_lats, delta_temp, meridtransp(col),tot_transp(col)
        enddo


        ! r_earth = 6371e3 !earth radius in m

        ! do i=1,ncols-1
        !     delta_T_edge(i) = tair_lowest_edges(i) - tair_lowest_edges(i+1)
        !     delta_x_edge(i) = x_edge(i) - x_edge(i+1)
        !     delta_y_edge(i) = 2. * r_earth * tand( (latbounds(i) - latbounds(i-1))/2. )
        !     h_scale = 7.5
        !     f_cor = -1.0 * abs(2. * 7.29e-5 * sind( latbounds(i) )) !check where this abs() should go NJE task
        !     beta = 2. * 7.29e-5 * cosd( latbounds(i) ) / r_earth
        !     gamma_d = 9.8
        !     d_mid(i) = h_scale * log( 1. - f_cor * abs(delta_T_edge(i)) / delta_y_edge(i) / ( h_scale * beta * &
        !         &( gamma_d + lapsecritcols(i) ) ) )
        !     d_trop(i) = wklm1cols(1,i) / wbrodlmcols(1,i) * Lv / ( cptot(1) * ( gamma_d + lapsecritcols(i) ) )
        !     meridtransp_edge(i) = delta_T_edge(i) / delta_x_edge(i) * (1.0 - (x_edge(i))**2.0) * d_vl(i)
        !     print*, i, latbounds(i), d_mid(i), d_trop(i), altzm(conv_trop_ind(i))/1000.
        ! enddo

        meridtransp_edge(0) = 0.0
        meridtransp_edge(ncols) = 0.0

        do col=1,ncols
            delta_meridtransp_edge(col) = (meridtransp_edge(col) - meridtransp_edge(col-1)) !* 5.0
            ! if (j==70) then
            !     print*, ddry(col), lambda(col), meridtransp_edge(col), delta_meridtransp_edge(col)
            ! end if
        enddo

        if (j==timesteps) then

            call writeoutputfile

            do col=1,ncols

                write(ccfracsfn,"(A84,I2)") 'Input Distributions/ccfracs col '&
                &, col-1      
                write(cctausfn,"(A84,I2)") 'Input Distributions/cctaus col ', &
                &col-1   
                write(ccaltsfn,"(A84,I2)") 'Input Distributions/ccalts col ', &
                &col-1   

                open(87,file=trim(ccfracsfn),form='formatted')
                open(88,file=trim(cctausfn),form='formatted')
                open(89,file=trim(ccaltsfn),form='formatted')

                open(92,file=('extra_clds'),form='formatted')

                read(92,*) extra_cld_tau
                read(92,*) extra_cld_frac
                read(92,*) extra_cld_alt
                read(92,*) extra_cld_latcol
                read(92,*) extra_cld_cldcol

                close(92)

                do cloudcol = 1,ncloudcols
                    read(87,*) ccfracs(cloudcol)
                    read(88,*) cctaus(cloudcol)
                    read(89,*) ccalts(cloudcol)
                end do

                close(87)
                close(88)
                close(89)

                do cloudcol = 1,ncloudcols
                    cloudcolfrac = ccfracs(cloudcol)
                    cloudcoltau = cctaus(cloudcol)
                    cloudcolalt = ccalts(cloudcol)
                    call writecloudstofile
                enddo

            enddo

            print*,
            print*, 'Reached maximum timestep'
            print*, 
            print*, ('----------------------------------------------')
            print*, 'Global Mean Temperature: ', tglobmean
            print*,
            transpcalled = 1
            stepssinceboxadj = 0

            do col=1,ncols
                boxnetradflux(col) = abs_sw(col)-olrcols(col)
                ! meridtransp(col) = (tglobmean - tboundmcols(col)) * mtranspfac
                ! meridtransp(col) = delta_meridtransp_edge(col)
                ddry(col) = ks * cptot(1) * eta**(3.0/5) * cos(phim)**(-4.0/5) * planet_radius**(-6.0/5) * pzm(0)*100.0 / gravity *&
                & planet_rotation**(-4.0/5) * ( (twarm-tcold) / twarm * tot_sol_abs_lh/(pzm(0)*100.0/gravity))**(3.0/5)
                tcels = twarm - 273.15
                t1_vl = tcels + 0.1
                t2_vl = tcels - 0.1
                delta_pv_star = (6.1094 * exp(17.625*t1_vl/(t1_vl+243.04)) - 6.1094 * exp(17.625*t2_vl/(t2_vl+243.04)))
                lambda(col) = (kl * Lv * mmwh2o * rel_hum(1)) / (ks * cptot(1) * mmwtot * pzm(0)*100.0) * delta_pv_star/0.1
                d_vl(col) = ddry(col) * (1.0 + lambda(col))

                print*, boxnetradflux(col), boxnetradflux_prev(col)
                ! if (boxnetradflux(col) / boxnetradflux_prev(col) < 0.0) then 
                !     ur_toafnet = ur_toafnet * 2.0
                !     print*, 'ur_toafnet increased to: ', ur_toafnet
                ! end if
                boxnetradflux_prev(col) = boxnetradflux(col)

                ! delta_x_lat = x_lat(col) - x_lat(col-1)

                ! meridtransp(col) = -1/(x_lat(col))

                if (mtranspon == 1) then
                    boxnettotflux(col) = boxnetradflux(col) + meridtransp(col)
                else
                    boxnettotflux(col) = boxnetradflux(col)
                end if
                tempchanges(col) = boxnettotflux(col) / ur_toafnet(col)
            enddo


            xglobsum = 0.
            xglobmean = 0.
            cossums = 0.
            do col=1,ncols
                xglobsum = xglobsum + boxnetradflux(col) * cosd(boxlats(col))
                cossums = cossums + cosd(boxlats(col))
            enddo
            xglobmean = tglobsum / cossums
            netradflux_globmean = xglobmean

            xglobsum = 0.
            xglobmean = 0.
            cossums = 0.
            do col=1,ncols
                xglobsum = xglobsum + meridtransp(col) * cosd(boxlats(col))
                cossums = cossums + cosd(boxlats(col))
            enddo
            xglobmean = tglobsum / cossums
            meridtransp_globmean = xglobmean

            xglobsum = 0.
            xglobmean = 0.
            cossums = 0.
            do col=1,ncols
                xglobsum = xglobsum + boxnettotflux(col) * cosd(boxlats(col))
                cossums = cossums + cosd(boxlats(col))
            enddo
            xglobmean = tglobsum / cossums
            nettotflux_globmean = xglobmean


            write(*,1106,advance='no') ''

            do col=1,ncols+1
                if (col < ncols+1) then
                    write(*,1104,advance='no') boxlats(col)
                else
                    write(*,*) 'AVG'
                endif
            enddo



            do col=1,ncols
                if (col < ncols) then
                    write(*,1105,advance='no') '-------------'
                else
                    write(*,1105) '-------------'
                endif
            enddo

            write(*,1106,advance='no') 'Box Net Radiative Flux | '
            do col=1,ncols+1
                if (col < ncols+1) then
                    write(*,1103,advance='no') boxnetradflux(col)
                else
                    write(*,1103) sum(boxnetradflux)/ncols, netradflux_globmean
                endif
            enddo

            write(*,1106,advance='no') 'Meridional Transport | '
            do col=1,ncols+1
                if (col < ncols+1) then
                    write(*,1103,advance='no') meridtransp(col)
                else
                    write(*,1103) sum(meridtransp)/ncols, meridtransp_globmean
                endif
            enddo

            write(*,1106,advance='no') 'Box Total Net Flux | '
            do col=1,ncols+1
                if (col < ncols+1) then
                    write(*,1103,advance='no') boxnettotflux(col)
                else
                    write(*,1103) sum(boxnettotflux)/ncols, nettotflux_globmean
                endif
            enddo

            write(*,1106,advance='no') 'Box Surface Temperature | '
            do col=1,ncols+1
                if (col < ncols+1) then
                    write(*,1103,advance='no') tboundmcols(col)
                else
                    write(*,1103) sum(tboundmcols)/ncols
                endif
            enddo

            write(*,1106,advance='no') 'Box OLR | '
            do col=1,ncols+1
                if (col < ncols+1) then
                    write(*,1103,advance='no') olrcols(col)
                else
                    write(*,1103) sum(olrcols)/ncols
                endif
            enddo

            write(*,1106,advance='no') 'Box Abs SW | '
            do col=1,ncols+1
                if (col < ncols+1) then
                    write(*,1103,advance='no') abs_sw(col)
                else
                    write(*,1103) sum(abs_sw)/ncols
                endif
            enddo

            write(*,1106,advance='no') 'Box Abs SW H2O | '
            do col=1,ncols+1
                if (col < ncols+1) then
                    write(*,1103,advance='no') abs_h2o_cols(col)
                else
                    write(*,1103) sum(abs_h2o_cols)/ncols
                endif
            enddo

            write(*,1106,advance='no') 'Box Abs SW O3 | '
            do col=1,ncols+1
                if (col < ncols+1) then
                    write(*,1103,advance='no') abs_o3_cols(col)
                else
                    write(*,1103) sum(abs_o3_cols)/ncols
                endif
            enddo

            write(*,1106,advance='no') 'Box Abs SFC | '
            do col=1,ncols+1
                if (col < ncols+1) then
                    write(*,1103,advance='no') abs_surf_cols(col)
                else
                    write(*,1103) sum(abs_surf_cols)/ncols
                endif
            enddo

            write(*,1106,advance='no') 'Box SEB | '
            do col=1,ncols+1
                if (col < ncols+1) then
                    write(*,1103,advance='no') sebcols(col)
                else
                    write(*,1103) sum(sebcols)/ncols
                endif
            enddo

            print*, ('----------------------------------------------')
            print*,

        endif

        transpcalled = 0

        if (j==400) then
            print*, 'j=400'
        end if


        ! Equilibrium check (eqbcheck)
        if (j > steps_before_first_eqbcheck .or. maxval(currentmaxhtrcols) < maxhtr ) then !NJE
            if ((maxval(currentmaxhtrcols) < maxhtr .and. stepssinceboxadj > 5) .or. stepssinceboxadj > steps_before_toa_adj)then
                ! if (stepssinceboxadj > steps_before_toa_adj) then
                print*, 
                print*, ('----------------------------------------------')
                print*, 'Global Mean Temperature: ', tglobmean
                print*,
                transpcalled = 1
                stepssinceboxadj = 0

                meridtransp_edge(0) = 0.0
                meridtransp_edge(ncols) = 0.0

                do col=1,ncols
                    ddry(col) = ks * cptot(1) * eta**(3.0/5) * cos(phim)**(-4.0/5) * planet_radius**(-6.0/5) * pzm(0)*100.0 / &
                    &gravity *planet_rotation**(-4.0/5) * ( (twarm-tcold) / twarm * tot_sol_abs_lh/(pzm(0)*100.0/gravity))&
                    &**(3.0/5.0)
                    tcels = twarm - 273.15
                    t1_vl = tcels + 0.1
                    t2_vl = tcels - 0.1
                    delta_pv_star = (6.1094 * exp(17.625*t1_vl/(t1_vl+243.04)) - 6.1094 * exp(17.625*t2_vl/(t2_vl+243.04)))
                    lambda(col) = (kl * Lv * mmwh2o * rel_hum(1)) / (ks * cptot(1) * mmwtot * pzm(0)*100.0) * delta_pv_star/0.2
                    d_vl(col) = ddry(col) * (1.0 + lambda(col))
                end do

                do edge = 1,ncols-1
                    delta_T_edge(edge) = tavelmcols(1,edge+1) - tavelmcols(1,edge)
                    delta_x_edge(edge) = x_lats(edge+1) - x_lats(edge)
                    meridtransp_edge(edge) = delta_T_edge(edge) / delta_x_edge(edge) * (1.0 - (x_edge(edge))**2.0) * d_vl(edge)
                end do

                do col=1,ncols
                    delta_meridtransp_edge(col) = meridtransp_edge(col) - meridtransp_edge(col-1)
                    if (mtransp_type ==  2) then
                        meridtransp(col) = delta_meridtransp_edge(col) * 1.0
                    end if
                end do



                ! delta_T_edge(0) = tavelmcols(1,1) - tair_lowest_edges(0)
                ! delta_T_edge(ncols) = tavelmcols(1,ncols) - tair_lowest_edges(ncols)
                ! delta_x_edge(0) = x_lats(1) - x_edge(0)
                ! delta_x_edge(ncols) = x_lats(ncols) - x_edge(ncols)
                do col=1,ncols
                    boxnetradflux(col) = abs_sw(col)-olrcols(col)
                    if (mtransp_type == 1) then
                        meridtransp(col) = (tglobmean - tboundmcols(col)) * mtranspfac
                    end if
                    ddry(col) = ks * cptot(1) * eta**(3.0/5) * cos(phim)**(-4.0/5) * planet_radius**(-6.0/5) * pzm(0)*100.0 / &
                    &gravity *planet_rotation**(-4.0/5) * ( (twarm-tcold) / twarm * tot_sol_abs_lh/(pzm(0)*100.0/gravity))&
                    &**(3.0/5.0)
                    tcels = twarm - 273.15
                    t1_vl = tcels + 0.1
                    t2_vl = tcels - 0.1
                    delta_pv_star = (6.1094 * exp(17.625*t1_vl/(t1_vl+243.04)) - 6.1094 * exp(17.625*t2_vl/(t2_vl+243.04)))
                    lambda(col) = (kl * Lv * mmwh2o * rel_hum(1)) / (ks * cptot(1) * mmwtot * pzm(0)*100.0) * delta_pv_star/0.2
                    d_vl(col) = ddry(col) * (1.0 + lambda(col))
                    ! if (col < ncols) then
                    !     delta_T_edge(col) = tair_lowest_edges(col+1) - tair_lowest_edges(col)
                    !     delta_x_edge(col) = x_edge(col+1) - x_edge(col)
                    ! end if
                    delta_y_edge(col) = r_earth * (latbounds(col) - latbounds(col-1) ) * 3.14 / 180.
                    h_scale = 7.5
                    f_cor = 2. * 7.29e-5 * sind( boxlats(col) ) !check where this abs() should go NJE task
                    beta = 2. * 7.29e-5 * cosd( boxlats(col) ) / r_earth
                    gamma_d = 9.8
                    d_mid(col) = h_scale * log( 1. - f_cor * delta_T_edge(col) / delta_y_edge(col) / ( h_scale * beta * &
                        &( gamma_d + lapsecritcols(col) ) ) )
                    ! d_mid(col) = d_mid(col) * 1.5
                    d_trop(col) = wklm1cols(1,col) / wbrodlmcols(1,col) * Lv / ( cptot(1) * ( gamma_d + lapsecritcols(col) ) )
                    ! meridtransp_edge(0) = delta_T_edge(0) / delta_x_edge(0) * (1.0 - (x_lats(1))**2.0) * d_vl(1)
                    ! if (col < ncols) then
                    !     meridtransp_edge(col) = delta_T_edge(col) / delta_x_edge(col) * (1.0 - (x_lats(col))**2.0) * d_vl(col)
                    ! end if
                    ! meridtransp_edge(0) = 0.0
                    ! meridtransp_edge(ncols) = 0.0
                    ! delta_meridtransp_edge(col) = (meridtransp_edge(col) - meridtransp_edge(col-1))

                    if (lapse_type == 1) then
                        lapsecritcols(col) = lapsecritcols(col) + (max(d_mid(col),d_trop(col))-&
                            &altzmcols(conv_trop_ind(col),col)/1000.) * 0.2
                    end if

                    if (mtranspon == 1) then
                        boxnettotflux(col) = boxnetradflux(col) + meridtransp(col)
                    else
                        boxnettotflux(col) = boxnetradflux(col)
                    end if

                    if (boxnettotflux_prev(col) .ne. 0 .and. boxnettotflux(col) / boxnettotflux_prev(col) < 0.0) then 
                        ! ur_toafnet(col) = ur_toafnet(col) * 2.0
                        ur_toafnet(col) = ur_toafnet(col) * 1.1
                        print*, 'ur_toafnet increased to: ', ur_toafnet(col), 'in col: ', col, boxnettotflux(col),&
                        & boxnettotflux_prev(col)
                    end if
                    boxnettotflux_prev(col) = boxnettotflux(col)


                    tempchanges(col) = (boxnetradflux(col) + meridtransp(col)*ur_mt) / ur_toafnet(col)
                ENDDO


                xglobsum = 0.
                xglobmean = 0.
                cossums = 0.
                do col=1,ncols
                    xglobsum = xglobsum + boxnetradflux(col) * cosd(boxlats(col))
                    cossums = cossums + cosd(boxlats(col))
                enddo
                xglobmean = xglobsum / cossums
                netradflux_globmean = xglobmean

                ! xglobsum = 0.
                ! xglobmean = 0.
                ! cossums = 0.
                ! do col=1,ncols
                !     xglobsum = xglobsum + meridtransp(col) * cosd(boxlats(col))
                !     cossums = cossums + cosd(boxlats(col))
                ! enddo
                ! xglobmean = xglobsum / cossums
                ! meridtransp_globmean = xglobmean
                meridtransp_globmean = sum(meridtransp) / ncols

                xglobsum = 0.
                xglobmean = 0.
                cossums = 0.
                do col=1,ncols
                    xglobsum = xglobsum + boxnettotflux(col) * cosd(boxlats(col))
                    cossums = cossums + cosd(boxlats(col))
                enddo
                xglobmean = xglobsum / cossums
                nettotflux_globmean = xglobmean

                xglobsum = 0.
                xglobmean = 0.
                cossums = 0.
                do col=1,ncols
                    xglobsum = xglobsum + olrcols(col) * cosd(boxlats(col))
                    cossums = cossums + cosd(boxlats(col))
                enddo
                xglobmean = xglobsum / cossums
                olr_globmean = xglobmean

                xglobsum = 0.
                xglobmean = 0.
                cossums = 0.
                do col=1,ncols
                    xglobsum = xglobsum + abs_sw(col) * cosd(boxlats(col))
                    cossums = cossums + cosd(boxlats(col))
                enddo
                xglobmean = xglobsum / cossums
                abs_sw_globmean = xglobmean

                xglobsum = 0.
                xglobmean = 0.
                cossums = 0.
                do col=1,ncols
                    xglobsum = xglobsum + abs_h2o_cols(col) * cosd(boxlats(col))
                    cossums = cossums + cosd(boxlats(col))
                enddo
                xglobmean = xglobsum / cossums
                abs_h2o_globmean = xglobmean

                xglobsum = 0.
                xglobmean = 0.
                cossums = 0.
                do col=1,ncols
                    xglobsum = xglobsum + abs_o3_cols(col) * cosd(boxlats(col))
                    cossums = cossums + cosd(boxlats(col))
                enddo
                xglobmean = xglobsum / cossums
                abs_o3_globmean = xglobmean

                xglobsum = 0.
                xglobmean = 0.
                cossums = 0.
                do col=1,ncols
                    xglobsum = xglobsum + abs_surf_cols(col) * cosd(boxlats(col))
                    cossums = cossums + cosd(boxlats(col))
                enddo
                xglobmean = xglobsum / cossums
                abs_surf_globmean = xglobmean

                xglobsum = 0.
                xglobmean = 0.
                cossums = 0.
                do col=1,ncols
                    xglobsum = xglobsum + sebcols(col) * cosd(boxlats(col))
                    cossums = cossums + cosd(boxlats(col))
                enddo
                xglobmean = xglobsum / cossums
                seb_globmean = xglobmean


                write(*,1106,advance='no') ''

                do col=1,ncols+1
                    if (col < ncols+1) then
                        write(*,1104,advance='no') boxlats(col)
                    else
                        write(*,*) 'AVG'
                    endif
                enddo



                do col=1,ncols
                    if (col < ncols) then
                        write(*,1105,advance='no') '-------------'
                    else
                        write(*,1105) '-------------'
                    endif
                enddo

                write(*,1106,advance='no') 'Box Net Radiative Flux | '
                do col=1,ncols+1
                    if (col < ncols+1) then
                        write(*,1103,advance='no') boxnetradflux(col)
                    else
                        write(*,1103) netradflux_globmean
                    endif
                enddo

                write(*,1106,advance='no') 'Meridional Transport | '
                do col=1,ncols+1
                    if (col < ncols+1) then
                        write(*,1103,advance='no') meridtransp(col)
                    else
                        write(*,1103) meridtransp_globmean
                    endif
                enddo

                write(*,1106,advance='no') 'Box Total Net Flux | '
                do col=1,ncols+1
                    if (col < ncols+1) then
                        write(*,1103,advance='no') boxnettotflux(col)
                    else
                        write(*,1103) nettotflux_globmean
                    endif
                enddo

                write(*,1106,advance='no') 'Box Surface Temperature | '
                do col=1,ncols+1
                    if (col < ncols+1) then
                        write(*,1103,advance='no') tboundmcols(col)
                    else
                        write(*,1103) tglobmean
                    endif
                enddo

                write(*,1106,advance='no') 'Box OLR | '
                do col=1,ncols+1
                    if (col < ncols+1) then
                        write(*,1103,advance='no') olrcols(col)
                    else
                        write(*,1103) olr_globmean
                    endif
                enddo

                write(*,1106,advance='no') 'Box Abs SW | '
                do col=1,ncols+1
                    if (col < ncols+1) then
                        write(*,1103,advance='no') abs_sw(col)
                    else
                        write(*,1103) abs_sw_globmean
                    endif
                enddo

                write(*,1106,advance='no') 'Box Abs SW H2O | '
                do col=1,ncols+1
                    if (col < ncols+1) then
                        write(*,1103,advance='no') abs_h2o_cols(col)
                    else
                        write(*,1103) abs_h2o_globmean
                    endif
                enddo

                write(*,1106,advance='no') 'Box Abs SW O3 | '
                do col=1,ncols+1
                    if (col < ncols+1) then
                        write(*,1103,advance='no') abs_o3_cols(col)
                    else
                        write(*,1103) abs_o3_globmean
                    endif
                enddo

                write(*,1106,advance='no') 'Box Abs SFC | '
                do col=1,ncols+1
                    if (col < ncols+1) then
                        write(*,1103,advance='no') abs_surf_cols(col)
                    else
                        write(*,1103) abs_surf_globmean
                    endif
                enddo

                write(*,1106,advance='no') 'Box SEB | '
                do col=1,ncols+1
                    if (col < ncols+1) then
                        write(*,1103,advance='no') sebcols(col)
                    else
                        write(*,1103) seb_globmean
                    endif
                enddo

                print*, ('----------------------------------------------')
                print*,

                ! toaeqbcheck
                if (maxval(abs(boxnettotflux)) < toa_precision) then                    

                    if(detailprint==1) then
                        print*, "Equilibrium reached in ",nint(j/ur_htr),"days"
                        print*,
                        print*, 'Tsgm: ',tglobmean
                    endif

                    call writeoutputfile


                    ! do i=1,ncols

                    !     write(ccfracsfn,"(A84,I2)") 'Input Distributions/ccfracs col ', col-1      
                    !     write(cctausfn,"(A84,I2)") 'Input Distributions/cctaus col ', col-1   
                    !     write(ccaltsfn,"(A84,I2)") 'Input Distributions/ccalts col ', col-1   

                    !     open(87,file=trim(ccfracsfn),form='formatted')
                    !     open(88,file=trim(cctausfn),form='formatted')
                    !     open(89,file=trim(ccaltsfn),form='formatted')

                    !     open(92,file=('extra_clds'),form='formatted')

                    !     read(92,*) extra_cld_tau
                    !     read(92,*) extra_cld_frac
                    !     read(92,*) extra_cld_alt
                    !     read(92,*) extra_cld_latcol
                    !     read(92,*) extra_cld_cldcol

                    !     close(92)

                    !     do cloudcol = 1,ncloudcols
                    !         read(87,*) ccfracs(cloudcol)
                    !         read(88,*) cctaus(cloudcol)
                    !         read(89,*) ccalts(cloudcol)
                    !     end do

                    !     close(87)
                    !     close(88)
                    !     close(89)

                    !     do cloudcol = 1,ncloudcols
                    !         cloudcolfrac = ccfracs(cloudcol)
                    !         cloudcoltau = cctaus(cloudcol)
                    !         cloudcolalt = ccalts(cloudcol)
                    !         call writecloudstofile
                    !     enddo

                    ! enddo

                    ! do col=1,ncols
                    !     do i=0,nlayersm
                    !         write(50,*) tzmcols(i,col)
                    !     enddo
                    ! enddo

                    EXIT
                endif
            endif
        endif

        stepssinceboxadj = stepssinceboxadj + 1

    enddo !Main do

    !    call printtime('END')
    call printtime

    tot_albedo = 1.0 - tot_sol_abs_lh / sol_inc

    ! if (detailprint==1) then
    !     write(*,'(A,F6.2)') "total LH albedo: ", tot_albedo
    ! endif

    write(62,1002) vol_mixco2,totuflum(OLR_layer)
    write(63,'(E8.2)') vol_mixco2

    do i=1,nlayersm
        write(63,'(F10.4)') totuflum(i)
    enddo

    write(63,*)

    !Write the non-layer variables to file
    ! write(50,1004) massatmo,totmolec,undrelax,lapsecrit,mo2,mixo2,timesteps,ipe,co2_inv,&
    ! air_inv,solar_constant,tot_albedo,sol_inc,surfacep,gravity,nlayersm,conv_trop_ind,&
    ! nlayersm 

    !Write general information for every layer to file
    ! do i=nlayersm,0,-1
    !     write(50,1000) i,altzm(i)/1000,pzm(i),pavelm(i),tzm(i),tavelm(i),totuflum(i),totdflum(i),fnetm(i),htrm(i),mixco2(i),&
    !     mixh2o(i),mixo3(i),wklm(1,i),wklm(2,i),htrh2o(i),rel_hum(i),mperlayr(i),wbrodlm(i),cmh2o(i),totcmh2o(i),solabsh2o(i),&
    !     laysolabsh2o(i),cmo3(i),totcmo3(i),solabso3vis(i),solabso3uv(i),solabso3tot(i),laysolabso3tot(i),htro3(i),es(i),&
    !     abspn(i),A_oz_l(i),htrlh(i),htro3_lh(i),u(i),cptot(i),wklm(3,i),band_flux_up(1,i),band_flux_up(2,i),band_flux_up(3,i),&
    !     band_flux_up(4,i),band_flux_up(5,i),band_flux_up(6,i),band_flux_up(7,i),band_flux_up(8,i),band_flux_up(9,i),&
    !     band_flux_up(10,i),band_flux_up(11,i),band_flux_up(12,i),band_flux_up(13,i),band_flux_up(14,i),band_flux_up(15,i),&
    !     band_flux_up(16,i),band_flux_down(1,i),band_flux_down(2,i),band_flux_down(3,i),band_flux_down(4,i),band_flux_down(5,i),&
    !     band_flux_down(6,i),band_flux_down(7,i),band_flux_down(8,i),band_flux_down(9,i),band_flux_down(10,i),band_flux_down(11,i),&
    !     band_flux_down(12,i),band_flux_down(13,i),band_flux_down(14,i),band_flux_down(15,i),band_flux_down(16,i), malr(i),&
    !     rsp_tot(i), cptot(i), theta(i),kappa(i),montgomery(i),exner(i),geopotential(i)
    ! enddo

    do i=nlayersm,1,-1
        write(61,1003) wklm(1,i),wklm(2,i),wklm(3,i),wklm(4,i),wklm(5,i),wklm(6,i),wklm(7,i),wklm(8,i),wklm(9,i)
    enddo

    !if (detailprint==1) then
    !        print*,
    !        print*, '--------------------------------------------'
    !        
    !        print*,
    !        write(*,'(A,F6.2)') "OLR: ",totuflum(OLR_layer)
    !        write(*,'(A,F6.2)') "Tg: ",tzm(0)
    !        print*,
    !        
    !        
    !        
    !        write(*,'(A,F6.2)') 'Absorption by H2O:      ',abs_h2o
    !        write(*,'(A,F6.2)') 'Absorption by O3:       ',abs_o3
    !        write(*,'(A,F6.2)') 'Absorption by surface:  ',abs_surf
    !        write(*,'(A,F6.2)') 'Total solar absorption: ',abs_h2o+abs_o3+abs_surf
    !        
    !        print*,
    !        
    !        
    !        print*, '--------------------------------------------'
    !        print*,
    !        
    !        print*,
    !endif

    print *, char(7)

    1000 FORMAT (I5,80E16.6)
    1001 FORMAT (I5,F8.3,7E12.3)
    1002 FORMAT (E8.2,F10.4)
    1003 FORMAT (9E12.3)
    1004 FORMAT (10E16.4,2I5,7E20.4,I7,3E16.4,I5,E16.4,2I5)
    1005 FORMAT (I5,2F12.4,E16.4)
    1080 FORMAT (6E24.16)
    1081 FORMAT (17E24.16)
    1082 FORMAT (39E11.4)
    1083 FORMAT (I3, I3, F5.3, F4.2, F4.2, A1)
    1100 FORMAT (i4,'-',i2.2,'-',i2.2,' ',i2.2,"_",i2.2,'_',i2.2)
    1101 FORMAT ('Max htr in column',I1)
    1102 FORMAT (I4,F6.2,' | ')
    1103 FORMAT (F10.2,' | ')
    1104 FORMAT ('LAT ',F6.2,' | ')
    1105 FORMAT (A13) 
    1106 FORMAT (A27)
    1107 FORMAT (I3,' ')
    1109 FORMAT (' ts ')
    1110 FORMAT (F12.4)

    ! close(50)
    close(51)
    close(52)
    close(53)
    close(56)
    close(57)
    close(58)
    close(60)
    close(61)
    close(62)
    close(63)
    close(74)
    close(76)
    close(81)
    close(82)
    close(83)
    close(84)
    close(85)
    close(90)
    close(92)

END subroutine wrapper
