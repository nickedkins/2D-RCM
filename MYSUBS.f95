MODULE MYSUBS

    !    use VARIABLES


    contains

    !         subroutine printtime(whencalled)
    subroutine printtime
	!nnn
        use VARIABLES
        ! character(8)  :: date
        ! character(10) :: time
        ! character(5)  :: zone
        ! integer,dimension(8) :: values
        ! character(len=*) :: whencalled

        call date_and_time(date,time,zone,values)
        call date_and_time(DATE=date,ZONE=zone)
        call date_and_time(TIME=time)
        call date_and_time(VALUES=values)
        write(*,*) trim(whencalled), ' Date: ',values(1),values(2),values(3), 'Time: ',values(5),values(6),values(7)
    end subroutine printtime

    ! subroutine levconvect(col,tzm,altzm,altlaym,lapsecrit,fixed_trop,conv,convecttype,malr)
    subroutine levconvect
        use VARIABLES
        ! real,dimension(0:MAXLAYM) :: altzm,tzm
        ! real,dimension(maxlaym) :: malr,altlaym
        ! real :: lapsecrit,fixed_trop
        ! integer :: convecttype,col,i
        ! real,dimension(maxlaym,ncols) :: conv

        select case(convecttype)
        case(0) !critical lapse rate
            do i=1,nlayersm
                if( (tzm(i) - tzm(i-1))/(altzm(i)-altzm(i-1))*1000.0 < lapsecrit .or. pzm(i) > fixed_trop(col)) then
                    tzm(i) = tzm(i-1) +lapsecrit*(altzm(i)-altzm(i-1))/1000.0
                    if (tzm(i) < t_min) tzm(i) = t_min
                    conv(i,col) = 1
                else
                    conv(i,col) = 0
                endif
            end do
        case(2) !moist adiabatic lapse rate
            do i=1,nlayersm
                if( (tavelm(i) - tzm(i-1))/(altlaym(i)-altzm(i-1))*1000.0 < (malr(i)*1000.0) .or. &
                    pavelm(i) > fixed_trop(col)) then
                    tzm(i) = tzm(i-1) + malr(i) * 1000.0*(altzm(i)-altzm(i-1))/1000.0
                    if (tzm(i) < t_min) tzm(i) = t_min
                    conv(i,col) = 1
                else
                    conv(i,col) = 0
                endif
            end do
        end select
    end subroutine levconvect

    ! subroutine calcswhtr(altzm,pzm,tzm,mixh2o,mmwtot,tau_cld,wklm,htrlh,htrlh_tot,htro3_lh,htro3_lh_tot,&
    !     abs_surf,R_g,tot_sol_abs_lh,abs_h2o,abs_o3,sol_inc,cloudcolp,cloudcoltau)
    subroutine calcswhtr
        use VARIABLES
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

        real :: sw_eps = 1.0e-32


        Ag1 = 0.0
        Ag2 = 0.0
        w_t = 0.0
        ylh = 0.0
        A_wv = 0.0
        Rbar_r = 0.0
        rbar_a = 0.
        rbarbarstar_a = 0.
        tau = 0.

        abs_surf = 0.
        abs_h2o = 0.
        abs_o3 = 0.
        tot_sol_abs_lh = 0.

        ! mu_0 = 0.5
        mag = 35.0/((1224.0*mu_0**2.0 + 1.0)**0.5)
        p_0 = 1013.0 !in hPa here
        t_0 = 273.0
        asym = 0.85

        pk(1) = 0.647
        pk(2) = 0.0698
        pk(3) = 0.1443
        pk(4) = 0.0584
        pk(5) = 0.0335
        pk(6) = 0.0225
        pk(7) = 0.0158
        pk(8) = 0.0087

        k_n(1) = 0.00004
        k_n(2) = 0.002
        k_n(3) = 0.035
        k_n(4) = 0.377
        k_n(5) = 1.95
        k_n(6) = 9.40
        k_n(7) = 44.6
        k_n(8) = 190.


        do i=1,nlayersm
            ! rho_w(i) = mixh2o(i) * (18.01/1000.0) / mmwtot * 1000.0 !times 1000 to convert to g/kg instead of kg/kg
            rho_w(i) = mixh2o(i) * (18.01/1000.0) / mmwtot  !times 1000 to convert to g/kg instead of kg/kg ! but I've removed that factor -- might be wrong
        enddo

        ! zlh=altzm(0:nlayersm)/1000.0

        ! plh=pzm(0:nlayersm)
        ! tlh=tzm(0:nlayersm)

        ! plh = plh*100.0

        ! zlh =  zlh(size(zlh):1:-1)
        ! plh =  plh(size(plh):1:-1)
        ! tlh =  tlh(size(tlh):1:-1)
        ! rho_w =  rho_w(size(rho_w):1:-1)
        ! rho_w(1) = rho_w(2)

        dk(1) = 4.0e-5

        do i=2,size(k_n)-1
            dk(i) = k_n(i) - k_n(i-1)
        enddo

        do i=2,nlayersm
            delta_z(i) = altzm(i) - altzm(i-1)
        enddo

        do i=1,nlayersm-1
            wl(i) = (rho_w(i)/1000.0*delta_z(i+1)*1000.0)*( (pzm(i)/p_0)**1.0 )*((t_0/tzm(i))**0.5)     
        enddo

        ! cloudindex = int(minloc(abs(pzm - cloudcolp),dim=1))
        cloudindex = int(minloc(abs(altzm/1000.0 - cloudcolalt),dim=1))
        tau_cld = 0.
        tau_cld(int(cloudindex)) = cloudcoltau

        tau = 0.

        temp_plh = plh

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        ! NJE normal method of including tau_cld in total tau calc
        ! do i=1,nlayersm
        !     do n=2,8
        !         tau(i,n) = tau_cld(i) + k_n(n) * wl(i)
        !     enddo
        ! enddo

        ! plh = temp_plh

        ! ssa = 0.

        ! temp = tau

        ! do i=1,nlayersm
        !     do n=2,8
        !         ssa(i,n) = tau_cld(i)/tau(i,n)
        !         if (ssa(i,n) .ge. 1.0-sw_eps) ssa(i,n) = 1.0-sw_eps
        !     enddo
        ! enddo

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! NJE remove the cloud optical thickness from the total tau calc, but leave it in for the SSA calc
        do i=1,nlayersm
            do n=2,8
                tau(i,n) = k_n(n) * wl(i)
            enddo
        enddo

        plh = temp_plh

        ssa = 0.

        temp = tau

        do i=1,nlayersm
            do n=2,8
                if (tau_cld(i) + tau(i,n) .ne. 0) ssa(i,n) = tau_cld(i) / ( tau_cld(i) + tau(i,n) )
                if (ssa(i,n) .ge. 1.0-sw_eps) ssa(i,n) = 1.0-sw_eps
            enddo
        enddo

        !!!!!!!!!!!!!!!!!!!!!!!!!!!        

        do i=1,nlayersm
            do n=2,8
                ulh(i,n) = ((1.0-asym*ssa(i,n))/(1.0-ssa(i,n)))**0.5
            enddo
        enddo

        do i=1,nlayersm
            do n=2,8
                trlh(i,n) = ((3.0*(1.0-ssa(i,n)*(1.0-asym*ssa(i,n))))**0.5)*tau(i,n)
            enddo
        enddo

        refl_lay_ind = 0

        do i=1,nlayersm-2
            ! if (tau_cld(i) .ge. 1e-4 .and. tau_cld(i) < 100.0 .and. altzm(i)/1000.0 < 10.0) then
            if (tau_cld(i) .ge. 1e-4) then
                refl_lay_ind = i
                exit
            endif
        enddo

        refl_lay_z = altzm(refl_lay_ind)

        cld = 0
        if (refl_lay_ind .ne. 0) cld = 1

        Translh = 0.

        do i=1,nlayersm
            do n=2,8
                if (tau_cld(i) .ge. 1e-4 ) then
                    Translh(i,n) = 4.0*ulh(i,n) / ( (ulh(i,n)+1.0)**2.0*exp(trlh(i,n)) - ((ulh(i,n)-1.0)**2.0*&
                        exp(-trlh(i,n) )))
                elseif (i < refl_lay_ind) then
                    Translh(i,n) = exp(-mag*tau(i,n))
                else
                    Translh(i,n) = exp(-5.0/3.0*tau(i,n))
                endif
            enddo
        enddo

        Rlh = 0.

        do i=1,nlayersm
            do n=2,8
                if(tau_cld(i) .ge. 1e-4) then
                    Rlh(i,n) = ((ulh(i,n)+1.0)*(ulh(i,n)-1.0)*(exp(trlh(i,n)) - exp(-trlh(i,n)) ))/&
                    ( (ulh(i,n)+1.0)**2.0* exp(trlh(i,n)) - (ulh(i,n)-1.0)**2.0*exp(-trlh(i,n))  )
                else
                    Rlh(i,n) = 0.0
                endif
            enddo
        enddo

        do i=1,nlayersm
            do n=2,8
                if (tau_cld(i) .ge. 1e-4) then
                    Rl_star(i,n) = Rlh(i,n)
                else
                    Rl_star(i,n) = 0.0
                endif
            enddo
        enddo

        Tl_star = 0.0

        do i=1,nlayersm
            do n=2,8
                if (tau_cld(i) .ge. 1e-4) then
                    Tl_star(i,n) = Translh(i,n)
                else
                    Tl_star(i,n) = exp(-5.0/3.0*tau(i,n))
                endif
            enddo
        enddo

        Llh = nlayersm - 1

        R_1l = 0.
        T_1l = 0.
        R_1l_star = 0.
        T_1l_star = 0.

        do n=2,8
            R_1l(1,n) = Rlh(1,n)
        enddo

        !adding going down
        do l=3,Llh+1 
            i = l-1
            do n=2,8
                if (l==3) then
                    Ra = Rlh(1,n)
                    Ta = Translh(1,n)
                    Ra_star = Rl_star(1,n)
                    Ta_star = Tl_star(1,n)
                else
                    Ra = R_1l(i-1,n)
                    Ta = T_1l(i-1,n)
                    Ra_star = R_1l_star(i-1,n)
                    Ta_star = T_1l_star(i-1,n)
                endif
                Rb = Rlh(i,n)
                Tb = Translh(i,n)
                Rb_star = Rl_star(i,n)
                Tb_star = Tl_star(i,n)
                if (Ra .ge. 1.0-sw_eps) Ra = 1.0-sw_eps
                if (Rb .ge. 1.0-sw_eps) Rb = 1.0-sw_eps
                if (Ta .ge. 1.0-sw_eps) Ta = 1.0-sw_eps
                if (Tb .ge. 1.0-sw_eps) Tb = 1.0-sw_eps
                if (Ra_star .ge. 1.0-sw_eps) Ra_star = 1.0-sw_eps
                if (Ta_star .ge. 1.0-sw_eps) Ta_star = 1.0-sw_eps
                if (Rb_star .ge. 1.0-sw_eps) Rb_star = 1.0-sw_eps
                if (Tb_star .ge. 1.0-sw_eps) Tb_star = 1.0-sw_eps
                if (Ra .le. sw_eps) Ra = sw_eps
                if (Rb .le. sw_eps) Rb = sw_eps
                if (Ta .le. sw_eps) Ta = sw_eps
                if (Tb .le. sw_eps) Tb = sw_eps
                if (Ra_star .le. sw_eps) Ra_star = sw_eps
                if (Ta_star .le. sw_eps) Ta_star = sw_eps
                if (Rb_star .le. sw_eps) Rb_star = sw_eps
                if (Tb_star .le. sw_eps) Tb_star = sw_eps
                R_1l(i,n) = Ra + Ta*Rb*Ta_star/(1.0-Ra_star*Rb)
                T_1l(i,n) = Ta*Tb/(1.0-Ra_star*Rb)
                R_1l_star(i,n) = Rb_star + Tb_star*Ra_star*Tb/(1.0-Ra_star*Rb)
                T_1l_star(i,n) = Tb_star*Ta_star/(1.0-Ra_star*Rb)
            enddo
        enddo

        R_lg = 0.0
        Ul= 0.0
        Dl = 0.0
        A_1l = 0.0
        Al = 0.0

        do n=2,8
            R_lg(Llh,n) = R_g
            T_1l(1,n) = Translh(1,n)
            R_1l(1,n) = Rlh(1,n)
            T_1l_star(1,n) = Tl_star(1,n)
            R_1l_star(1,n) = Rl_star(1,n)
        enddo

        !Adding going up
        do l = Llh+1,2,-1
            i=l-1
            do n=2,8
                if (i==Llh) then
                    Rb = R_g
                else
                    Rb = R_lg(i+1,n)
                endif
                Ra = Rlh(i,n)
                Ta = Translh(i,n)
                Ra_star = Rl_star(i,n)
                Ta_star = Tl_star(i,n)
                if (Ra .ge. 1.0-sw_eps) Ra = 1.0-sw_eps
                if (Rb .ge. 1.0-sw_eps) Rb = 1.0-sw_eps
                if (Ta .ge. 1.0-sw_eps) Ta = 1.0-sw_eps
                if (Tb .ge. 1.0-sw_eps) Tb = 1.0-sw_eps
                if (Ra_star .ge. 1.0-sw_eps) Ra_star = 1.0-sw_eps
                if (Ta_star .ge. 1.0-sw_eps) Ta_star = 1.0-sw_eps
                if (Rb_star .ge. 1.0-sw_eps) Rb_star = 1.0-sw_eps
                if (Tb_star .ge. 1.0-sw_eps) Tb_star = 1.0-sw_eps
                if (Ra .le. sw_eps) Ra = sw_eps
                if (Rb .le. sw_eps) Rb = sw_eps
                if (Ta .le. sw_eps) Ta = sw_eps
                if (Tb .le. sw_eps) Tb = sw_eps
                if (Ra_star .le. sw_eps) Ra_star = sw_eps
                if (Ta_star .le. sw_eps) Ta_star = sw_eps
                if (Rb_star .le. sw_eps) Rb_star = sw_eps
                if (Tb_star .le. sw_eps) Tb_star = sw_eps
                R_lg(i,n) = Ra + Ta*Rb*Ta_star/(1.0-Ra_star*Rb)
            enddo
        enddo

        do l=Llh,2,-1
            i=l-1
            do n=2,8
                Ul(i,n) = T_1l(i,n)*R_lg(i+1,n)/(1.0-R_1l_star(i,n)*R_lg(i+1,n))
                Dl(i,n) = T_1l(i,n)/(1.0-R_1l_star(i,n)*R_lg(i+1,n))
                A_1l(i,n) = pk(n)*(1.0-R_lg(1,n)+Ul(i,n)-Dl(i,n))
            enddo
        enddo

        do l=2,Llh+1
            i = l-1
            do n=2,8
                if(l==1) then
                    Al(i,n) = A_1l(1,n)
                else
                    Al(i,n) = A_1l(i,n)-A_1l(i-1,n)
                endif
                if (Al(i,n) .ne. Al(i,n)) Al(i,n) = 0.0 !Check for NaNs, since NaN is .ne. anything, even itself
            enddo
        enddo

        abspn = 0.0

        do l = 2,Llh
            i = l-1
            abspn(i) = sum(Al(i,2:8))
        enddo



        delta_p = 0.0
        htrlh = 0.0

        do i=1,nlayersm
            delta_p(i) = pzm(i-1) - pzm(i)
        enddo

        do l=2,Llh+1
            i=l-1
            htrlh(i) = sol_inc * abspn(i) * gravity / (delta_p(i)*100.0 * cplh)*60.0*60.0*24.0 !cplh is specific heat constant pressure dry air
            ! htrlh(i) = sol_inc / 3.0 * abspn(i) * gravity / (delta_p(i)*100.0 * cplh)*60.0*60.0*24.0 !NJE
            if (htrlh(i) .ne. htrlh(i)) htrlh(i) = 0.0            
        enddo


        !!!!!!!!!!!!!!!!!!!!!!!
        !Lacis and Hansen Ozone

        ! zlh =  zlh(nlayersm:1:-1)
        ! plh =  plh(size(plh):1:-1)
        ! delta_p = delta_p(nlayersm:1:-1)
        ! tlh =  tlh(size(tlh):1:-1)
        ! rho_w =  rho_w(size(rho_w):1:-1)
        ! htrlh = htrlh(size(htrlh):1:-1)

        hc = 8.0

        refl_lay_ind = 1

        do i=nlayersm-1,1,-1
            if (tau_cld(i) .ge. 1.0) then
                refl_lay_ind = i
                exit
            endif
        enddo

        tau_c = cloudcoltau

        u = 0.0
        x_o3 = 0.0
        A_oz_vis_x = 0.0
        xstar = 0.0
        A_oz_vis_xstar = 0.0
        A_oz_uv_x = 0.0
        A_oz_uv_xstar = 0.0
        A_oz_x = 0.0
        A_oz_xstar = 0.0
        Rbar_a = 0.0
        Rbarbarstar_a = 0.0
        Rbar = 0.0
        Rbarbar_a = 0.0

        do i=nlayersm-1,1,-1
            u(i) = wklm(3,i) / (2.6e19) + u(i+1)
            x_o3(i) = u(i) * mag
        enddo

        do i=1,nlayersm
            xstar(i) = u(refl_lay_ind)*mag + 1.9*(u(refl_lay_ind)-u(i))
            A_oz_vis_x(i) = 0.02188*x_o3(i)/(1.0+0.042*x_o3(i)+0.000323*x_o3(i)**2.0)
            A_oz_vis_xstar(i) = 0.02188*xstar(i)/(1.0+0.042*xstar(i)+0.000323*xstar(i)**2.0)
            A_oz_uv_x(i) = 1.082*x_o3(i)/((1.0+138.6*x_o3(i))**0.805) + 0.0658*x_o3(i)/(1.0+(103.6*x_o3(i))**3.0)
            A_oz_uv_xstar(i) = 1.082*xstar(i)/((1.0+138.6*xstar(i))**0.805) + 0.0658*xstar(i)/&
            (1.0+(103.6*xstar(i))**3.0)
            A_oz_x(i) = A_oz_vis_x(i) + A_oz_uv_x(i)
            A_oz_xstar(i) = A_oz_vis_xstar(i) + A_oz_uv_xstar(i)
        enddo

        if(cld==0) then
            Rbar_a = 0.219/(1+0.816*mu_0)
        else
            Rbar_a = 0.13*tau_c/(1.0+0.13*tau_c)
        endif

        if(cld==0) then
            Rbarbarstar_a = 0.144
        else
            Rbarbarstar_a = 0.13*tau_c/(1.0+0.13*tau_c)
        endif

        Rbar = Rbar_a + (1.0-Rbar_a)*(1.0-Rbarbarstar_a)*R_g/(1.0-Rbarbarstar_a*R_g)

        ! do i=1,nlayersm - 1
        do i=2,nlayersm - 1
            ! A_oz_l(i) =   mu_0* ( A_oz_x(i) - A_oz_x(i+1) + Rbar * (A_oz_xstar(i+1) - A_oz_xstar(i)) ) 
            A_oz_l(i) =   mu_0* ( A_oz_x(i-1) - A_oz_x(i) + Rbar * (A_oz_xstar(i) - A_oz_xstar(i-1)) ) 
            htro3_lh(i) = sol_inc  *  A_oz_l(i) * gravity / (delta_p(i)*100.0 * cplh)*60.0*60.0*24.0 !May need to add mu_0 and /2.0 (for night/day) in here
            htro3_lh(i) = htro3_lh(i) * 3.14 !  Check the validity of this line and the one above - do I have the correct incident flux?
            if (i < 4) htro3_lh(i) = 0.0
            if (htro3_lh(i) .ne. htro3_lh(i)) htro3_lh(i) = 0.0

        enddo



        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !Surface Absorption

        Ag1 = 0.0
        Ag2 = 0.0

        if (cld == 0) then
            w_t = wl(1)
            ylh = mag*w_t
            A_wv = 2.9*ylh / ((1.0+141.5*ylh)**0.635 + 5.925*ylh )
            Ag1 = mu_0*(0.353 - A_wv)*(1.0-R_g)
            Rbar_r = 0.28/(1.0+6.43*mu_0)
            Ag2 = mu_0*(0.647 - Rbar_r - A_oz_x(1))*(1.0-R_g)/(1.0-0.0685*R_g)  
        else
            Ag1 = 0.0
            Ag2 = 0.0
            do n=2,8
                if ( T_1l(nlayersm-2,n) .ne. T_1l(nlayersm-2,n)) T_1l(nlayersm-2,n) = 0.0
                ! Ag1 = Ag1 + T_1l(nlayersm-2,n)*pk(n)*(1.0-R_g)
                Ag1 = Ag1 + T_1l(1,n)*pk(n)*(1.0-R_g)
            enddo
            Ag2 = mu_0*(0.647-A_oz_x(refl_lay_ind))*(1.0-Rbar_a)*(1.0-R_g)/(1.0-Rbarbarstar_a*R_g)
        endif
        
        abs_surf_lh = 0.

        abs_surf_lh = (Ag1+Ag2)*sol_inc*2.0 !Unsure about that factor of 2.0 nje
        ! abs_surf_lh = (Ag1+Ag2)*sol_inc

        abs_surf_lhcols(col) = abs_surf_lh

        if (abs_surf_lh .ne. abs_surf_lh) then
            print*, "abs_surf is NaN:"
            call EXIT(11)
        endif

        abs_h2o = sum(abspn(:nlayersm-1))*sol_inc !total absorption of SW by H2O  
        abs_o3 = sum(A_oz_l(2:nlayersm))*sol_inc ! total absorption of SW by O3
        tot_sol_abs_lh = abs_h2o + abs_o3 + abs_surf_lh ! Total absorption of SW 

    end subroutine calcswhtr

    ! subroutine createcloudfile(cloudp,cloudcoltau)
    subroutine createcloudfile(cca,cct)
        use VARIABLES
        ! real,dimension(maxlaym) :: fracs,tau_cld
        ! integer :: irdcld,i
        INTEGER :: cloudindex
        real :: cloudalt

        cloudcolalt = cca
        cloudcoltau = cct

        irdcld = 25
        ! OPEN(IRDCLD,FILE='/Users/nickedkins/Dropbox/2D RCM Archive/2018-08-13/My IN_CLD_RRTM',FORM='FORMATTED')
        OPEN(IRDCLD,FILE='My IN_CLD_RRTM',FORM='FORMATTED')
        tau_cld = 0.
        fracs = 0.
        cloudindex = minloc(abs(altzm/1000.0 - cloudcolalt),dim=1)
        tau_cld(cloudindex) = cloudcoltau
        fracs(cloudindex) = 0.8

        if (col == extra_cld_latcol .and. cloudcol == extra_cld_cldcol ) then

            extra_cloudindex = int(minloc(abs(altzm/1000.0 - extra_cld_alt),dim=1))
            tau_cld(extra_cloudindex) = extra_cld_tau
            fracs(extra_cloudindex) = extra_cld_frac

        endif

        write(irdcld,*) "    0    1    1"



        do i=1,nlayersm
            if (tau_cld(nlayersm+1-i) < 10.0) then
                write(irdcld,'(A2,I3,F8.4,F8.4)') 'A  ', i, fracs(i), tau_cld(i)
            elseif (tau_cld(nlayersm+1-i) < 100.0) then 
                write(irdcld,'(A2,I3,F8.4,F8.3)') 'A  ', i, fracs(i), tau_cld(i)
            else 
                write(irdcld,'(A2,I3,F8.4,F8.2)') 'A  ', i, fracs(i), tau_cld(i)
            endif
        enddo

        write(irdcld,'(A)') '%'
        close(irdcld)

    end subroutine createcloudfile      

    subroutine writecloudstofile

        use VARIABLES
        INTEGER :: cloudindex
        real :: cloudalt

        tau_cld = 0.
        fracs = 0.
        cloudindex = minloc(abs(altzm/1000.0 - cloudcolalt),dim=1)
        tau_cld(cloudindex) = cloudcoltau
        fracs(cloudindex) = 0.8

        if (col == extra_cld_latcol .and. cloudcol == extra_cld_cldcol ) then

            extra_cloudindex = int(minloc(abs(altzm/1000.0 - extra_cld_alt),dim=1))
            tau_cld(extra_cloudindex) = extra_cld_tau
            fracs(extra_cloudindex) = extra_cld_frac

        endif

        do i=1,nlayersm
            write(50,*) fracs(i)
        enddo
        
        do i=1,nlayersm
            write(50,*) tau_cld(i)
        enddo

    end subroutine writecloudstofile

    subroutine writeoutputfile

        use VARIABLES

        do col=1,ncols
            do i=0,nlayersm
                write(50,*) tzmcols(i,col)
            enddo
        enddo

        do col=1,ncols
            do i=0,nlayersm
                write(50,*) altzmcols(i,col)
            enddo
        enddo

        do col=1,ncols
            do i=0,nlayersm
                write(50,*) totuflumcols(i,col)
            enddo
        enddo

        do col=1,ncols
            do i=0,nlayersm
                write(50,*) totdflumcols(i,col)
            enddo
        enddo

        do col=1,ncols
            do i=0,nlayersm
                write(50,*) htrmcols(i,col)
            enddo
        enddo

        do col=1,ncols
            do i=1,nlayersm
                write(50,*) htrh2ocols(i,col)
            enddo
        enddo

        do col=1,ncols
            do i=1,nlayersm
                write(50,*) htro3cols(i,col)
            enddo
        enddo

        do col=1,ncols
            do i=1,nlayersm
                write(50,*) wklm1cols(i,col)
            enddo
        enddo

        do col=1,ncols
            do i=1,nlayersm
                write(50,*) wklm2cols(i,col)
            enddo
        enddo

        do col=1,ncols
            do i=1,nlayersm
                write(50,*) wklm3cols(i,col)
            enddo
        enddo

        do col=1,ncols
            do i=1,nlayersm
                write(50,*) wbrodlmcols(i,col)
            enddo
        enddo

        do col=1,ncols
            do i=1,nlayersm
                write(50,*) abspncols(i,col)
            enddo
        enddo

        do col=1,ncols
            do i=1,nlayersm
                write(50,*) A_oz_lcols(i,col)
            enddo
        enddo

        do col=1,ncols
            do i=1,nlayersm
                write(50,*) abs_surf_lhcols(col)
            enddo
        enddo

        do col=1,ncols
            do i=1,nlayersm
                write(50,*) altlaymcols(i,col)
            enddo
        enddo

        do col=1,ncols
            do i=0,nlayersm
                write(50,*) pzmcols(i,col)
            enddo
        enddo

        do col=1,ncols
            do i=1,nlayersm
                write(50,*) pavelmcols(i,col)
            enddo
        enddo

        do col=1,ncols
            do i=1,nlayersm
                write(50,*) tavelmcols(i,col)
            enddo
        enddo

        do col=1,ncols
            do i=1,nlayersm
                write(50,*) tboundmcols(col)
            enddo
        enddo

        do col=1,ncols
            do i=1,nlayersm
                write(50,*) R_gcols(col)
            enddo
        enddo

        do col=1,ncols
            do i=1,nlayersm
                write(50,*) boxlats(col)
            enddo
        enddo


    end subroutine writeoutputfile
    
!     subroutine add_seb_to_tboundm
    
!       use VARIABLES

!             lhf = 10.
!             shf = 0.
!             ! sebfac = 0.01
!             seb = 0.
! !
!             seb = (totdflum(0) + abs_surf_lh - totuflum(0) - lhf - shf)
! !            tboundm = tboundm + seb * sebfac
            
!             tboundm = ((totdflum(0) + abs_surf_lh - lhf - shf)/5.67e-8)**0.25
            
! !            print*, ('----------------------------------------------')
! !                print*,
! !            print*, 'seb, totdflum(0),totuflum(0),abs_surf_lh, lhf ', seb, totdflum(0),totuflum(0),abs_surf_lh, lhf
! !            print*, 'tboundm: ', tboundm
! !            print*, ('----------------------------------------------')
! !            print*,         
!             tboundmcols(col-1) = tboundm
!             tzm(0) = tboundm
    
!     end subroutine add_seb_to_tboundm

END MODULE MYSUBS
