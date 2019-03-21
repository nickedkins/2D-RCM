! this program written by Roger Davies 27/8/2012ff
! calls RRTM
PROGRAM myprrtm
!use timeprogram
use VARIABLES
use MYSUBS
implicit none
!integer	::i


!call starttiming
	
	call wrapper
!	call rrtm



!call elapsedtiming
!call printtiming

! WRITE(*,9900)
! WRITE(*,9901)

!DO 3000 I = nlayersm, 0, -1
!!	write(*,9958) i,altzm(i)*0.001,pzm(i),tzm(i),totuflum(I),totdflum(I),fnetm(i),htrm(I)
!3000    CONTINUE
! 9958 FORMAT(1X,I3,F7.1,F10.2,F8.1,F10.3,F15.1,f12.1,f17.2)
! 9900 FORMAT(1X,'LEVEL alt   PRESSURE   temp UPWARD FLUX   DOWNWARD FLUX    NETFLUX        HEATING RATE')
! 9901 FORMAT(1X,'       km     mb         K      W/m2          W/m2           W/m2          degree/day')
! 9902 format(A13,F8.1,A9)

END program myprrtm
