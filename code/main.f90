      include "/home1/ccwang/workspace/general_codes/summary.f90"
	include "par_mod.f90"
	include "trendpar_mod.f90"
	include "tpar_mod.f90" 
	include "tparc_mod.f90" 
	include "ttpar_mod.f90" 

	
	program main_program
	use par_mod
	use tpar_mod
	use tparc_mod
        use ttpar_mod
	use trendpar_mod
	use gnuplot_int
        implicit none

        integer:: nobs,r,ierr,i,nins,ntot,bins,nraw
        integer,allocatable :: obs(:),obsraw(:)
        real*8 :: para_par(3),para_tpar(6),para_ttpar(10),para_trendpar(7)
        character(len=20) :: fn

        real*8 :: nllk,llk,objgrd(10),asycov(10,10),asycov_par(3,3),asycov_trendpar(7,7),asycov_tpar(6,6)
	
	!call test_asycov_tpar
	!call simulation_threshold_tpar
	
	!call test_asycov_tpar
	!stop
	!fn="ntrans5.dat"
	!fn="coli.txt"
	!nins=120
	!ntot=125
	
	!fn="ntrans5t.dat"
	!nins=920
	!ntot=1196
	!prd=92

	
	!fn="c2p3.dat"
	!nraw=32767
	!bins=40
	!ntot=nraw/bins
	!nins=ntot-100

	
	!fn="ntran1.dat"
	!nins=460
	!ntot=460

	fn="earthquakes.dat"
	nraw=111
	bins=1
	nins=100
	ntot=111

	!!sunspot data
	!fn="sunspots.dat"
	!nins=440
	!ntot=459
	!prd=11

        allocate(obsraw(nraw),obs(ntot))

        open(unit=10,file=fn,status="old",iostat=ierr)
        if(ierr.ne.0) then
                print *, "Filename error."
		stop
        end if

        read(10,*) obsraw
        close(10)
	!call gplot1(dble(obsraw))
	do i=1,ntot
		obs(i)=sum(obsraw(1+(i-1)*bins:i*bins))
	end do
	
	call gplot1(dble(obs))
	!asycov=0.0d0
	!para_par=(/0.2d0,0.20d0,0.2d0/)
!
	call epa11_by_R(para_par,obs(1:nins),llk)
	write(*,'(A,3F10.4)') "Estimated parameters by R:",para_par
!pause
	write(*,'(A,F10.4)') "Estimated LLK by R (averaged):",llk/nins
	!note that epa11 used another method to calculate llk, to ensure consistency
	! calculate by our subroutine
	call nllk_par(para_par,obs(1:nins),llk)
	llk=-llk
	write(*,'(A,F10.4)') "Estimated LLK by fortran (averaged):",llk
!pause

!	call epar_wg(para_par,obs(1:nins),llk)
	call epar(para_par,obs(1:nins),llk)
        call dc_par(para_par,obs(1:nins),asycov_par)
	print *, "Summary of results: PoiAR:"
	write(*,'(A)') "Estimated parameters and s.d."
	write(*,'(2F10.4)') (para_par(i),sqrt(asycov_par(i,i)),i=1,3)
!!pause
	print *, "averaged LLK=",llk
	print *, "AIC=", -2.0d0*llk*(nins-1)+2.0d0*3
	print *, "BIC=", -2.0d0*llk*(nins-1)+log(dble(nins)-1)*3
	call coverage_par(para_par,nins,obs)
pause

	para_tpar(1:3)=para_par
	para_tpar(4:6)=para_tpar(1:3)
	call etpar(para_tpar,r,obs(1:nins),llk,'U')
        call dc_tpar(para_tpar,r,obs,asycov_tpar)
	print *, "===Summary of results: TPoiAR==="
	print *, "parameter estimates:"
	write( *,'(A,I0)') "threshold=", r
	write(*,'(A)') "Estimated parameters and s.d."
	write(*,'(2F10.4)') (para_tpar(i),sqrt(asycov_tpar(i,i)),i=1,6)
	print *, "LLK=",llk
	print *, "AIC=", -2.0d0*llk*(nins-1)+2.0d0*7
	print *, "BIC=", -2.0d0*llk*(nins-1)+log(dble(nins-1))*7
        call coverage_tpar(para_tpar,r,nins,obs)

pause	
	para_tpar(1:3)=para_par
	para_tpar(4:6)=para_tpar(1:3)
	para_tpar(6)=0.0d0
	call etparc(para_tpar,r,obs(1:nins),llk,'U')
        call dc_tpar(para_tpar,r,obs,asycov_tpar)
	print *, "===Summary of results: TPoiARc==="
	print *, "parameter estimates:"
	write( *,'(A,I0)') "threshold=", r
	write(*,'(A)') "Estimated parameters and s.d."
	write(*,'(2F10.4)') (para_tpar(i),sqrt(asycov_tpar(i,i)),i=1,6)
	print *, "LLK=",llk
	print *, "AIC=", -2.0d0*llk*(nins-1)+2.0d0*7
	print *, "BIC=", -2.0d0*llk*(nins-1)+log(dble(nins-1))*7
        call coverage_tpar(para_tpar,r,nins,obs)

	print *,"END"
pause	

	para_trendpar=0.0d0
	para_trendpar(1:3)=para_par
	call etrendpar(para_trendpar,prd,obs(1:nins),llk)
        call dc_trendpar(para_trendpar,prd,obs(1:nins),asycov_trendpar)
	print *, "Summary of results: PoiAR with trend:"
	write(*,'(A)') "Estimated parameters and s.d."
	write(*,'(2F10.4)') (para_trendpar(i),sqrt(asycov_trendpar(i,i)),i=1,7)
	!print *, "averaged LLK=",llk
	!print *, "AIC=", -2.0d0*llk*(nins-1)+2.0d0*7
	!print *, "BIC=", -2.0d0*llk*(nins-1)+log(dble(nins-1))*7
        call coverage_trendpar(para_trendpar,prd,nins,obs)


	para_ttpar=0.0d0
	para_ttpar(1:6)=para_tpar
        call ettpar(para_ttpar,r,obs(1:nins),llk,'U')
	call dc_ttpar(para_ttpar,r,obs(1:nins),asycov)
        write(*,*) "Estimated threshold:",r
        write(*,*) "Estimated LLK:",llk
        write(*,*) "Estimated AIC:", -2.0d0*llk*(nins-1)+2.0d0*11
        write(*,*) "Estimated BIC:", -2.0d0*llk*(nins-1)+log(dble(nins-1))*11
        write(*,*) "Estimated paramaters and standard error:"
        write(*,'(2F10.3)') (para_ttpar(i),sqrt(asycov(i,i)/size(obs)),i=1,10)
        call coverage_ttpar(para_ttpar,r,nins,obs)
!!
	deallocate(obs)
        end program main_program
