!include "nvst.f90"

subroutine sim_tpar(r,para,obs)
        implicit none
        integer :: r
        real*8 :: para(6)
        integer :: obs(:)

        INTEGER :: MSEED, MSTATE, M
        PARAMETER (MSEED=1,MSTATE=633,M=1)
        INTEGER :: GENID, I, IFAIL, LSEED, LSTATE, SUBID
        INTEGER :: SEED(MSEED), STATE(MSTATE), X(1)
        
        integer :: n
        integer :: y(int(size(obs)*1.5))

        real*8 :: lambda(M)
        
        n=size(y)

        GENID = 3
        SUBID = 1
        LSTATE = MSTATE
        LSEED = MSEED
        IFAIL = 0
        call g05kgf(genid,subid,state,lstate,ifail)
        IF (IFAIL.NE.0) then 
                print *, "cannot initialize random number generator."
                pause
        end if

        y(1)=1
        lambda=1.0d0
        do i=2,n
                if(y(i-1) .le. r) then
                        lambda=para(1)+para(2)*lambda+para(3)*y(i-1)
                else
                        lambda=para(4)+para(5)*lambda+para(6)*y(i-1)
                end if
                CALL G05TKF(M,LAMBDA,STATE,X,IFAIL)
                y(i)=x(1)
        end do
        obs=y(n-size(obs)+1:n)
        end subroutine sim_tpar

	  program main
	  !use nvst_int
        integer,parameter :: repl=1,nobs=1000,mlag=5
        integer :: r,r0,j
        real*8 :: para(6),para0(6),epara(repl,6),er(repl),ac(mlag)
        integer :: obs(nobs)

        real*8 :: nllk,llk,objgrd(6)
        integer :: nstate,iuser(2),mode,i


        !r=7
        !para(1:3)=(/0.5d0,0.8d0,0.4d0/)
        !para(4:6)=(/0.2d0,0.7d0,0.1d0/)
        r=6
        para(1:3)=(/0.5d0,0.8d0,0.7d0/)
        para(4:6)=(/0.2d0,0.2d0,0.1d0/)
        r0=r
        para0=para

        do i=1,repl
            call sim_tpar(r0,para0,obs)
		do j=0,nobs
		    write( *, '(5F8.3)') obs(j)
		end do
	

		!call acf(dble(obs),mlag,ac)
		!write( *, '(5F8.3)') ac
pause
        end do
        end program 



