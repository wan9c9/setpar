        module tpar_mod
        use mat_ope_nag
        use nvst_int
        use gnuplot_int
        implicit none

        contains

!	function outer_product(x,y)
!	implicit none
!	real*8 :: x(:),y(:)
!	real*8 :: outer_product(size(x),size(y))
!
!	integer :: i
!
!	do i=1,size(x)
!		outer_product(i,:)=x(i)*y
!	end do
!	end function outer_product

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
        IFAIL = 1
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


        subroutine nllk_tpar(r,para,obs,nllk)
        !nllk=negative log-likilihood
        implicit none
        integer,intent(in) :: r
        real*8,intent(in) :: para(6)
        integer,intent(in) :: obs(:)
        real*8,intent(out) :: nllk

        integer :: i
        real*8 :: lambda

	if(obs(1).gt.0.0d0) then
		lambda=obs(1)
	else
		lambda=0.01d0
	end if

        nllk=0.0d0
        do i=2,size(obs)
                if(obs(i-1).le.r) then
                        lambda=para(1)+para(2)*lambda+para(3)*obs(i-1)
                else
                        lambda=para(4)+para(5)*lambda+para(6)*obs(i-1)
                end if
                nllk=nllk+lambda-obs(i)*log(lambda)
        end do
        nllk=nllk/(size(obs)-1)
        end subroutine nllk_tpar

        subroutine grd_tpar(r,para,obs,grd)
        implicit none
        integer,intent(in) :: r
        real*8,intent(in) :: para(6)
        integer,intent(in) :: obs(:)
        real*8,intent(out) :: grd(6)

        integer :: i
        real*8 :: lambda,jv(3),dl1(3),dl2(3)

	if(obs(1).gt.0.0d0) then
        	lambda=obs(1)
	else
		lambda=0.01d0
	end if

        grd=0.0d0
        dl1=0.0d0
        dl2=0.0d0
        do i=2,size(obs)
                jv=(/1.0d0,lambda,dble(obs(i-1))/)
                if(obs(i-1).le.r) then
                        dl1=jv+para(2)*dl1
                        dl2=para(2)*dl2
                        lambda=para(1)+para(2)*lambda+para(3)*obs(i-1)
                else
                        dl1=para(5)*dl1
                        dl2=jv+para(5)*dl2
                        lambda=para(4)+para(5)*lambda+para(6)*obs(i-1)
                end if
                grd(1:3)=grd(1:3)+(1.0d0-obs(i)/lambda)*dl1
                grd(4:6)=grd(4:6)+(1.0d0-obs(i)/lambda)*dl2
        end do
        grd=grd/(size(obs)-1)
        !print *,"GRADIENT=", grd
        end subroutine grd_tpar

        subroutine objfun_tpar(MODE,N,X,OBJF,OBJGRD,NSTATE,IUSER,RUSER)
        implicit none

        integer,intent(inout) :: mode
        integer,intent(in) :: n
        real*8, intent(in) :: x(n)
        real*8,intent(out) :: objf
        real*8,intent(out) :: objgrd(n)
        integer,intent(in) :: nstate
        integer,intent(in) :: iuser(2)
        real*8,intent(in) :: ruser(iuser(2))

        mode=2
        call nllk_tpar(iuser(1),x,int(ruser),objf)
        call grd_tpar(iuser(1),x,int(ruser),objgrd)
        end subroutine objfun_tpar

             
        subroutine etpar(para,r,obs,llk,th)  
        implicit none
        real*8 :: para(6)
        integer :: r
        integer :: obs(:)
        real*8 :: llk
        character*1 :: th

        real*8 :: a(6),bl(7),bu(7),c(1),cjac(1,1),clamda(7),&
                  objf,objgrd(6),work(203),rmtrx(6,6)
        integer :: iter,istate(7),iwork(19),ifail,iuser(2),i
       
        real*8 :: alpha(2),qv(2)
        integer :: q(2)

        real*8, allocatable :: parath(:,:),llkth(:)
        external e04udm

        alpha=(/0.1d0,0.9d0/)
        a=0.0d0
        a(5:6)=1.0d0
        bl(1:7)=0.001d0
        bu(1)=1.0d4
        bu(2)=1.0d1
        bu(3)=1.0d4
        bu(4)=1.0d4
        bu(5:6)=1.0d0
        bu(7)=1.0d0
        
        if(th .eq. 'K') then
                iuser(1)=r
                iuser(2)=size(obs)
                !call e04uef('Verify Level = 3')
                call e04uef('Print Level = 0')
                call e04uef('Major Iteration Limit = 2000')
		ifail=0
                call e04ucf(6,1,0,1,1,6,a,bl,bu,e04udm,&
                   objfun_tpar,iter,istate,c,cjac,clamda,objf,objgrd,rmtrx,para,&
                   iwork,19,work,203,iuser,dble(obs),ifail)
                llk=-objf
        else
                call g01amf(size(obs),dble(obs),2,alpha,qv,ifail)
                if(ifail.ne.0) then 
                        print *,"Error in finding empirical quantiles &
                                & of observations."
                        pause
                end if
                q=int(qv)
                write(*,*) "Since threshold is unknown, it will &
                             & be searched."
                write(*,*) "From ",q(1), " to ", q(2)
                allocate(parath(q(1):q(2),6))
                allocate(llkth(q(1):q(2)))
                write(*,*) "Threshold  Log-likilihood & 
                              &  Current_best_threshold" 
                do i=q(1),q(2)
                        parath(i,:)=para
                        iuser(1)=i
                        iuser(2)=size(obs)
                        !call e04uef('Verify Level = 3')
                        call e04uef('Print Level = 0')
                        call e04uef('Major Iteration Limit = 2000')
			ifail=0
                        call e04ucf(6,1,0,1,1,6,a,bl,bu,e04udm,&
                           objfun_tpar,iter,istate,c,cjac,clamda,objf,objgrd,rmtrx,parath(i,:),&
                           iwork,19,work,203,iuser,dble(obs),ifail)
                        llkth(i)=-objf
                        if( i.eq.q(1)) then
                                r=i
                        else
                                if( llkth(i).gt.llkth(r)) r=i
                        end if
                        write(*,'(4X,I, 5X,F12.5,2X, I )') ,i,llkth(i),r
                 end do
                 para=parath(r,:)
                 llk=llkth(r)
                 deallocate(parath,llkth)
        end if

        end subroutine etpar

        subroutine dc_tpar(para,r,obs,asycov)
        implicit none
        real*8 :: para(6)
        integer :: r
        integer :: obs(:)
       	real*8 :: asycov(6,6) 

        integer :: i
        real*8 :: lambda(size(obs)),rsdl(size(obs)), mrsdl,sdrsdl
	real*8 :: dlambda(size(obs),6),deri(6)
       	real*8 :: asy(6,6) 

	if(obs(1).gt.0.0d0) then
        	lambda(1)=obs(1)
	else
		lambda(1)=0.01d0
	end if
        rsdl=0.0d0
	dlambda=0.0d0
	asy=0.0d0
        do i=2,size(obs)
                if(obs(i-1).le.r) then
                        lambda(i)=para(1)+para(2)*lambda(i-1)+para(3)*obs(i-1)
			deri(1:3)=(/1.0d0,lambda(i-1),dble(obs(i-1))/)
			deri(4:6)=0.0d0
			dlambda(i,:)=deri+para(2)*dlambda(i-1,:)
                else
                        lambda(i)=para(4)+para(5)*lambda(i-1)+para(6)*obs(i-1)
			deri(1:3)=0.0d0
			deri(4:6)=(/1.0d0,lambda(i-1),dble(obs(i-1))/)
			dlambda(i,:)=deri+para(5)*dlambda(i-1,:)
                end if
		asy=asy+1.0d0/lambda(i)*outer_product(dlambda(i,:),dlambda(i,:))
                rsdl(i)=(obs(i)-lambda(i))/sqrt(lambda(i))
        end do
	asy=asy
        call sub_inv_pd_mat(asy, asycov,i)
	if(i.ne.0) print *,"Error in finding inv. of asycov" 
	!write(*, '(A,F10.3)') "Mean of (Y_t-lambda_t)^2", sum((obs-lambda)**2.0d0)/size(obs)
	!write(*, '(A,F10.3)') "Mean of (Y_t-lambda_t)^2/lambda_t", sum((obs-lambda)**2.0d0/lambda)/size(obs)
	!open(unit=11,file="fitting_summary.dat")
	!write(11,'(I5,I5,1X,F10.3,F10.3)') (i,obs(i),lambda(i),rsdl(i),i=1,size(obs))
	!close(11)
	!write(*,'(A)')  "Statistics of the Pearson residual:"
	!mrsdl=sum(rsdl)/size(obs)
	!sdrsdl=sqrt(sum((rsdl-mrsdl)**2.0d0)/size(obs))
	!write(*, '(A,F10.3)') "mean", mrsdl
	!write(*, '(A,F10.3)') "sd", sdrsdl
	!write(*, '(A,F10.3)') "skewness", sum((rsdl-mrsdl)**3.0d0)/size(obs)/sdrsdl**3.0d0
	!write(*, '(A,F10.3)') "ex. kur.", sum((rsdl-mrsdl)**4.0d0)/size(obs)/sdrsdl**4.0d0-3.0d0
        end subroutine dc_tpar



        subroutine simulation_threshold_tpar
        implicit none

        integer,parameter :: repl=1000,nobs=500,mlag=5
        integer :: r,r0
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
		call acf(dble(obs),mlag,ac)
		write( *, '(5F8.3)') ac
pause
                para=para0
		r=r0
                call etpar(para,r,obs,llk,'U')
                !pause  
                epara(i,:)=para
                er(i)=r
        end do
        print *,"sample size=",nobs
        print *,"replication=", repl
        print *,"True parameters"
        write(*,*) "Threshold:", r0
        write(*,*)"Other para:"
        write(*, '(3F12.3)') para0
        print *, "Average of estimated parameter:"
        write(*, '(3F12.3)') sum(epara,dim=1)/repl
        print *, "Covariance matrix:"
        call print_mat(covm(epara))

        write(*,*) "Estimated threshold:" 
        write(*,*) "Mean", mean(dble(er))
        write(*,*) "Variance", var(dble(er))
        !open(unit=11,file="er.dat",status="replace")
        !do i=1,repl
        !        write(11,*) er(i)
        !end do
        !close(11)
        !call gplot1(dble(er))
        end subroutine simulation_threshold_tpar


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        subroutine test_asycov_tpar
        implicit none
        integer,parameter :: repl=1000,nobs=3000
        integer :: r,r0
        real*8 :: para(6),para0(6),epara(repl,6),er(repl)
        integer :: obs(nobs)

        real*8 :: nllk,llk,objgrd(6),asycov(repl,6,6)
        integer :: nstate,iuser(2),mode,i


        !r=7
        !para(1:3)=(/0.5d0,0.7d0,0.2d0/)
        !para(4:6)=(/0.3d0,0.4d0,0.5d0/)
        
	r=6
        para(1:3)=(/0.5d0,0.8d0,0.7d0/)
        para(4:6)=(/0.2d0,0.2d0,0.1d0/)

	!r=6
	!para(1:3)=(/0.3d0,0.8d0,0.1d0/)
        !para(4:6)=(/0.1d0,0.7d0,0.25d0/)
        
	r0=r
        para0=para

        do i=1,repl
		!print *,"index repl=", i
                call sim_tpar(r0,para0,obs)
                para=para0
                !call etpar(para,r,obs,llk,'K')
                call etpar(para,r,obs,llk,'U')
		call dc_tpar(para,r,obs,asycov(i,:,:))
                epara(i,:)=para
                er(i)=r
        end do
        print *,"sample size=",nobs
        print *,"replication=", repl
        print *,"True parameters"
        write(*,*) "Threshold:", r0
        write(*,*) "Other para:"
        write(*, '(3F12.3)') para0
        print *, "Average of estimated parameter:"
        write(*, '(3F12.3)') sum(epara,dim=1)/repl
        print *, "Covariance matrix:"
	print *, "N*(Sample Variance):"
	!write(*, '(6F8.2)') (nobs*(covm(epara))(i,i),i=1,6)
        call print_mat(nobs*covm(epara))
	
	print *, "N*(Sample mean of estimated variance):"
	!write(*, '(6F8.2)') (nobs*(covm(epara))(i,i),i=1,6)
	call print_mat(nobs*sum(asycov,dim=1)/repl)
        
	write(*,*) "Estimated threshold:" 
        write(*,*) "Mean", mean(dble(er))
        write(*,*) "Variance", var(dble(er))
        open(unit=11,file="er.dat",status="replace")
        do i=1,repl
                write(11,'(10F10.3)') er(i), epara(i,:)
        end do
        close(11)
        !call gplot1(dble(er))
        end subroutine test_asycov_tpar


        subroutine simulation_tpar
        implicit none

        integer,parameter :: repl=1000,nobs=1000
        integer :: r
        real*8 :: para(6),para0(6),epara(repl,6)
        integer :: obs(nobs)

        real*8 :: nllk,llk,objgrd(6)
        integer :: nstate,iuser(2),mode,i


        !r=4
        !para(1:3)=(/0.3d0,0.9d0,0.4d0/)
        !para(4:6)=(/1.0d0,0.5d0,0.1d0/)
        r=8
        para(1:3)=(/0.30d0,0.88d0,0.05d0/)
        para(4:6)=(/0.14d0,0.73d0,0.23d0/)


        para0=para
        do i=1,repl
                call sim_tpar(r,para0,obs)
                para=para0
                call etpar(para,r,obs,llk,'K')
                !pause  
                epara(i,:)=para
        end do
        print *,"sample size=",nobs
        print *,"replication=", repl
        print *,"True parameters"
        write(*, '(30E12.3)') para0
        print *, "Average of estimated parameter:"
        write(*, '(30E12.3)') sum(epara,dim=1)/repl
        print *, "Covariance matrix:"

        call print_mat(covm(epara))
        end subroutine simulation_tpar
        
        subroutine coverage_tpar(para,r,nins,obs)
	implicit none
	integer,intent(in) :: r
        real*8,intent(in) :: para(6)
	integer :: nins
        integer,intent(in) :: obs(:)

        integer :: i
        real*8 :: lambda(size(obs))
	
	if(obs(1).gt.0.0d0) then
        	lambda(1)=obs(1)
	else
		lambda(1)=0.01d0
	end if

        do i=2,size(obs)
                if(obs(i-1).le.r) then
                        lambda(i)=para(1)+para(2)*lambda(i-1)+para(3)*obs(i-1)
                else
                        lambda(i)=para(4)+para(5)*lambda(i-1)+para(6)*obs(i-1)
                end if
        end do
	open(unit=11,file="tpar_data_and_fitted.dat")
	write(11,'(I5,1X,F10.3)') (obs(i),lambda(i),i=1,size(obs))
	close(11)
	write(*,'(A)') "===OUTPUT OF COVERAGE_TPAR==="
	write(*,'(A)') "For the TPoiAR model:"
	write(*,'(A,I0)') "IN SAMPLE NO. OF OBS:", nins
	write(*,'(A,F15.5)') "in-sample mse:",sum((obs(1:nins)-lambda(1:nins))**2.0d0)/dble(nins)
	write(*,'(A,I0)') "OUT-OF-SAMPLE SAMPLE NO. OF OBS:", size(obs)-nins
	write(*,'(A,F15.5)') "out-of-sample mse:",sum((obs(nins+1:)-lambda(nins+1:))**2.0d0)/dble(size(obs)-nins)
	write(*,'(A)') "END===OUTPUT OF COVERAGE_TPAR===END"
        end subroutine coverage_tpar
        
       

	subroutine rd_tpar
        implicit none

        integer,allocatable :: obs(:)
        integer :: r,ierr
        real*8 :: para(6)
        character(len=20) :: fn

	integer :: i,nins,ntot
        real*8 :: nllk,llk,objgrd(6),asycov(6,6),parapar(3)


        !write(*,*) "Please type name of data file."
        !read(*,*) fn
	!fn="ntrans5.dat"

	!fn="ntrans5t.dat"
	!nins=920
	!ntot=1380
	!parapar=(/2.4958937d0,0.5276162d0,0.4038537d0/)
	
	!fn="sunspots.dat"
	!nins=440
	!ntot=459
	!prd=11

	!fn="blowfly.dat"
	!nins=340
	!ntot=361

	!fn="poliom.dat"
	!nins=160
	!ntot=168

	!fn="firearm_homocide.dat"
	!nins=300
	!ntot=313


	!fn="coli.txt"
	!nins=120
	!ntot=125
	!parapar=(/0.71d0,0.22d0,0.41d0/)
	
        open(unit=10,file=fn,status="old",iostat=ierr)
        if(ierr.ne.0) then
                print *, "Filename error."
                return
        end if
        allocate(obs(ntot))
        read(10,*) obs
        close(10)

        para(1:3)=(/0.71d0,0.22d0,0.41d0/)
        para(4:6)=(/0.71d0,0.22d0,0.41d0/)
        call etpar(para,r,obs(1:nins),llk,'U')
	print *, "Summary of results: TPoiAR:"
	print *, "LLK=",llk
	print *, "AIC=", -2.0d0*llk*(nins-1)+2.0d0*7
	call dc_tpar(para,r,obs(1:nins),asycov)
        write(*,*) "Estimated threshold:",r
        write(*,*) "Estimated paramaters and standard error:"
        write(*,'(2F10.3)') (para(i),sqrt(asycov(i,i)/size(obs)),i=1,6)
        
	call coverage_tpar(para,r,nins,obs)

	parapar=(/0.22d0,0.726d0,0.198d0/)
	!call coverage_par(parapar,nins,obs)
        deallocate(obs)
        end subroutine rd_tpar

        end module tpar_mod
