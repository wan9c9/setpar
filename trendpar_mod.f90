        module trendpar_mod
        use mat_ope_nag
        use nvst_int
        use gnuplot_int
        implicit none

contains
        subroutine nllk_trendpar(para,prd,obs,nllk)
        !nllk=negative log-likilihood
        implicit none
        real*8,intent(in) :: para(7)
        integer,intent(in) :: prd
        integer,intent(in) :: obs(:)
        real*8,intent(out) :: nllk

        integer :: i
        real*8 :: lambda

        nllk=0.0d0
        lambda=obs(1)
	if(obs(1).eq.0.0d0) lambda=0.01d0

        do i=2,size(obs)
                lambda=para(1)+para(2)*lambda+para(3)*obs(i-1)
		lambda=lambda+para(4)*cos(2.0d0*pi()/prd*i)+para(5)*sin(2.0d0*pi()/prd*i)&
			+para(6)*cos(2.0d0*pi()/prd*2.0d0*i)+para(7)*sin(2.0d0*pi()/prd*2.0d0*i)
                nllk=nllk+lambda-obs(i)*log(lambda)
        end do
        nllk=nllk/(size(obs)-1)
        end subroutine nllk_trendpar

        subroutine grd_trendpar(para,prd,obs,grd)
        implicit none
        real*8,intent(in) :: para(7)
	integer,intent(in) :: prd
        integer,intent(in) :: obs(:)
        real*8,intent(out) :: grd(7)

        integer :: i
        real*8 :: lambda,jv(3),dl1(3),dl3(4),jvp(4)

        lambda=obs(1)
	if(obs(1).eq.0.0d0) lambda=0.01d0

        dl1=0.0d0
	dl3=0.0d0
        grd=0.0d0
        do i=2,size(obs)
                jv=(/1.0d0,lambda,dble(obs(i-1))/)
		jvp=(/cos(2.0d0*pi()/prd*i),sin(2.0d0*pi()/prd*i),&
                   cos(2.0d0*pi()/prd*2.0d0*i),sin(2.0d0*pi()/prd*2.0d0*i)/)
		dl1=jv+para(2)*dl1
		dl3=jvp+para(2)*dl3
		lambda=para(1)+para(2)*lambda+para(3)*obs(i-1)+dot_product(para(4:7),jvp)
                grd(1:3)=grd(1:3)+(1.0d0-obs(i)/lambda)*dl1
		grd(4:7)=grd(4:7)+(1.0d0-obs(i)/lambda)*dl3
        end do
        grd=grd/(size(obs)-1)
        !print *,"GRADIENT=", grd
        end subroutine grd_trendpar

        !subroutine objfun_ttpar(MODE,N,X,OBJF,OBJGRD,NSTATE,IUSER,RUSER)
        !implicit none
!
!        integer,intent(inout) :: mode
!        integer,intent(in) :: n
!        real*8, intent(in) :: x(n)
!        real*8,intent(out) :: objf
!        real*8,intent(out) :: objgrd(n)
!        integer,intent(in) :: nstate
!        integer,intent(in) :: iuser(1)
!        real*8,intent(in) :: ruser(iuser(1))
!
!        mode=2
!        call nllk_trendpar(x,int(ruser),objf)
!        call grd_trendpar(x,int(ruser),objgrd)
!        end subroutine objfun_ttpar

	subroutine objfun_trendpar(N,X,OBJF,OBJGRD,IUSER,RUSER)
        implicit none
        integer,intent(in) :: n
        real*8, intent(in) :: x(n)
        real*8,intent(out) :: objf
        real*8,intent(out) :: objgrd(n)
        integer,intent(in) :: iuser(2)
        real*8,intent(in) :: ruser(iuser(2))

        call nllk_trendpar(x,iuser(1),int(ruser),objf)
        call grd_trendpar(x,iuser(1),int(ruser),objgrd)
        end subroutine objfun_trendpar


             
        subroutine dc_trendpar(para,prd,obs,asycov)
        implicit none
        real*8 :: para(7)
        integer :: prd
        integer :: obs(:)
       	real*8 :: asycov(7,7) 

        integer :: i
        real*8 :: lambda(size(obs)),rsdl(size(obs))
	real*8 :: dlambda(size(obs),7),deri(7),jvp(4)
      	real*8 :: asy(7,7) 
	real*8 :: mrsdl, sdrsdl

	lambda(1)=obs(1)
	if(obs(1).eq.0.0d0) lambda(1)=0.01d0
	rsdl=0.0d0
	dlambda=0.0d0
	asy=0.0d0
	do i=2,size(obs)
		jvp=(/cos(2.0d0*pi()/prd*i),sin(2.0d0*pi()/prd*i),&
		    cos(2.0d0*pi()/prd*2.0d0*i),sin(2.0d0*pi()/prd*2.0d0*i)/)
		lambda(i)=para(1)+para(2)*lambda(i-1)+para(3)*obs(i-1)+dot_product(para(4:7),jvp)
		deri(1:3)=(/1.0d0,lambda(i-1),dble(obs(i-1))/)
		deri(4:7)=jvp
		dlambda(i,:)=deri+para(2)*dlambda(i-1,:)
		asy=asy+1.0d0/lambda(i)*outer_product(dlambda(i,:),dlambda(i,:))
		rsdl(i)=(obs(i)-lambda(i))/sqrt(lambda(i))
	end do
	asy=asy
        call sub_inv_pd_mat(asy, asycov,i)
	if(i.ne.0) print *,"Error in finding inv. of asycov" 
!	open(unit=11,file="ttpar_summary.dat")
!	write(11,'(I5,1X,F10.3)') (obs(i),lambda(i),i=1,size(obs))
!	close(11)
!	write(*,'(A)')  "Statistics of the Pearson residual:"
!	mrsdl=sum(rsdl)/size(obs)
!	sdrsdl=sqrt(sum((rsdl-mrsdl)**2.0d0)/size(obs))
!	write(*, '(A,F10.3)') "mean", mrsdl
!	write(*, '(A,F10.3)') "sd", sdrsdl
!	write(*, '(A,F10.3)') "skewness", sum((rsdl-mrsdl)**3.0d0)/size(obs)/sdrsdl**3.0d0
!	write(*, '(A,F10.3)') "ex. kur.", sum((rsdl-mrsdl)**4.0d0)/size(obs)/sdrsdl**4.0d0-3.0d0
        end subroutine dc_trendpar

        subroutine coverage_trendpar(para,prd,nins,obs)
        implicit none
        real*8,intent(in) :: para(7)
        integer,intent(in) :: prd,nins
        integer,intent(in) :: obs(:)

        integer :: i
        real*8 :: lambda(size(obs)),rsdl(size(obs)),nllk
!	real*8 :: dlambda(size(obs),10),deri(10),jvp(4)
!       	real*8 :: asy(10,10) 
!	real*8 :: mrsdl, sdrsdl
!
	lambda(1)=obs(1)
	if(obs(1).eq.0.0d0) lambda(1)=0.01d0

	do i=2,size(obs)
                lambda(i)=para(1)+para(2)*lambda(i-1)+para(3)*obs(i-1)
		lambda(i)=lambda(i)+para(4)*cos(2.0d0*pi()/prd*i)+para(5)*sin(2.0d0*pi()/prd*i)&
			+para(6)*cos(2.0d0*pi()/prd*2.0d0*i)+para(7)*sin(2.0d0*pi()/prd*2.0d0*i)
        end do
	nllk=sum(lambda(2:nins)-obs(2:nins)*log(lambda(2:nins)))
!	!write(*, '(A,F10.3)') "Mean of (Y_t-lambda_t)^2:", sum((obs-lambda)**2.0d0)/size(obs)
!	!write(*, '(A,F10.3)') "Mean of (Y_t-lambda_t)^2/lambda_t:", sum((obs-lambda)**2.0d0/lambda)/size(obs)
!	open(unit=11,file="ttpar_data_and_fitted.dat")
!	write(11,'(I5,1X,F10.3)') (obs(i),lambda(i),i=1,size(obs))
!	close(11)
	write(*,'(A)') "===OUTPUT OF COVERAGE_PAR WITH TREND==="
	write(*,'(A)') "For PoiAR with Fourier trend:"
	write(*,'(A,F15.5)') "LLK=", -nllk
	write(*,'(A,F15.5)') "AIC=", 2.0d0*nllk+2.0d0*7
	write(*,'(A,F15.5)') "BIC=", 2.0d0*nllk+log(dble(nins-1))*7
	write(*,'(A,F15.5)') "in-sample mse:",sum((obs(1:nins)-lambda(1:nins))**2.0d0)/dble(nins)
	write(*,'(A,F15.5)') "out-of-sample mse:",sum((obs(nins+1:)-lambda(nins+1:))**2.0d0)/dble(size(obs)-nins)
	write(*,'(A)') "END===OUTPUT OF COVERAGE_PAR WITH TREND===END"
        end subroutine coverage_trendpar

	subroutine etrendpar(p,prd,obs,llk)
	implicit none
	real*8,intent(inout) :: p(7)
	integer,intent(in) :: prd
	integer,intent(in) :: obs(:)
	real*8, intent(out) :: llk

	real*8 :: bl(7),bu(7),nllk,gc(7),w(98)
	integer :: ia(1),iw(9),ifail
	
	bl=1.0d-5
	bl(4:7)=-1.0d2
	bu=1.0d0
	bu(1)=1.0d2
	bu(4:7)=1.0d2

	ifail=-1
	ia(1)=prd
	ia(2)=size(obs)
	call e04kzf(7,0,objfun_trendpar,bl,bu,&
	            p,nllk,gc,iw,9,&
		    w,98,ia,dble(obs),ifail)
	llk=-nllk
	print *, "ESTIMATED PARAMETERS USING ANALYTIC GRADIENTS:"
	print *, p
	write(*,'(A,F12.5)') "averaged LLK=",llk
	!if(ifail.ne.0) pause
	end subroutine etrendpar
!!


!        subroutine rd_ttpar
!        implicit none
!
!        integer:: nobs
!        
!        integer,allocatable :: obs(:)
!        integer :: r,ierr
!        real*8 :: para(10)
!        character(len=20) :: fn
!
!	integer :: i,nins,ntot
!        real*8 :: nllk,llk,objgrd(10),asycov(10,10)
!
!
!	!fn="ntrans5.dat"
!	!fn="coli.txt"
!	!nins=120
!	!ntot=125
	
!!	!fn="ntrans5t.dat"
!	!nins=920
!	!ntot=1380
!	!prd=92
!
!	!!sunspot data
!	fn="sunspots.dat"
!	nins=440
!	ntot=459
!	prd=11
!
!        open(unit=10,file=fn,status="old",iostat=ierr)
!        if(ierr.ne.0) then
!                print *, "Filename error."
!                return
!        end if
!        allocate(obs(ntot))
!        read(10,*) obs
!        close(10)
!
!	asycov=0.0d0
!	para=0.0d0
!        para(1:3)=(/0.5d0,0.4d0,0.2d0/)
!        para(4:6)=(/0.5d0,0.4d0,0.2d0/)
!
!        call ettpar(para,r,obs(1:nins),llk,'U')
!	call dc_ttpar(para,r,obs(1:nins),asycov)
!        write(*,*) "Estimated threshold:",r
!        write(*,*) "Estimated LLK:",llk
!        write(*,*) "Estimated AIC:", -2.0d0*llk*(nins-1)+2.0d0*11
!        write(*,*) "Estimated paramaters and standard error:"
!        write(*,'(2F10.3)') (para(i),sqrt(asycov(i,i)/size(obs)),i=1,10)
!        call coverage_ttpar(para,r,nins,obs)
!
!        end subroutine rd_ttpar

        end module trendpar_mod

