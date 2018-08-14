        module ttpar_mod
        use mat_ope_nag
        use nvst_int
        use gnuplot_int
        implicit none

	integer :: prd
        
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

        subroutine nllk_ttpar(r,para,obs,nllk)
        !nllk=negative log-likilihood
        implicit none
        integer,intent(in) :: r
        real*8,intent(in) :: para(10)
        integer,intent(in) :: obs(:)
        real*8,intent(out) :: nllk

        integer :: i
        real*8 :: lambda

        nllk=0.0d0
        lambda=obs(1)
	if(obs(1).eq.0.0d0) lambda=0.01d0

        do i=2,size(obs)
                if(obs(i-1).le.r) then
                        lambda=para(1)+para(2)*lambda+para(3)*obs(i-1)
                else
                        lambda=para(4)+para(5)*lambda+para(6)*obs(i-1)
                end if
		lambda=lambda+para(7)*cos(2.0d0*pi()/prd*i)+para(8)*sin(2.0d0*pi()/prd*i)&
			+para(9)*cos(2.0d0*pi()/prd*2.0d0*i)+para(10)*sin(2.0d0*pi()/prd*2.0d0*i)
                nllk=nllk+lambda-obs(i)*log(lambda)
        end do
        nllk=nllk/(size(obs)-1)
        end subroutine nllk_ttpar

        subroutine grd_ttpar(r,para,obs,grd)
        implicit none
        integer,intent(in) :: r
        real*8,intent(in) :: para(10)
        integer,intent(in) :: obs(:)
        real*8,intent(out) :: grd(10)

        integer :: i
        real*8 :: lambda,jv(3),dl1(3),dl2(3),dl3(4),jvp(4)

        grd=0.0d0
        lambda=obs(1)
	if(obs(1).eq.0.0d0) lambda=0.01d0
        dl1=0.0d0
        dl2=0.0d0
	dl3=0.0d0
        do i=2,size(obs)
                jv=(/1.0d0,lambda,dble(obs(i-1))/)
		jvp=(/cos(2.0d0*pi()/prd*i),sin(2.0d0*pi()/prd*i),&
                   cos(2.0d0*pi()/prd*2.0d0*i),sin(2.0d0*pi()/prd*2.0d0*i)/)
                if(obs(i-1).le.r) then
                        dl1=jv+para(2)*dl1
                        dl2=para(2)*dl2
			dl3=jvp+para(2)*dl3
                        lambda=para(1)+para(2)*lambda+para(3)*obs(i-1)+dot_product(para(7:10),jvp)
                else
                        dl1=para(5)*dl1
                        dl2=jv+para(5)*dl2
			dl3=jvp+para(5)*dl3
                        lambda=para(4)+para(5)*lambda+para(6)*obs(i-1)+dot_product(para(7:10),jvp)
                end if
                grd(1:3)=grd(1:3)+(1.0d0-obs(i)/lambda)*dl1
                grd(4:6)=grd(4:6)+(1.0d0-obs(i)/lambda)*dl2
		grd(7:10)=grd(7:10)+(1.0d0-obs(i)/lambda)*dl3
        end do
        grd=grd/(size(obs)-1)
        !print *,"GRADIENT=", grd
        end subroutine grd_ttpar

        subroutine objfun_ttpar(MODE,N,X,OBJF,OBJGRD,NSTATE,IUSER,RUSER)
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
        call nllk_ttpar(iuser(1),x,int(ruser),objf)
        call grd_ttpar(iuser(1),x,int(ruser),objgrd)
        end subroutine objfun_ttpar

             
        subroutine dc_ttpar(para,r,obs,asycov)
        implicit none
        real*8 :: para(10)
        integer :: r
        integer :: obs(:)
       	real*8 :: asycov(10,10) 

        integer :: i
        real*8 :: lambda(size(obs)),rsdl(size(obs))
	real*8 :: dlambda(size(obs),10),deri(10),jvp(4)
       	real*8 :: asy(10,10) 
	real*8 :: mrsdl, sdrsdl

        lambda(1)=obs(1)
	if(obs(1).eq.0.0d0) lambda(1)=0.01d0
        rsdl=0.0d0
	dlambda=0.0d0
	asy=0.0d0
	
        do i=2,size(obs)
		jvp=(/cos(2.0d0*pi()/prd*i),sin(2.0d0*pi()/prd*i),&
                    cos(2.0d0*pi()/prd*2.0d0*i),sin(2.0d0*pi()/prd*2.0d0*i)/)
                if(obs(i-1).le.r) then
                        lambda(i)=para(1)+para(2)*lambda(i-1)+para(3)*obs(i-1)+dot_product(para(7:10),jvp)
			deri(1:3)=(/1.0d0,lambda(i-1),dble(obs(i-1))/)
			deri(4:6)=0.0d0
			deri(7:10)=jvp
			dlambda(i,:)=deri+para(2)*dlambda(i-1,:)
                else
                        lambda(i)=para(4)+para(5)*lambda(i-1)+para(6)*obs(i-1)+dot_product(para(7:10),jvp)
			deri(1:3)=0.0d0
			deri(4:6)=(/1.0d0,lambda(i-1),dble(obs(i-1))/)
			deri(7:10)=jvp
			dlambda(i,:)=deri+para(5)*dlambda(i-1,:)
                end if
		asy=asy+1.0d0/lambda(i)*outer_product(dlambda(i,:),dlambda(i,:))
                rsdl(i)=(obs(i)-lambda(i))/sqrt(lambda(i))
        end do
	asy=asy/dble(size(obs)-1)
        call sub_inv_pd_mat(asy, asycov,i)
	write(*, '(A,F10.3)') "Mean of (Y_t-lambda_t)^2:", sum((obs-lambda)**2.0d0)/size(obs)
	write(*, '(A,F10.3)') "Mean of (Y_t-lambda_t)^2/lambda_t:", sum((obs-lambda)**2.0d0/lambda)/size(obs)
	if(i.ne.0) print *,"Error in finding inv. of asycov" 
	open(unit=11,file="ttpar_summary.dat")
	write(11,'(I5,1X,F10.3)') (obs(i),lambda(i),i=1,size(obs))
	close(11)
	write(*,'(A)')  "Statistics of the Pearson residual:"
	mrsdl=sum(rsdl)/size(obs)
	sdrsdl=sqrt(sum((rsdl-mrsdl)**2.0d0)/size(obs))
	write(*, '(A,F10.3)') "mean", mrsdl
	write(*, '(A,F10.3)') "sd", sdrsdl
	write(*, '(A,F10.3)') "skewness", sum((rsdl-mrsdl)**3.0d0)/size(obs)/sdrsdl**3.0d0
	write(*, '(A,F10.3)') "ex. kur.", sum((rsdl-mrsdl)**4.0d0)/size(obs)/sdrsdl**4.0d0-3.0d0
        end subroutine dc_ttpar

        subroutine coverage_ttpar(para,r,nins,obs)
        implicit none
        real*8 :: para(10)
        integer :: r,nins
        integer :: obs(:)

        integer :: i
        real*8 :: lambda(size(obs)),rsdl(size(obs))
	real*8 :: dlambda(size(obs),10),deri(10),jvp(4)
       	real*8 :: asy(10,10) 
	real*8 :: mrsdl, sdrsdl

	lambda(1)=obs(1)
	if(obs(1).eq.0.0d0) lambda=0.01d0
        do i=2,size(obs)
                if(obs(i-1).le.r) then
                        lambda(i)=para(1)+para(2)*lambda(i-1)+para(3)*obs(i-1)
                else
                        lambda(i)=para(4)+para(5)*lambda(i-1)+para(6)*obs(i-1)
                end if
		lambda(i)=lambda(i)+para(7)*cos(2.0d0*pi()/prd*i)+para(8)*sin(2.0d0*pi()/prd*i)&
			+para(9)*cos(2.0d0*pi()/prd*2.0d0*i)+para(10)*sin(2.0d0*pi()/prd*2.0d0*i)
		
        end do
	!write(*, '(A,F10.3)') "Mean of (Y_t-lambda_t)^2:", sum((obs-lambda)**2.0d0)/size(obs)
	!write(*, '(A,F10.3)') "Mean of (Y_t-lambda_t)^2/lambda_t:", sum((obs-lambda)**2.0d0/lambda)/size(obs)
	open(unit=11,file="ttpar_data_and_fitted.dat")
	write(11,'(I5,1X,F10.3)') (obs(i),lambda(i),i=1,size(obs))
	close(11)
	write(*,'(A)') "For the TPoiAR with Fourier trend:"
	write(*,'(A,D15.5)') "in-sample mse:",sum((obs(1:nins)-lambda(1:nins))**2.0d0)/dble(nins)
	write(*,'(A,D15.5)') "out-of-sample mse:",sum((obs(nins+1:)-lambda(nins+1:))**2.0d0)/dble(size(obs)-nins)
	!write(*,'(A)')  "Statistics of the Pearson residual:"
	!mrsdl=sum(rsdl)/size(obs)
	!sdrsdl=sqrt(sum((rsdl-mrsdl)**2.0d0)/size(obs))
	!write(*, '(A,F10.3)') "mean", mrsdl
	!write(*, '(A,F10.3)') "sd", sdrsdl
	!write(*, '(A,F10.3)') "skewness", sum((rsdl-mrsdl)**3.0d0)/size(obs)/sdrsdl**3.0d0
	!write(*, '(A,F10.3)') "ex. kur.", sum((rsdl-mrsdl)**4.0d0)/size(obs)/sdrsdl**4.0d0-3.0d0
        end subroutine coverage_ttpar



        subroutine ettpar(para,r,obs,llk,th)  
        implicit none
        real*8 :: para(10)
        integer :: r
        integer :: obs(:)
        real*8 :: llk
        character*1 :: th

        real*8 :: a(10),bl(11),bu(11),c(1),cjac(1,1),clamda(11),&
                  objf,objgrd(10),work(411),rmtrx(10,10)
        integer :: iter,istate(11),iwork(31),ifail,iuser(2),i
       
        real*8 :: alpha(2),qv(2)
        integer :: q(2)

        real*8, allocatable :: parath(:,:),llkth(:)
        external e04udm

        alpha=(/0.25d0,0.75d0/)
        a=0.0d0
        a(5:6)=1.0d0
        bl(1:6)=0.001d0
	bl(7:10)=-1.0d2
	bl(11)=1.0d-2
        bu(1)=1.0d2
        bu(2:3)=1.0d1
        bu(4)=1.0d2
        bu(5:6)=1.0d0
	bu(7:10)=1.0d2
        bu(11)=1.0d0
        
        if(th .eq. 'K') then
                iuser(1)=r
                iuser(2)=size(obs)
                !call e04uef('Verify Level = 3')
                call e04uef('Print Level = 0')
                call e04uef('Major Iteration Limit = 2000')
		ifail=0
                call e04ucf(10,1,0,1,1,10,a,bl,bu,e04udm,&
                   objfun_ttpar,iter,istate,c,cjac,clamda,objf,objgrd,rmtrx,para,&
                   iwork,31,work,411,iuser,dble(obs),ifail)
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
                allocate(parath(q(1):q(2),10))
                allocate(llkth(q(1):q(2)))
                write(*,*) "Threshold  Log-likilihood & 
                              &  Current_best_threshold" 
                do i=q(1),q(2)
			if(i.eq.q(1)) then
                        	parath(i,:)=para
			else
				parath(i,:)=parath(i-1,:)
			end if
                        iuser(1)=i
                        iuser(2)=size(obs)
                        !call e04uef('Verify Level = 3')
                        call e04uef('Print Level = 0')
                        call e04uef('Major Iteration Limit = 2000')
			ifail=0
                	call e04ucf(10,1,0,1,1,10,a,bl,bu,e04udm,&
                           objfun_ttpar,iter,istate,c,cjac,clamda,objf,objgrd,rmtrx,parath(i,:),&
                           iwork,31,work,411,iuser,dble(obs),ifail)
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

        end subroutine ettpar

        subroutine rd_ttpar
        implicit none

        integer:: nobs
        
        integer,allocatable :: obs(:)
        integer :: r,ierr
        real*8 :: para(10)
        character(len=20) :: fn

	integer :: i,nins,ntot
        real*8 :: nllk,llk,objgrd(10),asycov(10,10)


	!fn="ntrans5.dat"
	!fn="coli.txt"
	!nins=120
	!ntot=125
	
	!fn="ntrans5t.dat"
	!nins=920
	!ntot=1380
	!prd=92

	!!sunspot data
	fn="sunspots.dat"
	nins=440
	ntot=459
	prd=11

        open(unit=10,file=fn,status="old",iostat=ierr)
        if(ierr.ne.0) then
                print *, "Filename error."
                return
        end if
        allocate(obs(ntot))
        read(10,*) obs
        close(10)

	asycov=0.0d0
	para=0.0d0
        para(1:3)=(/0.5d0,0.4d0,0.2d0/)
        para(4:6)=(/0.5d0,0.4d0,0.2d0/)

        call ettpar(para,r,obs(1:nins),llk,'U')
	call dc_ttpar(para,r,obs(1:nins),asycov)
        write(*,*) "Estimated threshold:",r
        write(*,*) "Estimated LLK:",llk
        write(*,*) "Estimated AIC:", -2.0d0*llk*(nins-1)+2.0d0*11
        write(*,*) "Estimated paramaters and standard error:"
        write(*,'(2F10.3)') (para(i),sqrt(asycov(i,i)/size(obs)),i=1,10)
        call coverage_ttpar(para,r,nins,obs)

        end subroutine rd_ttpar

	subroutine coli_ttpar
        implicit none
        
	integer:: nobs
        integer,allocatable :: obs(:)
        integer :: r,ierr
        real*8 :: para(10)
        character(len=20) :: fn

	integer :: i,nins,ntot
        real*8 :: nllk,llk,objgrd(10),asycov(10,10)

	!fn="ntrans5.dat"
	fn="coli.txt"
	nins=120
	ntot=125
	prd=12
        open(unit=10,file=fn,status="old",iostat=ierr)
        if(ierr.ne.0) then
                print *, "Filename error."
                return
        end if
        allocate(obs(ntot))
        read(10,*) obs
        close(10)

	asycov=0.0d0
	para=0.0d0
        para(1:3)=(/0.5d0,0.4d0,0.2d0/)
        para(4:6)=(/0.5d0,0.4d0,0.2d0/)

        call ettpar(para,r,obs(1:nins),llk,'U')

	
	call dc_ttpar(para,r,obs(1:nins),asycov)
        write(*,*) "Estimated threshold:",r
        write(*,*) "Estimated LLK:",llk
        write(*,*) "Estimated AIC:", -2.0d0*llk*(nins-1)+2.0d0*11
        write(*,*) "Estimated paramaters and standard error:"
        write(*,'(2F10.3)') (para(i),sqrt(asycov(i,i)/size(obs)),i=1,10)
        call coverage_ttpar(para,r,nins,obs)

        end subroutine coli_ttpar


        end module ttpar_mod

