        module tparc_mod
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
        
        subroutine nllk_tparc(r,para,obs,nllk)
        !nllk=negative log-likilihood
        implicit none
        integer,intent(in) :: r
        real*8,intent(in) :: para(6)
        integer,intent(in) :: obs(:)
        real*8,intent(out) :: nllk

        integer :: i
        real*8 :: lambda

        nllk=0.0d0

	if(obs(1).gt.0.0d0) then
        	lambda=obs(1)
	else
		lambda=0.01d0
	end if


        do i=2,size(obs)
                if(obs(i-1).le.r) then
                        lambda=para(1)+para(2)*lambda+para(3)*obs(i-1)
                else
                        lambda=para(4)+para(5)*lambda+para(6)*obs(i-1)
                end if
                nllk=nllk+lambda-obs(i)*log(lambda)
        end do
        nllk=nllk/(size(obs)-1)
        end subroutine nllk_tparc

        subroutine grd_tparc(r,para,obs,grd)
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
        end subroutine grd_tparc

        subroutine objfun_tparc(MODE,N,X,OBJF,OBJGRD,NSTATE,IUSER,RUSER)
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
        call nllk_tparc(iuser(1),x,int(ruser),objf)
        call grd_tparc(iuser(1),x,int(ruser),objgrd)
        end subroutine objfun_tparc

             
        subroutine dc_tparc(para,r,obs,asycov)
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
	write(*, '(A,F10.3)') "Mean of (Y_t-lambda_t)^2", sum((obs-lambda)**2.0d0)/size(obs)
	write(*, '(A,F10.3)') "Mean of (Y_t-lambda_t)^2/lambda_t", sum((obs-lambda)**2.0d0/lambda)/size(obs)
	if(i.ne.0) print *,"Error in finding inv. of asycov" 
	open(unit=11,file="summary.dat")
	write(11,'(I,1X,F10.3)') (obs(i),lambda(i),i=1,size(obs))
	close(11)
	write(*,'(A)')  "Statistics of the Pearson residual:"
	mrsdl=sum(rsdl)/size(obs)
	sdrsdl=sqrt(sum((rsdl-mrsdl)**2.0d0)/size(obs))
	write(*, '(A,F10.3)') "mean", mrsdl
	write(*, '(A,F10.3)') "sd", sdrsdl
	write(*, '(A,F10.3)') "skewness", sum((rsdl-mrsdl)**3.0d0)/size(obs)/sdrsdl**3.0d0
	write(*, '(A,F10.3)') "ex. kur.", sum((rsdl-mrsdl)**4.0d0)/size(obs)/sdrsdl**4.0d0-3.0d0
        end subroutine dc_tparc


        subroutine etparc(para,r,obs,llk,th)  
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
        bl(1:7)=1.0d-4
	bl(6)=0.0d0
        bu(1)=1.0d4
        bu(2:3)=1.0d1
        bu(4)=1.0d4
        bu(5)=1.0d0
        bu(6)=0.0d0
        bu(7)=1.0d0
        
        if(th .eq. 'K') then
                iuser(1)=r
                iuser(2)=size(obs)
                !call e04uef('Verify Level = 3')
                call e04uef('Print Level = 0')
                call e04uef('Major Iteration Limit = 2000')
		ifail=0
                call e04ucf(6,1,0,1,1,6,a,bl,bu,e04udm,&
                   objfun_tparc,iter,istate,c,cjac,clamda,objf,objgrd,rmtrx,para,&
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
                           objfun_tparc,iter,istate,c,cjac,clamda,objf,objgrd,rmtrx,parath(i,:),&
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

        end subroutine etparc
        end module tparc_mod
