        include "../general_codes/summary.f90"

        module gtpar_mod
        use mat_ope_nag
        use nvst_int
        use gnuplot_int
        implicit none

        contains

	function outer_product(x,y)
	implicit none
	real*8 :: x(:),y(:)
	real*8 :: outer_product(size(x),size(y))

	integer :: i

	do i=1,size(x)
		outer_product(i,:)=x(i)*y
	end do
	end function outer_product

        subroutine sim_gtpar(r,p,q,para,obs)
        implicit none
        integer :: r,p,q
        real*8 :: para(2*(1+p+q))
        integer :: obs(:)

        INTEGER :: MSEED, MSTATE, M
        PARAMETER (MSEED=1,MSTATE=633,M=1)
        INTEGER :: GENID, I, IFAIL, LSEED, LSTATE, SUBID
        INTEGER :: SEED(MSEED), STATE(MSTATE), X(1)
        
        integer :: n,j
        integer :: y(int(size(obs)*1.5))

        real*8 :: lambda(int(size(obs)*1.5))
        
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

        y(1:max(p,q))=10
        lambda=1.0d0
	j=p+q+1
        do i=max(p,q)+1,n
                if(y(i-1).le.r) then
                        lambda(i)=para(1)+dot_product(para(2:p+1),lambda(i-1:i-p:-1))+&
			  dot_product(para(p+2:j),y(i-1:i-q:-1))
                else
			lambda(i)=para(j+1)+dot_product(para(j+2:j+p+1),lambda(i-1:i-p:-1))+&
                                  dot_product(para(j+p+2:2*j),y(i-1:i-q:-1))
                end if
                CALL G05TKF(M,Lambda(i:i),STATE,X,IFAIL)
                y(i)=x(1)
        end do
        obs=y(n-size(obs)+1:n)
	open(10,file="sim_data.dat")
	write(10,'(I10, F10.3)') (y(i),lambda(i),i=n-size(obs)+1,n)
	close(10)
        end subroutine sim_gtpar


	subroutine nllk_gtpar(r,p,q,para,obs,nllk)
	!nllk=negative log-likilihood
	implicit none
	integer,intent(in) :: r,p,q
	real*8,intent(in) :: para(2*(1+p+q))
	integer,intent(in) :: obs(:)
	real*8,intent(out) :: nllk

	integer :: i,idx,pmax
	real*8 :: jv(1+p+q)
	real*8 :: lambda(size(obs))

	pmax=max(p,q)

	where(obs(1:pmax).eq.0) 
		lambda(1:pmax)=0.1d0
	elsewhere
		lambda(1:pmax)=dble(obs(1:pmax))
	end where

	idx=p+q+1
	do i=pmax+1,size(obs)
                jv=(/1.0d0,lambda(i-1:i-p:-1),dble(obs(i-1:i-q:-1))/)
		if(obs(i-1).le.r) then
			lambda(i)=dot_product(para(1:idx),jv)
		else 
			lambda(i)=dot_product(para(idx+1:),jv)
		end if
	end do
	nllk=sum(lambda(pmax+1:)-obs(pmax+1:)*log(lambda(pmax+1:)))/dble(size(obs)-pmax)

	!!! The following is for testing only.
	!open(10,file="nllk_cal_detail.dat")
	!write(10,'(I10, F10.3)') (obs(i),lambda(i),i=1,size(obs))
	!close(10)
	!print *,"nllk_gtpar output:"
	!print *, "NLLK=",nllk
	end subroutine nllk_gtpar

        subroutine grd_gtpar(r,p,q,para,obs,grd)
        implicit none
        integer,intent(in) :: r,p,q
	real*8,intent(in) :: para(2*(1+p+q))
        integer,intent(in) :: obs(:)
	real*8,intent(out) :: grd(2*(1+p+q))

        integer :: i,j,idx,pmax
        real*8 :: lambda(size(obs))
	real*8 :: dlambda(size(obs),2*(1+p+q))
	real*8,dimension(1+p+q) :: jv

	!real*8 :: nllk

	pmax=max(p,q)
	idx=1+p+q

	where(obs(1:pmax).eq.0) 
		lambda(1:pmax)=0.1d0
	elsewhere
		lambda(1:pmax)=dble(obs(1:pmax))
	end where

        grd=0.0d0
	dlambda=0.0d0
	do i=pmax+1,size(obs)
                jv=(/1.0d0,lambda(i-1:i-p:-1),dble(obs(i-1:i-q:-1))/)
		if(obs(i-1).le.r) then
			lambda(i)=dot_product(para(1:idx),jv)
			dlambda(i,1:idx)=jv
			do j=1,p
				dlambda(i,:)=dlambda(i,:)+para(1+j)*dlambda(i-j,:)
			end do
		else 
			lambda(i)=dot_product(para(idx+1:),jv)
			dlambda(i,idx+1:)=jv
			do j=1,p
				dlambda(i,:)=dlambda(i,:)+para(idx+1+j)*dlambda(i-j,:)
			end do
		end if
		grd=grd+(1.0d0-obs(i)/lambda(i))*dlambda(i,:)
		!write( *,'(I5,20D10.2)') i, grd
		!if(modulo(i,100).eq.0) pause
	end do
        grd=grd/(size(obs)-pmax)
	!!!! The following is for testing only.
	!nllk=sum(lambda(pmax+1:)-obs(pmax+1:)*log(lambda(pmax+1:)))/dble(size(obs)-pmax)
	!print *, "grd_gtpar output:"
	!print *, "NLLK=",nllk
        !print *,"GRADIENT=", grd
	!open(10,file="grd_cal_detail.dat")
	!write(10,'(I10, F10.3,10F10.3)') (obs(i),lambda(i),dlambda(i,:), i=1,size(obs))
	!close(10)
        end subroutine grd_gtpar

        subroutine objfun_gtpar(MODE,N,X,OBJF,OBJGRD,NSTATE,IUSER,RUSER)
        implicit none

        integer,intent(inout) :: mode
        integer,intent(in) :: n
        real*8, intent(in) :: x(n)
        real*8,intent(out) :: objf
        real*8,intent(out) :: objgrd(n)
        integer,intent(in) :: nstate
        integer,intent(in) :: iuser(4)
        real*8,intent(in) :: ruser(iuser(4))

        mode=2
        call nllk_gtpar(iuser(1),iuser(2),iuser(3),x,int(ruser),objf)
        call grd_gtpar(iuser(1),iuser(2),iuser(3),x,int(ruser),objgrd)
        end subroutine objfun_gtpar

             
        subroutine dc_gtpar(r,p,q,para,obs,asycov,ninsample)
        implicit none
        integer :: r,p,q
        real*8 :: para(2*(1+p+q))
        integer :: obs(:)
       	real*8 :: asycov(2*(1+p+q),2*(1+p+q))
	integer,optional :: ninsample


        integer :: i,j,pmax,idx,nins
        real*8 :: lambda(size(obs)),rsdl(size(obs)), mrsdl,sdrsdl
	real*8 :: dlambda(size(obs),2*(1+p+q)),deri(2*(1+p+q)),jv(1+p+q)
       	real*8 :: asy(2*(1+p+q),2*(1+p+q))


	if(present(ninsample).and. ninsample .lt. size(obs)) then
		nins=ninsample
	else
		nins=size(obs)
	end if

	pmax=max(p,q)
	idx=1+p+q

	where(obs(1:pmax).eq.0) 
		lambda(1:pmax)=0.1d0
	elsewhere
		lambda(1:pmax)=dble(obs(1:pmax))
	end where

!        rsdl=0.0d0
	dlambda=0.0d0
	asy=0.0d0

	do i=pmax+1,size(obs)
                jv=(/1.0d0,lambda(i-1:i-p:-1),dble(obs(i-1:i-q:-1))/)
		if(obs(i-1).le.r) then
			lambda(i)=dot_product(para(1:idx),jv)
			dlambda(i,1:idx)=jv
			do j=1,p
				dlambda(i,:)=dlambda(i,:)+para(1+j)*dlambda(i-j,:)
			end do
		else 
			lambda(i)=dot_product(para(idx+1:),jv)
			dlambda(i,idx+1:)=jv
			do j=1,p
				dlambda(i,:)=dlambda(i,:)+para(idx+1+j)*dlambda(i-j,:)
			end do
		end if
		asy=asy+1.0d0/lambda(i)*outer_product(dlambda(i,:),dlambda(i,:))
		!grd=grd+(1.0d0-obs(i)/lambda(i))*dlambda(i,:)
		!write( *,'(I5,20D10.2)') i, grd
		!if(modulo(i,100).eq.0) pause
	end do

	asy=asy/dble(size(obs)-pmax)
        call sub_inv_pd_mat(asy, asycov,i)

	write(*,'(A)') "The standard errors of parameter estimates:"
	write(*,'(I,2F10.3)') (i,para(i),sqrt(asycov(i,i)),i=1,2*(1+p+q))

	if(nins.lt.size(obs)) then
		write(*,'(A,I0,A)') "The first ",nins, " observations are in-sample:"
		write(*, '(A,F10.3)') "Average of (Y_t-lambda_t)^2", sum((obs(1:nins)-lambda(1:nins))**2.0d0)/nins
		write(*,'(A,I0,A)') "The last ",nins, " observations are out-of-sample:"
		write(*, '(A,F10.3)') "Average of (Y_t-lambda_t)^2", sum((obs(nins+1:)-lambda(nins+1:))**2.0d0)/(size(obs)-nins)
	else
		write(*,'(A)') "All data are used to estimate parameters, only the in-sample MSE of forecasts will calculated."
		write(*, '(A,F10.3)') "Average of (Y_t-lambda_t)^2", sum((obs(1:nins)-lambda(1:nins))**2.0d0)/nins

	end if

!	if(i.ne.0) print *,"Error in finding inv. of asycov" 
!	open(unit=11,file="summary.dat")
!	write(11,'(I,1X,F10.3)') (obs(i),lambda(i),i=1,size(obs))
!	close(11)
!	write(*,'(A)')  "Statistics of the Pearson residual:"
!	mrsdl=sum(rsdl)/size(obs)
!	sdrsdl=sqrt(sum((rsdl-mrsdl)**2.0d0)/size(obs))
!	write(*, '(A,F10.3)') "mean", mrsdl
!	write(*, '(A,F10.3)') "sd", sdrsdl
!	write(*, '(A,F10.3)') "skewness", sum((rsdl-mrsdl)**3.0d0)/size(obs)/sdrsdl**3.0d0
!	write(*, '(A,F10.3)') "ex. kur.", sum((rsdl-mrsdl)**4.0d0)/size(obs)/sdrsdl**4.0d0-3.0d0
        end subroutine dc_gtpar


        subroutine egtpar(r,p,q,para,obs,llk,th)  
        implicit none
        integer :: r,p,q
        real*8 :: para(2*(1+p+q))
        integer :: obs(:)
        real*8 :: llk
        character*1 :: th

	integer :: idx
        real*8 :: a(2*(1+p+q)),bl(2*(1+p+q)+1),bu(2*(1+p+q)+1),c(1),cjac(1,1),clamda(2*(1+p+q)+1),&
                  objf,objgrd(2*(1+p+q)),rmtrx(2*(1+p+q),2*(1+p+q))
        integer :: iter,istate(2*(1+p+q)),iwork(6*(1+p+q)+1),ifail,iuser(4),i
      	integer :: lwork
	real*8,dimension(:),allocatable :: work
        real*8 :: alpha(2),qv(2)
        integer :: qtl(2)

        real*8, allocatable :: parath(:,:),llkth(:)
        external e04udm


	idx=1+p+q
	lwork=8*(1+p+q)**2+40*(1+p+q)+11
        alpha=(/0.2d0,0.80d0/)
        a=0.0d0
        a(p+q+2:)=1.0d0
        bl=0.001d0
        bu(1)=1.0d2
        bu(2:idx)=1.0d1
        bu(idx+1)=1.0d2
        bu(idx+2:2*idx)=1.0d0
        bu(2*idx+1)=1.0d0
        allocate(work(lwork)) 
       
        if(th .eq. 'K') then
                iuser(1)=r
                iuser(2)=p
		iuser(3)=q
		iuser(4)=size(obs)
                !call e04uef('Verify Level = 3')
                !call e04uef('Print Level = 3')
                call e04uef('Major Iteration Limit = 2000')
		ifail=0
                call e04ucf(2*(1+p+q),1,0,1,1,&
		   2*(1+p+q),a,bl,bu,e04udm,&
                   objfun_gtpar,iter,istate,c,cjac,clamda,objf,objgrd,rmtrx,para,&
                   iwork,6*(1+p+q)+1,work,lwork,iuser,dble(obs),ifail)
                llk=-objf
        else
                call g01amf(size(obs),dble(obs),2,alpha,qv,ifail)
                if(ifail.ne.0) then 
                        print *,"Error in finding empirical quantiles &
                                & of observations."
                        pause
                end if
                qtl=int(qv)
                write(*,*) "Since threshold is unknown, it will &
                             & be searched."
                write(*,*) "From ",qtl(1), " to ", qtl(2)
                allocate(parath(qtl(1):qtl(2),2*(1+p+q)))
                allocate(llkth(qtl(1):qtl(2)))
                write(*,*) "Threshold Log-likilihood & 
                              &  Current_best_threshold" 
                do i=qtl(1),qtl(2)
                        parath(i,:)=para
                        iuser(1)=i
                        iuser(2)=p
			iuser(3)=q
			iuser(4)=size(obs)
                        !call e04uef('Verify Level = 3')
                        call e04uef('Print Level = 0')
                        call e04uef('Major Iteration Limit = 2000')
			ifail=0
                        call e04ucf(2*(1+p+q),1,0,1,1,&
			   2*(1+p+q),a,bl,bu,e04udm,&
			   objfun_gtpar,iter,istate,c,cjac,clamda,objf,objgrd,rmtrx,para,&
			   iwork,6*(1+p+q)+1,work,lwork,iuser,dble(obs),ifail)
                        llkth(i)=-objf
                        if( i.eq.qtl(1)) then
                                r=i
                        else
                                if(llkth(i).gt.llkth(r)) r=i
                        end if
                        write(*,'(4X,I, 5X,F12.5,2X, I )') ,i,llkth(i),r
                 end do
                 para=parath(r,:)
                 llk=llkth(r)
                 deallocate(parath,llkth)
        end if
	
	deallocate(work)
        end subroutine egtpar

	subroutine simulation1
	implicit none

        integer,parameter :: p=2,q=2,repl=1000,nobs=3000
	integer,parameter :: idx=1+p+q
	integer,parameter :: npara=2*(1+p+q)
        integer :: r,r0
        real*8 :: para(npara),para0(npara),epara(repl,npara),er(repl)
        integer :: obs(nobs)

        real*8 :: nllk,llk,objgrd(npara)
        integer :: nstate,iuser(3),mode,i


        !r=7
        !para(1:3)=(/0.5d0,0.8d0,0.4d0/)
        !para(4:6)=(/0.2d0,0.7d0,0.1d0/)
        !r=6
        !para(1:3)=(/0.33d0,0.91d0,0.001d0/)
        !para(4:6)=(/0.11d0,0.83d0,0.15d0/)
	r=3
        para(1:idx)=(/0.6d0,0.3d0,0.1d0,0.2d0,0.1d0/)
        para(idx+1:)=(/0.8d0,0.4d0,0.2d0,0.2d0,0.1d0/)

        r0=r
        para0=para

        do i=1,repl
		write(*,'(A,I0)') "Replication number:",i
                call sim_gtpar(r0,p,q,para0,obs)

		
		!!For testing only
		!call nllk_gtpar(r,p,q,para,obs,nllk)
        	!call grd_gtpar(r,p,q,para,obs,objgrd)
		!call gplot1(dble(obs))
		!!

                para=para0
		r=r0
                call egtpar(r,p,q,para,obs,llk,'K')
                epara(i,:)=para
                er(i)=r
        end do
	write(*,*) "SUMMARY OF SIMULATION RESULTS:"
        print *,"sample size=",nobs
        print *,"replication=", repl
        print *,"True parameters"
        write(*,*) "Threshold:", r0
        write(*,*) "Other para:"
        write(*, '(5F12.3)') para0
        print *, "Average of estimated parameter:"
        write(*, '(5F12.3)') sum(epara,dim=1)/repl
        print *, "Covariance matrix:"
        call print_mat(covm(epara))

        write(*,*) "Estimated threshold:" 
        write(*,*) "Mean", mean(dble(er))
        write(*,*) "Variance", var(dble(er))
        open(unit=11,file="er.dat",status="replace")
        do i=1,repl
                write(11,*) er(i)
        end do
        close(11)
        !call gplot1(dble(er))
        end subroutine simulation1
        
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       ! subroutine test_asycov
       ! implicit none
       ! integer,parameter :: repl=1000,nobs=3000
       ! integer :: r,r0
       ! real*8 :: para(6),para0(6),epara(repl,6),er(repl)
       ! integer :: obs(nobs)
!
!        real*8 :: nllk,llk,objgrd(6),asycov(repl,6,6)
!        integer :: nstate,iuser(2),mode,i
!
!
!        r=7
!        para(1:3)=(/0.5d0,0.8d0,0.4d0/)
!        para(4:6)=(/0.2d0,0.7d0,0.1d0/)
!        
!	!r=6
!	!para(1:3)=(/0.3d0,0.8d0,0.1d0/)
!        !para(4:6)=(/0.1d0,0.7d0,0.25d0/)
!        
!	r0=r
!        para0=para
!
!        do i=1,repl
!		print *,"index repl=", i
!                call sim_tpar(r0,para0,obs)
!                para=para0
!                call etpar(para,r,obs,llk,'U')
!		call dc_tpar(para,r,obs,asycov(i,:,:))
!                epara(i,:)=para
!                er(i)=r
!        end do
!        print *,"sample size=",nobs
!        print *,"replication=", repl
!        print *,"True parameters"
!        write(*,*) "Threshold:", r0
!        write(*,*)"Other para:"
!        write(*, '(3F12.3)') para0
!        print *, "Average of estimated parameter:"
!        write(*, '(3F12.3)') sum(epara,dim=1)/repl
!        print *, "Covariance matrix:"
!	print *, "T*SV:"
!        call print_mat(nobs*covm(epara))
!	print *, "SM of G:"
!	call print_mat(sum(asycov,dim=1)/repl)
!        
!	write(*,*) "Estimated threshold:" 
!        write(*,*) "Mean", mean(dble(er))
!        write(*,*) "Variance", var(dble(er))
!        open(unit=11,file="er.dat",status="replace")
!        do i=1,repl
!                write(11,'(10F10.3)') er(i), epara(i,:)
!        end do
!        close(11)
!        !call gplot1(dble(er))
!        end subroutine test_asycov
!

!        subroutine simulation
!        implicit none
!
!        integer,parameter :: repl=1000,nobs=1000
!        integer :: r
!        real*8 :: para(6),para0(6),epara(repl,6)
!        integer :: obs(nobs)
!
!        real*8 :: nllk,llk,objgrd(6)
!        integer :: nstate,iuser(2),mode,i
!
!
!        !r=4
!        !para(1:3)=(/0.3d0,0.9d0,0.4d0/)
!        !para(4:6)=(/1.0d0,0.5d0,0.1d0/)
!        r=8
!        para(1:3)=(/0.30d0,0.88d0,0.05d0/)
!        para(4:6)=(/0.14d0,0.73d0,0.23d0/)
!
!
!        para0=para
!        do i=1,repl
!                call sim_tpar(r,para0,obs)
!                para=para0
!                call etpar(para,r,obs,llk,'K')
!                !pause  
!                epara(i,:)=para
!        end do
!        print *,"sample size=",nobs
!        print *,"replication=", repl
!        print *,"True parameters"
!        write(*, '(30E12.3)') para0
!        print *, "Average of estimated parameter:"
!        write(*, '(30E12.3)') sum(epara,dim=1)/repl
!        print *, "Covariance matrix:"
!
!        call print_mat(covm(epara))
!        end subroutine simulation
!        
!        subroutine coverage_tpar(para,r,nins,obs)
!	implicit none
!	integer,intent(in) :: r
!        real*8,intent(in) :: para(6)
!	integer :: nins
!!        integer,intent(in) :: obs(:)
!
!        integer :: i
!        real*8 :: lambda(size(obs))
!	
!	if(obs(1).gt.0.0d0) then
!        	lambda(1)=obs(1)
!	else
!		lambda(1)=0.01d0
!	end if
!
!        do i=2,size(obs)
!                if(obs(i-1).le.r) then
!                        lambda(i)=para(1)+para(2)*lambda(i-1)+para(3)*obs(i-1)
!                else
!                        lambda(i)=para(4)+para(5)*lambda(i-1)+para(6)*obs(i-1)
!                end if
!        end do
!	open(unit=11,file="tpar_out_of_sample_forecasting.dat")
!	write(11,'(I5,1X,F10.3)') (obs(i),lambda(i),i=1,size(obs))
!	close(11)
!	write(*,'(A)') "For the TPoiAR model:"
!	write(*,'(A,D15.5)') "in-sample mse:",sum((obs(1:nins)-lambda(1:nins))**2.0d0)/dble(nins)
!	write(*,'(A,D15.5)') "out-of-sample mse:",sum((obs(nins+1:)-lambda(nins+1:))**2.0d0)/dble(size(obs)-nins)
!        end subroutine coverage_tpar
!        
       

	subroutine rd_gtpar
        implicit none

	integer,parameter :: p=2,q=2
        integer,allocatable :: obs(:)
        integer :: r,ierr
        real*8 :: para(2*(1+p+q))
        character(len=20) :: fn

	integer :: i,nins,ntot
        real*8 :: nllk,llk,objgrd(2*(1+p+q)),asycov(2*(1+p+q),2*(1+p+q)),parapar(3)


        !write(*,*) "Please type name of data file."
        !read(*,*) fn
	fn="ntrans5t.dat"
	!fn="coli.txt"
	nins=920
	ntot=1380
	
        open(unit=10,file=fn,status="old",iostat=ierr)
        if(ierr.ne.0) then
                print *, "Filename error."
                return
        end if
        allocate(obs(ntot))
        read(10,*) obs
        close(10)
!
        para(1:3)=(/0.71d0,0.22d0,0.41d0/)
        para(4:6)=(/0.71d0,0.22d0,0.41d0/)
	!! use the first nins data to estimate parameters
        call egtpar(r,p,q,para,obs(1:nins),llk,'U')
	print *, "LLK=",llk
	call dc_gtpar(r,p,q,para,obs,asycov,ninsample=nins)
	write(*,*) "Estimated threshold:",r
	
        deallocate(obs)
        end subroutine rd_gtpar
!
        end module gtpar_mod

        program main_program
        use gtpar_mod
        implicit none
        
	!call test_asycov
        !call simulation_threshold
        !call simulation1
        call rd_gtpar
        end program main_program

