	module par_mod
	use mat_ope_nag
	use nvst_int
	implicit none

	!integer,parameter :: repl=1000
	!integer,parameter :: nobs=1000,l=8
	!integer :: obs(nobs)
	!real*8 :: p0(3)
	!save

	contains
		
	!subroutine simpa11(p)
	!implicit none
	!real*8,intent(in) :: p(3)
	!integer :: yt,i
	!!real*8 :: lambda,d,a,b
	!integer,external :: zbqlpoi
	!external zbqlini
!
!	d=p(1)
!!	a=p(2)
!	b=p(3)
!
!	call zbqlini(0)
!	
!	lambda=d/(1.0d0-a-b)
!	yt=int(lambda)
!	do i=1,nobs/2
!		lambda=d+a*lambda+b*yt
!		yt=zbqlpoi(lambda)
!	end do
!
!	obs(1)=yt
!	do i=2,nobs
!		lambda=d+a*lambda+b*obs(i-1)
!		obs(i)=zbqlpoi(lambda)
!	end do
!	end subroutine simpa11

!	subroutine simpa
!	implicit none
!	integer,parameter :: np=1,nq=2
!	real*8 :: p(0:np+nq)
!	integer :: yt(int(1.5*nobs)),i
!	real*8 :: lambda(int(1.5*nobs))
!
!	integer,external :: zbqlpoi
!	external zbqlini
!
!	call zbqlini(0)
!
!		
!	!p=(/0.3d0,0.4d0,0.5d0/)
!	p=(/0.2d0,0.4d0,0.3d0,0.1d0/)
!	lambda(1:1+np+nq)=p(0)
!	yt(1:1+np+nq)=3
!	do i=1+np+nq+1,int(1.5*nobs)
!		lambda(i)=p(0)+dot_product(p(1:np),lambda(i-1:i-np:-1))+dot_product(p(np+1:np+nq),yt(i-1:i-nq:-1))
!		yt(i)=zbqlpoi(lambda(i))
!	end do
!	obs=yt(int(1.5*nobs)-nobs+1:)
!	end subroutine simpa

	!p=(d,a,b)
	subroutine nllk_par(p,obs,nllk)
	implicit none
	real*8, intent(in) :: p(:)
	integer,intent(in) :: obs(:)
	real*8, intent(out) :: nllk

	real*8 :: lambda(size(obs))
	integer :: i,nobs

	nobs=size(obs)

	lambda(1)=obs(1)
	do i=2,nobs
		lambda(i)=p(1)+p(2)*lambda(i-1)+p(3)*obs(i-1)
	end do
	nllk=-sum(log(lambda(2:))*dble(obs(2:))-lambda(2:))/dble(size(obs)-1)

	end subroutine nllk_par
	
	subroutine nllk_gc_par(p,obs,nllk,gc)
	implicit none
	real*8, intent(in) :: p(3)
	integer,intent(in) :: obs(:)
	real*8, intent(out) :: nllk,gc(3)

	real*8 :: lambda(size(obs)),dl(3),jv(3)
	integer :: i
	
	lambda(1)=obs(1)

	gc=0.0d0
	dl=0.0d0
	do i=2,size(obs)
		lambda(i)=p(1)+p(2)*lambda(i-1)+p(3)*obs(i-1)
		jv=(/1.0d0,lambda(i-1),dble(obs(i-1))/)
		dl=jv+p(2)*dl
		gc=gc+(dble(obs(i))/lambda(i)-1.0d0)*dl
	end do
	gc=-gc/(size(obs)-1)
	nllk=-sum(log(lambda(2:))*dble(obs(2:))-lambda(2:))/(size(obs)-1)
	end subroutine nllk_gc_par

	subroutine obj_par_wg(n,xc,fc,gc,ia,ra) !p,obs,nllk,gc)
	implicit none
	integer, intent(in) :: n
	real*8,intent(in) :: xc(n)
	real*8, intent(out) :: fc,gc(3)
	integer :: ia(1)
	real*8 :: ra(ia(1))

	call nllk_gc_par(xc,int(ra),fc,gc)
	end subroutine obj_par_wg

        subroutine dc_par(para,obs,asycov)
        implicit none
        real*8 :: para(3)
        integer :: obs(:)
       	real*8 :: asycov(3,3) 

        integer :: i
        real*8 :: lambda(size(obs)),rsdl(size(obs)), mrsdl,sdrsdl
	real*8 :: dlambda(size(obs),3),deri(3)
       	real*8 :: asy(3,3) 

	if(obs(1).gt.0.0d0) then
        	lambda(1)=obs(1)
	else
		lambda(1)=0.01d0
	end if
        rsdl=0.0d0
	dlambda=0.0d0
	asy=0.0d0
        do i=2,size(obs)
		lambda(i)=para(1)+para(2)*lambda(i-1)+para(3)*obs(i-1)
		deri=(/1.0d0,lambda(i-1),dble(obs(i-1))/)
		dlambda(i,:)=deri+para(2)*dlambda(i-1,:)
		asy=asy+1.0d0/lambda(i)*outer_product(dlambda(i,:),dlambda(i,:))
                rsdl(i)=(obs(i)-lambda(i))/sqrt(lambda(i))
        end do
	asy=asy
        call sub_inv_pd_mat(asy, asycov,i)
	if(i.ne.0) print *,"Error in finding inv. of asycov" 
        end subroutine dc_par


!	Diagnostic checking of PA11
!	subroutine dc_pa11(p,l,alpha,ac,rslt,ifail)
!	implicit none
!		real*8, intent(in) :: p(3)
!		real*8, intent(in) :: alpha
!		integer,intent(in) :: l !=lag
!		real*8, intent(out) :: ac(l)
!		character*1,intent(out) :: rslt
!		integer,optional :: ifail
!
!		real*8 :: lambda(nobs),res(nobs),dl(nobs,3),jv(3),c(l)
!		real*8 :: g(3,3),g_inv(3,3),asymv(l,l),asymv_inv(l,l),x(l,3),t
!		integer :: i,j
!		real*8,external :: g01ecf
!
!		g=0.0d0
!		x=0.0d0
!		lambda(1)=p(1)/(1.0d0-p(2)-p(3))
!		res(1)=(obs(1)-lambda(1))/sqrt(lambda(1))
!		dl=0.0d0
!		do i=2,nobs
!			lambda(i)=p(1)+p(2)*lambda(i-1)+p(3)*obs(i-1)
!			res(i)=(obs(i)-lambda(i))/sqrt(lambda(i))
!			jv=(/1.0d0,lambda(i-1),dble(obs(i-1))/)
!			dl(i,:)=jv+p(2)*dl(i-1,:)
!			g=g+vxvt(dl(i,:))/lambda(i)
!			do j=1,l
!				if(i-j>0) x(j,:)=x(j,:)+res(i-j)*dl(i,:)/sqrt(lambda(i))
!			end do
!		end do
!		g=g/dble(nobs-1)
!		call sub_inv_pd_mat(g,g_inv,ifail)
!		g=g_inv
!		!call print_mat(g)
!		!pause
!		do i=1,l
!			x(i,:)=x(i,:)/dble(nobs-i)
!		end do
!		call acf(res,l,ac)
!		asymv=matmul(matmul(x,g),transpose(x))
!		asymv=unitmat(l)-asymv
!		!call print_mat(asymv)
!		!pause
!		!asymv= inv_pd_mat(asymv)
!		call sub_inv_pd_mat(asymv,asymv_inv,ifail)
!		asymv=asymv_inv
!		t=dot_product(matmul(asymv,ac),ac)*dble(nobs)
!		t=g01ecf("U",t,dble(l),i)
!		if(i .ne.0) print *,"ERROR IN CALCULATING CHI-SQUARE PROB."
!		if(t<alpha) then
!			rslt="F"
!		else
!			rslt="T"
!		end if
!	end subroutine dc_pa11
!
	!Diagnostic checking of PA11
!	subroutine dc_pa11(p,l,alpha,ac,rslt,ifail)
!	implicit none
!		real*8, intent(in) :: p(3)
!		real*8, intent(in) :: alpha
!		integer,intent(in) :: l !=lag
!		real*8, intent(out) :: ac(l)
!		character*1,intent(out) :: rslt
!		integer,optional :: ifail
!
!		real*8 :: lambda(nobs),res(nobs),dl(nobs,3),jv(3),c(l)
!		real*8 :: g(3,3),g_inv(3,3),asymv(l,l),asymv_inv(l,l),x(l,3),t
!		integer :: i,j
!		real*8,external :: g01ecf
!
!		g=0.0d0
!		x=0.0d0
!		lambda(1)=p(1)/(1.0d0-p(2)-p(3))
!		res(1)=(obs(1)-lambda(1))/sqrt(lambda(1))
!		dl=0.0d0
!		do i=2,nobs
!			lambda(i)=p(1)+p(2)*lambda(i-1)+p(3)*obs(i-1)
!			res(i)=(obs(i)-lambda(i))/sqrt(lambda(i))
!			jv=(/1.0d0,lambda(i-1),dble(obs(i-1))/)
!			dl(i,:)=jv+p(2)*dl(i-1,:)
!			g=g+vxvt(dl(i,:))/lambda(i)
!			do j=1,l
!				if(i-j>0) x(j,:)=x(j,:)+res(i-j)*dl(i,:)/sqrt(lambda(i))
!			end do
!		end do
!		g=g/dble(nobs-1)
!		call sub_inv_pd_mat(g,g_inv,ifail)
!		g=g_inv
!		!call print_mat(g)
!		!pause
!		do i=1,l
!			x(i,:)=x(i,:)/dble(nobs-i)
!		end do
!		call acf(res,l,ac)
!		asymv=matmul(matmul(x,g),transpose(x))
!		asymv=unitmat(l)-asymv
!		!call print_mat(asymv)
!		!pause
!		!asymv= inv_pd_mat(asymv)
!		call sub_inv_pd_mat(asymv,asymv_inv,ifail)
!		asymv=asymv_inv
!		t=dot_product(matmul(asymv,ac),ac)*dble(nobs)
!		t=g01ecf("U",t,dble(l),i)
!		if(i .ne.0) print *,"ERROR IN CALCULATING CHI-SQUARE PROB."
!		if(t<alpha) then
!			rslt="F"
!		else
!			rslt="T"
!		end if
!	end subroutine dc_pa11
!
	SUBROUTINE obj_par(N, XC, FC, IUSER, rUSER)
	implicit none
	integer :: n
	real*8 :: xc(n),fc
	integer :: iuser(1)
	real*8 :: ruser(iuser(1))

	call nllk_par(xc,int(ruser),fc)

	end subroutine obj_par

	!estimate pa11
	subroutine epar(para,obs,llk)
	implicit none

	real*8,intent(inout) :: para(3)
	integer,intent(in) :: obs(:)
	real*8,intent(out) :: llk

	real*8 :: bl(3),bu(3),w(39)
	integer :: iuser(1),iw(5),ifail

	bl=1.0d-5
	bu=1.0d0
	bu(1)=1.0d2

	ifail=-1
	iuser=size(obs)
	call e04jyf(3,0,obj_par,bl,bu,&
	            para,llk,iw,5,w,39,iuser,dble(obs),ifail)
	llk=-llk
	end subroutine epar

!estimate pa11 with gradient
	subroutine epar_wg(p,obs,llk)
	implicit none
	real*8,intent(inout) :: p(3)
	integer,intent(in) :: obs(:)
	real*8, intent(out) :: llk

	real*8 :: bl(3),bu(3),nllk,gc(3),w(30)
	integer :: ia(1),iw(5),ifail
	
	bl=1.0d-5
	bu=1.0d0
	bu(1)=1.0d2

	p=0.1d0
	ifail=-1
	ia=size(obs)
	call e04kzf(3,0,obj_par_wg,bl,bu,&
	            p,nllk,gc,iw,5,&
		    w,30,ia,dble(obs),ifail)
	llk=-nllk
	print *, "ESTIMATED PARAMETERS USING ANALYTIC GRADIENTS:"
	print *, p
	write(*,'(A,F12.5)') "averaged LLK=",llk
	!if(ifail.ne.0) pause
	end subroutine epar_wg
!!

!	subroutine main_par
!	implicit none
!	
!	real*8, parameter :: alpha=0.05
!	character*1 :: rslt
!	integer :: irslt(repl)
!	integer :: ia(1),i,ifail
!	real*8 :: p(3),v,ra(1),ac(repl,l)
!
!	real*8 :: pr(repl,3)
!	
!	open(unit=10,file="par_hat.dat")
!	open(unit=11,file="ac_res.dat")
!
!	p0=(/0.3d0,0.4d0,0.5d0/)
!	irslt=0
!	i=1
!	do
!		print *,"i=",i
!		!call simpa11(p0)
!		call simpa
!		print *, "DATA SIMULATED."
!		!open(unit=10,file="sim.dat")
!		!write(10, fmt='(I)') obs
!		!close(10)
!		p=(/0.2d0,0.40d0,0.3d0/)
!		call epa11(p,ifail)
!		print *, "ESTMIATION FINISHED."
!		print *, p
!		if(ifail.ne.0) cycle
!		pr(i,:)=p
!		call dc_pa11(p,l,alpha,ac(i,:),rslt,ifail)
!		if(ifail.ne.0) cycle
!		if(rslt.eq."F") irslt(i)=1
!		print *, "DIAGNOSITIC CHECKING FINISHED:",rslt
!		i=i+1
!		if(i>repl) exit
!	end do
!	close(10)
!	close(11)
!
!	!call print_mat(nobs*covm(pr))
!
!
!	print *,"NOBS=",nobs
!	print *, "REPL=", repl
!
!	print *, "SAMPLE MEAN OF ESTIMATED PARAMETERS."
!	print *, meanv(pr)
!	print *, "SAMPLE COVARIANCE MATRIX OF ESTIMATED PARAMETERS."
!	call print_mat(nobs*covm(pr))
!
!	print *, "SAMPLE MEAN OF SACV OF ESTIMATED RESIDUALS."
!	print *, meanv(ac)
!	print *, "SAMPLE COVARIANCE MATRIX OF SACV OF ESTIMATED RESIDUALS."
!	call print_mat(nobs*covm(ac))
!
!	print *, "EMPIRICAL SIZE:", sum(irslt)/dble(repl)
!
!
!	end subroutine main_par
!
	subroutine epa11_by_R(p,obs,llk)
	implicit none
	real*8,intent(inout) :: p(3)
	integer,intent(in) :: obs(:)
	real*8,intent(out) :: llk

	open(unit=10,file="obs.dat")
	write(10,'(I0)') obs
	close(10)

	call system("R --vanilla <epa11.R> epa11.out")
	open(unit=11,file="epa11_par.dat")
	read(11,*) p
	close(11)

	open(unit=11,file="epa11_llk.dat")
	read(11,*) llk
	close(11)
	
	end subroutine epa11_by_R

	subroutine coverage_par(para,nins,obs)
	implicit none
        real*8,intent(in) :: para(3)
	integer,intent(in) :: nins
        integer,intent(in) :: obs(:)

        integer :: i
        real*8 :: lambda(size(obs)),llk
	
	if(obs(1).gt.0) then
        	lambda(1)=obs(1)
	else
		lambda(1)=0.01d0
	end if

        do i=2,size(obs)
                        lambda(i)=para(1)+para(2)*lambda(i-1)+para(3)*obs(i-1)
        end do
	llk=sum(-lambda(2:nins)+dble(obs(2:nins))*log(lambda(2:nins)))/dble(nins-1)
	open(unit=11,file="PoiAR_obs_forecasts_residuals.dat")
	write(11,'(I5,1X,2F10.3)') (obs(i),lambda(i),(obs(i)-lambda(i))/sqrt(lambda(i)),i=1,size(obs))
	close(11)
	write(*,'(A)') "===OUTPUT OF COVERAGE_PAR==="
	write(*,'(A)') "For the PoiAR model:"
	write(*,'(A,F15.5)') "LLK=",llk
	write(*,'(A,F20.5)') "AIC=", -2.0d0*llk*(nins-1)+2.0d0*3
	write(*,'(A,F20.5)') "BIC=", -2.0d0*llk*(nins-1)+log(dble(nins-1))*3
	write(*,'(A,I0)') "in-sample no. of obs:",nins
	write(*,'(A,F15.5)') "in-sample mse:",sum((obs(1:nins)-lambda(1:nins))**2.0d0)/dble(nins)
	write(*,'(A,I0)') "out-of-sample no. of obs:",size(obs)-nins
	write(*,'(A,F15.5)') "out-of-sample mse:",sum((obs(nins+1:)-lambda(nins+1:))**2.0d0)/dble(size(obs)-nins)
	write(*,'(A)') "END===OUTPUT OF COVERAGE_PAR===END"
	end subroutine coverage_par
 
end module par_mod
