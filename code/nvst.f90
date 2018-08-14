!!!!!!!!!!!!!!
!Date   : 2011.03.13
!Author : Wang, Chao Charles 

module nvst_int
implicit none

type stock
        integer::code
        character(len=12)::name
end type stock

type date
        integer::y,m,d
end type date

type stkquote
        character(len=8) :: date
        character(len=6) :: time
        real*8 :: julian
        real*8 :: price
        integer :: volume
end type stkquote
 
        
type daily
        !type(date) :: dat
        integer :: julian
        real*8 :: ope, hig, low, clo
        !integer :: vol
        real*8 :: vol
        real*8 :: adj
end type daily

type daily_pointer
        type(daily),pointer :: pointer_to_daily
end type daily_pointer

type daily_node
        type(daily) :: detail
        type(daily_node), pointer :: next
end type daily_node

!Specification of directories
character(len=100),parameter :: intraday_dir="/home/chao/workspace/investment/intraday_data/"
character(len=100),parameter :: daily_dir="/home/chao/workspace/investment/daily_data/"
character(len=100),parameter :: alert_dir="/home/chao/workspace/investment/alert/"
character(len=100),parameter :: code_dir="/home/chao/workspace/investment/code/"
character(len=100),parameter :: bak_dir=intraday_dir//"bak_files"

character(len=30),parameter :: quote_fmt='(A8,A6,1X,F20.10,1X,F10.4,I15)'
character(len=27),parameter :: alert_record_fmt='(A8,A6,1X,A4,F10.4,1X,A1)'
character(len=40),parameter :: alert_record_log_fmt=&
                       '(A8,A6,1X,A4,F10.4,1X,A1,F10.4,F7.2,"%")'
integer, parameter :: rknd=kind(1.0d0)
type single_record
                integer :: stkcode,stk_id
                real :: price
                integer :: shares
                character :: trad_type
                integer :: time,date
                integer :: trad_curr
        end type single_record

interface check_order
        module procedure check_order_d, check_order_i
end interface       

contains
	

        subroutine count_no_lines(f,no_lines)
        implicit none

                character(len=*) :: f
                integer :: no_lines

                character(len=80) :: sys_cmd=''

                write(sys_cmd,fmt='("wc ",A," > tmp_wc_file")') trim(f)
                call system(sys_cmd)
                open(unit=10,file="tmp_wc_file",status="old")
                read (10,*) no_lines
                close(10)
                call system("rm -f tmp_wc_file")
        end subroutine count_no_lines

        subroutine read_trade_record(f,ncrd,rcd)
        implicit none

                character(len=*) :: f
                integer :: ncrd
                type(single_record) :: rcd(ncrd)

                !integer, parameter :: max_no=3000000
                !character(len=40) :: sys_cmd
                !type(single_record) :: t_rcd(max_no)
                !integer :: n,ios,no_lines

                1000 format(I5,I9,F8.3,I11,A1,I6,I8,I2)

                open(unit=10,file=f,status="old")
                read(10, 1000) rcd
                close(10)
        end subroutine read_trade_record

        ! analyse realized volatility
        subroutine arv(qt)
        implicit none
        type(stkquote) :: qt(:)
        
        integer :: n,idx1(size(qt))

      
        idx1=int(qt%julian)
        n=int(qt(size(qt))%julian)-int(qt(1)%julian)
        end subroutine arv


        !one-day realized-volatility by ZMA adjusted
        subroutine odrvzma(prices,rv,std)
        implicit none
        real*8,intent(in) :: prices(:)
        real*8,intent(out) :: rv,std

        integer :: n,k,i,j,nb
        real*8 :: rva,c,veps,rq
        real*8 :: r(size(prices)-1),p(size(prices))
        real*8,allocatable :: rvsg(:)
        
        n=size(prices)
        p=prices
        r=log(p(2:n))-log(p(1:n-1))
        n=n-1
        rva=sum(r**2.0d0)
        print *,"rva=",rva
        veps=rva/(2.0d0*n)
        print *,"veps=",veps
        rq=sum(r**4.0d0)*n/3.0d0
        print *,"rq=",rq
        c=(rq/12.0d0/veps**2.0d0)**(-1.0d0/3.0d0)
        print *,"c=",c
        k=int(c*n**(2.0d0/3.0d0))
        print *,"k=",k
        allocate(rvsg(k))
        rvsg=0.0d0
        do i=1,k
                do j=i+k,n,k
                        rvsg(i)=rvsg(i)+(log(p(j))-log(p(j-k)))**2.0d0
                end do
        end do
        print *,"rvsg=",rvsg

        nb=(n-k+1)/dble(k)
        rv=(sum(rvsg)/dble(k)-dble(nb)/dble(n)*rva)/(1.0d0-dble(nb)/dble(n))
        std=8.0d0/c**2.0d0*veps**2.0d0+c*4.0d0/3.0d0*rq
        std=sqrt(std)/(1.0d0-dble(nb)/dble(n))
        deallocate(rvsg)
        print *,"rv_all=",rva
        print *, "rv_zma=",rv
        print *, "sd=",std
        end subroutine odrvzma

        !check alerts per day to make statistics about the validicity of such the
        !strategy
        subroutine check_alert_record
        implicit none
        integer :: ios
        character(len=8) :: ad !alert date
        character(len=6) :: at !alert time
        character(len=4) :: stkcode
        real*8 :: ap    !alert price
        real*8 :: cp    !close price
        character :: actn       !action
        type(stkquote) :: cq !current quote
        
        !open alert_record
        open(unit=10,file=trim(alert_dir)//"alert_record.dat",status="old",iostat=ios)
        print *, "alert_record.dat iso=",ios
        !if alert_record for today exists 
        if( ios .eq. 0) then 
                open(unit=11,file=trim(alert_dir)//"alert_log.dat",position="append",iostat=ios)
                print *,"alert_log.dat ios=",ios
                do 
                        read(10,fmt=alert_record_fmt,iostat=ios) ad,at,stkcode,cp,actn
                        if(ios .eq. 0) then
                                call get_quote(stkcode,cq)
                                write(11,fmt=alert_record_log_fmt) &
                                ad,at,stkcode,cp,actn,cq%price,100*(cq%price-ap)/ap
                        else
                                exit
                        end if
                end do
                close(10)
                close(11)
        end if

        end subroutine check_alert_record

        subroutine track_quote
        implicit none
                character(len=4) :: stkcode
                integer :: unitnumber

                type(stkquote) :: lq !last quote
                type(stkquote) :: cq !current quote

                integer :: nstk,i
                character(len=4),dimension(:),allocatable:: stklist
                character(len=10) :: ctime
                real*8 :: t
                integer :: ios
                 
                open(unit=100,file=trim(intraday_dir)//"tracklist.dat",status="old",iostat=ios)
                call count_no_lines(trim(intraday_dir)//"tracklist.dat",nstk)
                allocate(stklist(nstk))
                read(100,*) stklist
                close(100)
                do
                        do i=1,nstk
                                stkcode=stklist(i)
                                call get_quote(stkcode,cq)
                                unitnumber=100+i
                                open(unit=unitnumber,file=trim(intraday_dir)//stkcode//"hk.qhf",position="append",iostat=ios)
                                backspace(unitnumber,iostat=ios)
                                if(ios.eq.0) then
                                        read(unitnumber,fmt=quote_fmt,iostat=ios) lq%date,&
                                        lq%time,lq%julian,lq%price,lq%volume
                                        if(ios.ne.0) then
                                                backspace(unitnumber)
                                                read(unitnumber,*,iostat=ios) lq
                                        end if
                                        print *,"last quote is "
                                        write(*,fmt=quote_fmt) lq%date,lq%time,lq%julian,lq%price,lq%volume
                                        if (cq%date//cq%time .ne. lq%date//lq%time) then
                                                write(unitnumber,fmt=quote_fmt) &
                                                cq%date,cq%time,cq%julian,cq%price,cq%volume
                                                if(cq%price .ne. lq%price) call check_price(stkcode,cq%price)
                                       end if
                                else
                                       write(unitnumber,fmt=quote_fmt) &
                                       cq%date,cq%time,cq%julian,cq%price,cq%volume
                                end if
                                close(unitnumber)
                        end do
                        call date_and_time(time=ctime)
                        read(ctime,*) t
                        if(t> 163000.0d0) exit
                        call sleep(120)
                end do 

                deallocate(stklist)
        end subroutine track_quote



        subroutine check_price(stkcode,price)
        implicit none
        character(len=4),intent(in) :: stkcode
        real(kind=rknd),intent(in) :: price
        
        character(len=8) :: cdate
        character(len=12) :: ctime
        character(len=4) :: stkc
        real(kind=rknd) :: qtl(9) 
        character(len=25) :: as !alert subject ,e.g.: 2000-00-00: buy 0175
        character(len=15) :: afn !alert file name ,e.g.: yyyymmdd_b_0175
        integer :: irecord,ios
        character :: actn

        call system("cd "//code_dir)
        open(unit=10,file=trim(daily_dir)//stkcode//"_frcst.dat")
        read(10,fmt='(A4,1X,I5)') stkc,irecord
        read(10,*) qtl
        close(10)
        
        call date_and_time(cdate,ctime)

        if(price .le. qtl(2) .or. price .gt.qtl(8)) then
                if(price .le. qtl(2)) then
                        !issue letter to buy
                        actn="B"
                else
                        actn="S"
                end if
                write(afn,fmt='(A8,A1,A4)') cdate,actn,stkcode
                open(unit=10,file=trim(alert_dir)//afn,status="new",iostat=ios)
                if(ios.eq.0) then !No alert has been issued.
                        as=afn
                        !write(as,fmt='(A8,A1,A4)') cdate,actn,stkcode
                        write(10,fmt='(A,F9.4)') "Current price:",price
                        if(actn.eq. "B") then
                                write(10,fmt='(A,F9.4)') "0.10 quantile:",qtl(2)
                                write(10,fmt='(A,F9.4)') "price with 1% return:",&
                                price*1.01
                                write(10,fmt='(A,F9.4)') "price with 2% return:",&
                                price*1.02
                        else
                                write(10,fmt='(A,F9.4)') "Current price:",price
                                write(10,fmt='(A,F9.4)') "0.90 quantile:",qtl(8)
                        endif
                        write(10,fmt='(A,F9.4)') "Estimated close price:",qtl(5)
                        write(10,fmt='(A,F9.4)') "Number of observations:",irecord
                        write(10,fmt='(A,/,5F9.4,/,F9.4,/,5F9.4)') "Quantiles: ",qtl
                        close(10)
                        call system('echo "Latest record in data file:" >> '//trim(alert_dir)//afn) 
                        call system("head -1 "//trim(daily_dir)//stkcode//"hk.dat >>"//trim(alert_dir)//afn)
                        call system("tail -6 "//trim(daily_dir)//stkcode//"_frcst.dat >>"//trim(alert_dir)//afn)
                        call system("mail -s "//trim(as)//" wan9c9@gmail.com < "//trim(alert_dir)//afn)
                        !keep a record of the alert in the file alert_record.dat 
                        open(unit=11,file=trim(alert_dir)//"alert_record.dat",position="append")
                        write(11,fmt=alert_record_fmt) cdate,ctime(1:6),stkcode,price,actn
                        close(11)
                        !call system("cat "//trim(alert_dir)//afn)
                end if
        end if
        end subroutine check_price

        !subroutine write_quote(cq, unitnumber)
        !implicit none
        !        type(stkquote) :: cq
        !        integer :: unitnumber
        !        write(unitnumber,fmt='(A14,1X,D16.8,1X,D16.8,I15)') &
        !        cq%date//cq%time,cq%julian,cq%price,cq%volume
        !end subroutine write_quote

        subroutine get_quote(stkcode,quote)
        implicit none

                character(len=4),intent(in) :: stkcode
                type(stkquote),intent(out) :: quote
                
                character(len=200) :: get_price
                character(len=10) :: datee,timee
                integer :: ios
                write(get_price, fmt='(A,A,A)') &
                                'wget -q -O tmp.txt &
                                 &"http://hk.rd.yahoo.com/finance/quotes&
                                &/internal/summary/download/*http://hk.finance.yahoo.com/d/quotes.csv?s=',&
                                stkcode,'.HK&f=d1t1l1v&e=.csv"'
                call system(get_price)
                call system("cat tmp.txt")

                open(unit=12,file="tmp.txt",status="old",iostat=ios)
                if(ios .eq.0) then
                        read(12,*) datee,timee,quote%price,quote%volume
                        print *, datee,timee,quote%price,quote%volume
                        call standardize_datestring(datee,quote%date)
                        call standardize_timestring(timee,quote%time)
                        quote%julian=datetime2julian(quote%date,quote%time)
                        !print *, quote
                        close(12)
                        call system("rm tmp.txt")
                else
                        write(*,*) "Cannot download quote."
                end if
        end subroutine get_quote
        
        subroutine format_quote(uff,ff)
        implicit none

                character(len=*),intent(in) :: uff
                character(len=*),intent(in) :: ff
                
                type(stkquote) :: q
                
                character(len=10) :: datee,timee
                integer :: ios 

                open(unit=10,file=uff,status="old")
                open(unit=11,file=ff)
               
                do  
                        read(10,*,iostat=ios) q%price,datee,timee,q%volume
                        if(ios .eq. 0) then
                                call standardize_datestring(datee,q%date)
                                call standardize_timestring(timee,q%time)
                                q%julian=datetime2julian(q%date,q%time)
                                write(11,fmt=quote_fmt) &
                                       q%date,q%time,q%julian,q%price,q%volume
                        else
                                exit
                        end if
                end do
                close(10)
                close(11)
        end subroutine format_quote

        subroutine standardize_datestring(datestring,sdate)
        implicit none
                character(len=*),intent(in) :: datestring
                character(len=len(datestring)) :: copy
                character(len=8),intent(out) :: sdate
               
                
                integer :: i,n,j,m,d,y
                n=len(datestring)

                copy=datestring
                j=1
                do i=1,n
                        if (datestring(i:i) .eq. "/") then
                                copy(i:i)=" "  
                                j=j+1
                                if(j>2) exit
                        endif
                enddo
                read(copy,*) m,d,y 
                write(sdate,fmt='(I4.4,I2.2,I2.2)') y,m,d
        end subroutine standardize_datestring

        subroutine standardize_timestring(timestring,stime)
        implicit none
                character(len=*),intent(in) :: timestring
                character(len=6),intent(out) :: stime
               
                integer :: i,n,h,m,pm
                character(len=len(trim(timestring))) :: copy
                
                n=len(trim(timestring))
                
                if(timestring(n-1:n) .eq."pm") then
                        pm=1
                else
                        if(timestring(n-1:n).eq."am") then
                                pm=0
                        else
                                Print *, "Error in time"
                        end if 
                endif
                
                copy=timestring
                copy(n-1:n)="  " 
                do i=1,n
                        if (timestring(i:i) .eq.":") then
                                copy(i:i)=" "  
                        endif
                enddo
                read(copy,*) h,m
                write(stime,fmt='(I2.2,I2.2,I2.2)') h+12*pm,m,0
        end subroutine standardize_timestring

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!! check if there exists 0 in data

        subroutine check_zero (x,cnt)
        implicit none
        type(daily),intent(in) :: x(:)
        integer,intent(out) :: cnt(size(x))
 
        real*8 :: tmp(size(x))

        
        tmp=x%ope*x%hig*x%low*x%clo*x%vol*x%adj
        where(tmp.eq.0) 
                cnt=1
        elsewhere
                cnt=0
        end where

        end subroutine check_zero
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!count number of lines/record in a specific file
!        subroutine count_no_lines(f,no_lines)
!        implicit none
!                character(len=*) :: f
!                integer :: no_lines
!
!                !write(sys_cmd,fmt='("wc ",A," > tmp_wc_file")') trim(f)
!                !call system(sys_cmd)
!                call system("wc "//trim(f)//" > tmp_wc_file")
!                open(unit=10,file="tmp_wc_file",status="old")
!                read (10,*) no_lines
!                close(10)
!                call system("rm -f tmp_wc_file")
!        end subroutine count_no_lines

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!       calculate base date in julian date
        
        subroutine cal_base_date( x, bd)
        implicit none
        type(daily),intent(in) :: x(:)
        integer,intent(out) :: bd

        
        type(daily) :: y(size(x))
        integer :: order,i,n
        integer :: ind_0(size(x))

        n=size(x)
        order=check_order(x%julian)
        if(order.eq. 1) then
                y=x
        else
                if(order.eq.-1) then
                        y=x(n:1:-1)
                else
                        print *, "error in data."
                end if
        end if
        
        call check_zero(y, ind_0)

        if( sum(ind_0).eq.0) then
                bd=y(1)%julian
        else
                do i=1,n
                        if(ind_0(i).eq.0) then
                                bd=y(i)%julian
                                return
                        end if
                end do
!                do i=n,1,-1
!                        if(ind_0(i).eq. 1) then
!                                bd=y(i)%julian
!                                j=max(1,i-2)
!                                if( any(ind_0(1:i-1).eq.0)) then
!                                        write(*, *) "possible trading &
!                        &exists before the BASE DATE."
!                                end if
!!                                return
!                        end if
!                end do
        end if

        end subroutine cal_base_date
!



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!   check_order checks the order of the array x   
!!!!   if x is non-decreasing, i=1                    
!!!!   if x is non-increasing, i=-1                     
!!!!   if x is not monotone,   i=0
!!!! 

        function check_order_d(x)
        implicit none
        real*8,intent(in) :: x(:)
        integer :: check_order_d

        integer :: n
        
        n=size(x)

        if ( all(x(2:n)-x(1:n-1).ge.0.0d0)) then
                check_order_d=1
        else
                if ( all(x(2:n)-x(1:n-1) .le.0.0d0)) then
                       check_order_d=-1
                else
                        check_order_d=0
                end if
        end if
        end function check_order_d

        function check_order_i(x)
        implicit none
        integer, intent(in) :: x(:)
        integer :: check_order_i

        integer :: n
        
        n=size(x)

        if ( all(x(2:n)-x(1:n-1).ge. 0)) then
                check_order_i=1
        else
                if ( all(x(2:n)-x(1:n-1) .le. 0)) then
                        check_order_i=-1
                else
                        check_order_i=0
                end if
        end if
        end function check_order_i
!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


        subroutine details2logincr(x,y)
        !this subroutine try to find the log-increment of x, stored in y
        implicit none
        type(daily), intent(in) :: x(:)
        type(daily), intent(out) :: y(size(x)-1)

        integer :: n
        type(daily) :: t(size(x))

        n=size(x)

        if(check_order(x%julian).eq.1) then
                t=x
        else
                t=x(n:1:-1)
        end if

        y%julian=t(2:n)%julian
        y%ope=log(t(2:n)%ope)-log(t(1:n-1)%adj)
        y%hig=log(t(2:n)%hig)-log(t(1:n-1)%adj)
        y%low=log(t(2:n)%low)-log(t(1:n-1)%adj)
        y%clo=log(t(2:n)%clo)-log(t(1:n-1)%adj)
        y%vol=log(t(2:n)%vol)-log(t(1:n-1)%vol)
        y%adj=log(t(2:n)%adj)-log(t(1:n-1)%adj)

        return

        end subroutine details2logincr

!******************************
!f=filename of data stored.
!x=output variable for stockdata
!nrec=output for number of record
!nrec0 = number of records with 0 in values of prices or volumes
!nl= number of lines of records read, optional 
!nl= nrec+nrec0 ,optional 

!******************************
!f=filename(including directory) of data stored.
!x=output variable for stockdata
!nrec=output for number of record
subroutine datafmttd_input(f,x,nrec)
        implicit none
        character(len=*), intent(in) :: f
        type(daily), dimension(:), intent(out) :: x
        integer, intent(out) :: nrec
        integer :: i,n
        integer :: stat

        n=size(X,1)
        open(unit=1,file=trim(f),status='old')
                
        do i=1,n   
                read(1,fmt='(I4,2I2.2,6D14.6)',iostat=stat) x(i)
                !x(i)%dat,x(i)%ope,x(i)%hig,x(i)%low,x(i)%clo,x(i)%vol,x(i)%adj
                if(stat .lt. 0) then
                        print *, "end_of_file occurs."
                        nrec=i-1
                        write(*,fmt='("Total number of records:",I5)') nrec
                        close(1)
                        return
                end if
                if(stat .gt. 0) then
                        print *, "An error occurs during reading."
                        print *, "The last record is the",i-1, "th"
                        print *, "The last record is:"
                        write( *,fmt='(I10,6E12.5)') x(i-1)
                        nrec=i-1
                        close(1)
                        return
                end if
        end do
        nrec=n
        close(1)
        !write(*,*) "Data_input finished." 
end subroutine datafmttd_input

subroutine data_output(x, flnm,ifmt)
implicit none
        type(daily), dimension(:), intent(in) :: x
        character(len=*), intent(in) :: flnm
        integer :: ifmt
        !ifmt=0, sequence input
        !ifmt=1, fmtted input
        !integer :: ndct 

        integer ::i,n,iost
        write(*,*) "******Enter into data_output.*********" 
        write(*,'(A,A)') "The file name is:",flnm
        n=size(x)
        write(*,'(A,I0)') "Number of record:",n
        open(unit=100, file=trim(flnm), status="replace",iostat=iost)
        if( iost .ne. 0) then
                write(*,*) "File with name", flnm," cannot be openned."
                write(*,*) "Data_output failed."
                return
        end if
        do i=1,n
                if(ifmt.eq.0) then
                        write(100,*) x(i)
                end if
                if(ifmt.eq.1) then
                        write(100, fmt='(I4.4,1X,I2.2,1X,I2.2,6D14.6)') x(i)
                end if
        end do
        close(100)
        write(*,*) "*********Data_output finished.*********"
end subroutine data_output

subroutine data_append(x, flnm,ifmt)
implicit none
        type(daily), dimension(:), intent(in) :: x
        character(len=*), intent(in) :: flnm
        integer :: ifmt
        !ifmt=0, sequence input
        !ifmt=1, fmtted input
        !integer :: ndct 

        integer ::i,n,iost
        write(*,*) "Enter into data_output."
        open(unit=100, file=trim(flnm), status="old",position="append",iostat=iost)
        if( iost .ne. 0) then
                write(*,*) "File with name", flnm," cannot be openned."
                write(*,*) "Data_output failed."
                return
        end if
        n=size(x)
        do i=1,n
                if(ifmt.eq.0) then
                        write(100,*) x(i)
                end if
                if(ifmt.eq.1) then
                        write(100, fmt='(I4,2I2.2,6D14.6)') x(i)
                end if
        end do
        close(100)
        write(*,*) "DATA APPEND finished."
end subroutine data_append

subroutine updatedata(stockcode)
implicit none
        integer,intent(in) :: stockcode

        type(daily), dimension(8000) :: x_t
        integer :: nrec 
        character(len=20) :: flnm="./data/0001hk.dat"
        character(len=24) :: flnm_tmp="./data/0001hk.dat_tmp"
 
        write(flnm_tmp,fmt='(A,I4.4,A)') "./data/",stockcode,"hk.dat_tmp"
        call data_input(flnm_tmp, x_t,nrec)
        
        if(nrec.gt.0) then
                write(flnm,fmt='(A,I4.4,A)') "./data/",stockcode,"hk.dat"
                call data_append(x_t(nrec:1:-1),flnm,1)
        else
                write(*,*) "Number of record = 0, nothing to do"
        end if
        
        return
end subroutine updatedata

subroutine fmthd(flnm)
implicit none
        character(len=*),intent(in) :: flnm

        type(daily), dimension(8000) :: x_t
        integer :: nrec 
        
        call data_input(flnm, x_t,nrec)
        
        if(nrec.gt.0) then
                write(*,'(A,I0)') "Number of record:",nrec
                write(*,*) "The first record is ", x_t(1)
                write(*,*) "The last record is ", x_t(nrec)
                if(check_order(x_t(1:nrec)%julian) .eq. 1) then
                	call data_output(x_t(1:nrec),flnm,1)
                else
                        call data_output(x_t(nrec:1:-1),flnm,1)
                end if
        else
                write(*,*) "Number of record = 0, nothing to do"
        end if
        return
end subroutine fmthd

!*********
!mean 
real*8 function mean(X)
implicit none
        REAL*8,dimension(:) :: x
        mean=sum(x)/size(x)
end function mean

!variance
real*8 function var(x)
implicit none
	real*8::x(:)
        var=sum((x-mean(x))**2)/(size(x)-1)
end function var

!covariance
real*8 function cov(x,y)
implicit none
	real*8 :: x(:),y(:)
        if(size(x).ne.size(y)) then 
                print *, "size is not consistent in COV"
        end if
	cov=sum((x-mean(x))*(y-mean(y)))/(size(x)-1)
end function cov

!correlation
real*8 function corr(x,y)
implicit none
	real*8::x(:),y(:)
        if(size(x).ne.size(y)) then 
                print *, "size is not consistent in COV"
        end if
	corr=cov(x,y)/sqrt(var(x)*var(y))
end function corr

!mean vector
function meanv(X)
implicit none
REAL*8,dimension(:,:) :: x
real*8 :: meanv(size(x,2))
	meanv=sum(x,dim=1)/size(x,1)
end function meanv

!covariance matrix
function covm(x)
implicit none
	real*8 :: x(:,:)
	real*8 :: covm(size(x,2),size(x,2))
	
	integer :: i,j,n
        n=size(x,2)
	do i=1,n
                do j=i,n
			covm(i,j)=cov(x(:,i),x(:,j))
			covm(j,i)=covm(i,j)
		end do
	end do
	
end function covm

!correlation matrix
function corrm(x)
implicit none
	real*8 :: x(:,:)
	real*8 :: corrm(size(x,2),size(x,2))
        
        integer :: i,j,n

        corrm=covm(x)

        n=size(x,2)
        do i=1,n
                do j=i+1,n
                        corrm(i,j)=corrm(i,j)/sqrt(corrm(i,i)*corrm(j,j))
                        corrm(j,i)=corrm(i,j)
                end do
        end do
end function corrm

!auto-correlation function
subroutine acf(x,mlag,r)
implicit none
        real*8,intent(in) :: x(:)
        integer, intent(in) :: mlag
        real*8, intent(out) :: r(mlag)
        
        integer :: i,n
        real*8 :: v
        
        v=var(x)
        n=size(x)
        do i=1,mlag
                r(i)=cov(x(1:n-i),x(i+1:n))/v
        end do

end subroutine acf
!cross auto-correlation function
subroutine xacf(x,mlag,r)
implicit none
        real*8,intent(in) :: x(:,:)
        integer, intent(in) :: mlag
        real*8, intent(out) :: r(0:mlag,size(x,2),size(x,2))
        
        integer :: i,j,k,n
       
        n=size(x,1) 
        do i=0,mlag
        do j=1,size(x,2)
        do k=1,size(x,2)           
                r(i,j,k)=corr(x(1+i:n,j),x(1:n-i,k)) 
        end do
        end do
        end do

end subroutine xacf
!cross auto-covariance
subroutine xacv(x,mlag,r)
implicit none
        real*8,intent(in) :: x(:,:)
        integer, intent(in) :: mlag
        real*8, intent(out) :: r(0:mlag,size(x,2),size(x,2))
        
        integer :: i,j,k,n
       
        n=size(x,1) 
        do i=0,mlag
        do j=1,size(x,2)
        do k=1,size(x,2)           
                r(i,j,k)=cov(x(1+i:n,j),x(1:n-i,k)) 
        end do
        end do
        end do
end subroutine xacv

real*8 function sd(x)
implicit none
	real*8::x(:)
	sd=sqrt(var(x))
end function sd

!log difference
function ldiff(x)
implicit none
        real*8 :: x(:)
        real*8 :: ldiff(size(x)-1)
                
        ldiff=log(x(2:size(x)))-log(x(1:size(x)-1))
end function ldiff

!log difference column-wisely
function ldiffc(x)
implicit none
        real*8 :: x(:,:)
        real*8 :: ldiffc(size(x,1)-1,size(x,2))
                
        ldiffc=log(x(2:size(x),:))-log(x(1:size(x)-1,:))
end function ldiffc

subroutine summary_sign_change(x, p2p,p2n,n2n,n2p)
implicit none
real*8, intent(in) :: x(:)
integer, intent(out) :: p2p, p2n, n2n, n2p

integer :: i
p2p=0
p2n=0
n2n=0
n2p=0
write(*,*) "number of observation:",size(x)
do i=2,size(x)
        if((x(i) > 0.0d0) .and. (x(i-1) > 0.0d0)) p2p=p2p+1 
        if((x(i) > 0.0d0) .and. (x(i-1) < 0.0d0)) n2p=n2p+1 
        if((x(i) < 0.0d0) .and. (x(i-1) > 0.0d0)) p2n=p2n+1 
        if((x(i) < 0.0d0) .and. (x(i-1) < 0.0d0)) n2n=n2n+1 
end do

write(*,*) "Summary of sign change:"
write(*,*) "p2p=",p2p
write(*,*) "p2n=",p2n
write(*,*) "n2p=",n2p
write(*,*) "n2n=",n2n
end subroutine summary_sign_change

!
! Calendar module
!
! Conversions to and from Julian dates and find day_of_the_week.
!
! Author: Robert Iles, March 1994
!
! There is no warranty on this code ........
!
!--------------------------------------------------------------
!--------------------------------------------------------------
!
! Return the day of the week 
!
!   Input
!     julian :: Integer, Optional
!               If present the julian day for which the weekday is
!               required, if absent the day of the week TODAY is
!               returned
!   Output
!     weekday :: Integer, optional
!               The day of the week, 0 = sunday
!     day :: Character*(*), optional
!               The name of the day of the week, e.g. 'Sunday'
!               Minimum length = 9
!     ierr :: Integer
!               Error return, 0=correct
!                            -1=invalid Julian day
!                            -2=neither day nor weekday specified
!
      subroutine day_of_week(julian, weekday, day, ierr)

      implicit none
      integer,intent(in),optional  :: julian
      integer,intent(out),optional :: weekday
      integer, intent(out)         :: ierr
      integer                      :: iweekday
      character*(*),intent(out),optional:: day
 
      ierr = 0
      if(present(julian)) then   
        if(julian < 0) then
          ierr = -1
          return
        endif
        iweekday = mod(julian+1, 7)
      else
        iweekday = date_to_julian(ierr=ierr)
        if(ierr/=0) return
      endif

      if(.not.present(day).and. .not.present(weekday)) then
        ierr=-2
        return
      endif

      if(present(day)) then
        if(iweekday.eq.0) then
          day = 'Sunday'
        else if(iweekday==1) then
          day = 'Monday'
        else if(iweekday==2) then
          day = 'Tuesday'
        else if(iweekday==3) then
          day = 'Wednesday'
        else if(iweekday==4) then
          day = 'Thursday'
        else if(iweekday==5) then
          day = 'Friday'
        else
          day = 'Saturday'
        endif
      endif
      if(present(weekday)) weekday=iweekday

      end subroutine day_of_week


        function datetime2julian(date,time) result(julian)
        !!!!!!!!!!!!!!!!!!
        !date="YYYYMMDD"
        !time="HHMMSS"
        !!!!!!!!!!!!!!!!!
        implicit none
                integer,parameter :: knd=kind(1.0d0)
                character(len=8),intent(in):: date
                character(len=6),optional,intent(in) :: time
                real(knd) :: julian,t
                !real :: datetime2julian

                integer :: day, month, year, hh,mm,ss,ierr

                read(date,fmt='(I4,I2,I2)') year,month,day

                ierr=0
                julian=date_to_julian(day,month,year,ierr=ierr)
                if(ierr.ne.0) then
                        print *, "error in transforming date & time to julian"
                        print *, "Error code: ",ierr
                        julian=0.0d0
                else
                        if(present(time)) then
                                read(time,fmt='(I2,I2,I2)') hh,mm,ss
                                t=ss/60.0d0
                                t=(t+mm)/60.0d0
                                t=(t+hh)/24.0d0
                        else
                                t=0.0d0
                        endif
                        julian=julian+t
                endif
        end function datetime2julian

!-------------------------------------------------------------------------
!
! Convert a Julian day to a day/month/year
!
!   Input
!     julian :: Integer
!               The julian day to convert
!   Output
!     day   :: Integer, optional
!               The day of the month
!     month :: Integer, optional
!               The day of the month
!     year  :: Integer, optional
!               The day of the month
!     values :: Integer(:),optional
!               Array, minimum size=3, values(1)=year
!                                      values(2)=month
!                                      values(3)=day
!     ierr :: Integer
!               Error return, 0=correct
!                            -1=invalid Julian day
!                            -2=no output specified
!

function julian2date_f(julian)
implicit none
integer :: julian
integer :: values(3),ierr
character(len=10) :: julian2date_f

	call julian_to_date(julian, values=values,ierr=ierr)
	write(julian2date_f, fmt='(I4.4,"-",I2.2,"-",I2.2)') values

end function 	

      subroutine julian_to_date(julian,day,month,year,values,ierr)

      implicit none
      integer,parameter      :: igreg=2299161
      integer,parameter      :: k = kind(0.0d0)
      integer,intent(in)     :: julian
      integer,intent(out),optional    :: day, month, year
      integer,intent(out),optional    :: values(:)
      integer,intent(out)             :: ierr
      integer                :: ia, ja, jb, jc, jd, je, iday, imonth, iyear
      real(kind=k)           :: xc

      if(julian < 0) then
        ierr = -1
        return
      else
        ierr = 0
      endif
      if(.not.present(values).and..not.present(day).and. &
         .not.present(month).and..not.present(year)) then
        ierr=-2
        return
      endif

      if (julian.ge.igreg) then
        ia = (real(julian-1867216,k)-0.25d0)/36524.25d0
        ja = julian + 1+ia-int(0.25*ia)
      else
        ja = julian
      end if

      jb = ja + 1524
      xc = (real(jb-2439870,k)-122.1d0)/365.25d0
      jc = 6680.0d0 + xc
      jd = 365*jc + int(0.25d0*real(jc,k))
      je = int(real(jb-jd,k)/30.6001d0)

      iday = jb - jd - int(30.6001d0*real(je,k))

      imonth = je - 1
      if (imonth.gt.12) imonth = imonth - 12

      iyear = jc - 4715
      if (imonth.gt.2) iyear = iyear - 1
      if (iyear.le.0) iyear = iyear - 1
!
! Assign output values
!
      if(present(values)) then
        values(1) = iyear
        values(2) = imonth
        values(3) = iday
      endif
      if(present(year)) year=iyear
      if(present(month)) month=imonth
      if(present(day)) day=iday

      end subroutine julian_to_date
!-----------------------------------------------------------
!
! Convert a day/month/year to a Julian day 
!
!  If VALUES and one of DAY, MONTH and YEAR are missing this
!  will return the Julian day for TODAY
!
!   Input
!     day   :: Integer, optional
!               The day of the month
!     month :: Integer, optional
!               The day of the month
!     year  :: Integer, optional
!               The day of the month
!     values :: Integer(:),optional
!               Array, minimum size=3, values(1)=year
!                                      values(2)=month
!                                      values(3)=day
!   Output
!     julian :: Integer, Function Result
!               The julian day to convert
!     ierr :: Integer
!               Error return, 0=correct
!                            -1=invalid year
!                            -2=invalid month
!                            -3=invalid day
!                            -4=invalid date (29th Feb, non leap-year)
!
      function date_to_julian(day, month, year, values, ierr) result(julian)
!
!  gregorian started midday oct 15 1582
!
      implicit none
      integer,parameter :: igreg=15+31* (10+(12*1582))
      integer,parameter      :: k = kind(0.0d0)
      integer,parameter :: limit(12) = &
                           (/31,29,31,30,31,30,31,31,30,31,30,31/)
      real(kind=kind(0.0d0))::xi
      integer,intent(in),optional :: day, month, year
      integer,intent(in),optional :: values(:)
      integer:: julian, v(9)
   
      integer :: iy, ja, jm, jy, ierr
      integer :: iday,imonth,iyear

      if(present(values)) then
        iyear = values(1)
        imonth = values(2)
        iday = values(3)      
      else if(present(day).and.present(month).and.present(year)) then
        iyear = year ; imonth = month ; iday = day
      else
        call date_and_time(values=v)
        iyear = v(1)
        imonth = v(2)
        iday = v(3)      
      endif

      if (iyear==0) then
        ierr=-1
        return
      endif
      if (imonth <= 0 .or. imonth > 12) then
        ierr=-2
        return
      endif
      if (iday > limit(imonth)) then
        ierr=-3
        return
      endif
      if(imonth==2.and.iday==29.and.iyear>1582)then
        if(mod(iyear,4)/=0) then
          ierr=-4
          return
        endif
      endif
      iy = iyear
      if (iyear.lt.0) iy = iy + 1
      if (imonth.gt.2) then
        jy = iy
        jm = imonth + 1
      else
        jy = iy - 1
        jm = imonth + 13
      end if
!
!  note: SLIGHTLY STRANGE CONSTRUCTIONS REQUIRED TO AVOID PROBLEMS WITH
!        OPTIMISATION OR GENERAL ERRORS UNDER VMS!
!
        julian = iday + int(30.6001d0*real(jm,k))
        if (jy.lt.0) then
          xi = 365.25d0*real(jy,k)
          if(int(xi).ne.xi)xi = xi -1
          julian = julian + int(xi)
        else
          julian = julian + 365*jy
          julian = julian + int(0.25*real(jy,k))
        end if
        julian = julian + 1720995

        if (iday+ (31* (imonth + (12*iy))) >= igreg) then
          ja=jy/100
          julian = julian - ja
          julian = julian + 2
          julian = julian + ja/4
        end if

      end function date_to_julian












!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


        subroutine augment_details(x,y)
        !this subroutine try to find the log-return of x, stored in y
        implicit none
        type(daily), intent(in) :: x(:)
	real*8 :: y(size(x),13)
         
        integer :: n
        type(daily) :: t(size(x))

        n=size(x)


        !need to consider the order to data in time ascending	
        if(check_order(x%julian).eq.1) then
                t=x
        else
                t=x(n:1:-1)
        end if

	
        y(:,1)=t%julian
	y(:,2)=t%ope
	y(:,3)=t%hig
	y(:,4)=t%low
	y(:,5)=t%clo
	y(:,6)=t%vol
	y(:,7)=t%adj
	y(1,8)=log(y(1,2))
	y(2:n,8)=log(y(2:n,2))-log(y(1:n-1,7))
	y(:,9)=log(y(:,3))-log(y(:,2))
	y(:,10)=log(y(:,4))-log(y(:,2))
	y(:,11)=log(y(:,5))-log(y(:,2))
	y(:,12)=log(y(:,7))-log(y(:,2))
	y(1,13)=log(y(1,6))
	y(2:n,13)=log(y(2:n,6))-log(y(1:n-1,6))
	
        return

        end subroutine augment_details


subroutine t2wr(x,eps,y) !Transform data to the Whole Real line
        !this subroutine try to transform the calculated stock price data
        !the original x store the tranformed daily detail as
        ! julian_date, OPN, UPR, DNR(<0),RCL,VOL
        ! which is calculate in the subroutine augment_details.
        !in the output, y
        ! they are transformed as
        ! julian_date, OPN, log(UPR+eps),log(-DNR+eps), arctan((RCL-1/2)/(pi/2+eps)), VOL
        ! so that except, the date, all the values have the possibility to varies on the whole real line
        ! although the transformation is not uniform, and may cause un-preceived difficulties.
        implicit none
        real*8,intent(in) :: x(:,:)
        real*8, intent(in) :: eps
	real*8, intent(out) :: y(size(x,1),size(x,2))
     
        y(:,1)=x(:,1)
        y(:,2)=x(:,2)
        y(:,3)=log(x(:,3)*1.0d4+eps)
        y(:,4)=log(-x(:,4)*1.0d4+eps)
        y(:,5)=tan((x(:,5)-.5d0)*(3.14159265358/2.0d0+eps))
        y(:,6)=x(:,6)
 
        return
        
        end subroutine t2wr
        
! this subroutine estimates the 
!

subroutine data_input(f,x,nrec,nrec0,nl)
        implicit none
        character(len=*), intent(in) :: f
        type(daily), dimension(:), intent(out) :: x
        integer, intent(out) :: nrec
        integer, intent(out), optional :: nrec0
        integer, intent(out), optional :: nl

        character(len=60) :: header(8)
        integer :: i,n
        real*8 :: prod
        integer :: dd,mm,yy,ierr, stat, tmp
        character :: ans,warn
        character(len=9) :: day
        


        n=size(X,1)
        open(unit=1,file=trim(f),status='old')
        read(1,*,iostat=stat) tmp
        if(stat.lt.0) then
                print *,"error in reading the first row"
                print *,"set number of record=0"
                nrec=0
                return
        end if
        if(stat.eq.0) then
                write(*,*) "The first row is data."
                rewind(1)
        else
                !stat.gt.0
                rewind(1)
                read(1,*,iostat=stat) header       
                if( stat .eq. 0) then
                        !write(*,*) "Headers are found."
                                if( trim(header(1)).ne."Date") then
                                        print *,trim(header(1)), ".ne. Date"
                                        print *,"please check file."
                                        nrec=0
                                        return
                                end if
                                if( trim(header(2)).ne."Open") then
                                        print *,trim(header(2)), ".ne. Open"
                                        print *,"please check file."
                                        nrec=0
                                        return
                                end if
                                if( trim(header(3)).ne."High") then
                                        print *,trim(header(3)), ".ne. High"
                                        print *,"please check file."
                                        nrec=0
                                        return
                                end if
                                if( trim(header(4)).ne."Low") then
                                        print *,trim(header(4)), ".ne. Low"
                                        print *,"please check file."
                                        nrec=0
                                        return
                                end if
                                if( trim(header(5)).ne."Close") then 
                                        print *,trim(header(5)), ".ne. Close"
                                        print *,"please check file."
                                        nrec=0
                                        return
                                end if
                                if( trim(header(6)).ne."Volume") then 
                                        print *,trim(header(6)), ".ne. Volume"
                                        print *,"please check file."
                                        nrec=0
                                        return
                                end if
                                if( trim(header(7)).ne."Adj") then 
                                        print *,trim(header(7)), ".ne. Adj"
                                        print *,"please check file."
                                        nrec=0
                                        return
                                end if
                                if( trim(header(8)).ne."Close") then 
                                        print *,trim(header(8)), ".ne. Close"
                                        print *,"please check file."
                                        nrec=0
                                        return
                                end if
                end if
        end if
       
        i=0 
        nl=0
        nrec0=0
        ans='n'
        warn='y'
        do  
                i=i+1
                if( i.gt. n) then
                        print *, "Maximal number of record read,&
                         & but not to the end of the data file."
                        nrec=i-1
                        write(*,fmt='("Total number of records:",I5)') nrec
                        close(1)
                        return
                end if
 
                read(1,fmt=*,iostat=stat) &
                yy,mm,dd,x(i)%ope,x(i)%hig,x(i)%low,x(i)%clo,x(i)%vol,x(i)%adj
                
                x(i)%julian=date_to_julian(day=dd,month=mm,year=yy,ierr=ierr)
                prod=x(i)%ope*x(i)%hig*x(i)%low*x(i)%clo*x(i)%vol*x(i)%adj
                               !x(i)%dat,x(i)%ope,x(i)%hig,x(i)%low,x(i)%clo,x(i)%vol,x(i)%adj
                !x(i)%vol=real(m)
                !read(1,fmt=*,iostat=stat) x(i)
                if(stat .lt. 0) then
                        print *, "Read till the end of the file."
                        nrec=i-1
                        write(*,fmt='("Total number of records read:",I5)') nrec
                        close(1)
                        return
                end if
                if(stat .gt. 0) then
                        print *, "An error occurs during reading."
                        print *, "The last record is the",i-1
                        print *, "The last record is:"
                        write( *,fmt='(I10,6E12.5)') x(i-1)
                        nrec=i-1
                        close(1)
                        return
                end if
                if( prod .eq. 0.0d0 ) then
                        !check the abnormal data
                         nrec0=nrec0+1 
                         !price does not change
                         if(x(i)%ope.eq.x(i)%hig .and. x(i)%hig.eq.x(i)%low &
                         .and. x(i)%low.eq.x(i)%clo .and. x(i)%clo.eq.x(i)%adj .and. x(i)%vol.eq.0.0d0) then
                         	i=i-1
                         else 
                                 !price changes, but volume is zero
                         	if(ans.ne.'Y') then
                                        if( ans .eq. 'N') then
                                                i=i-1
                                        else
                                        if( warn.eq.'y') then
                                                write(*, *) "Data with zero trade, occured at date", julian2date_f(x(i)%julian)
                                                call day_of_week(julian=x(i)%julian, day=day, ierr=ierr)
                                                write(*, fmt='(I4.4, I3.2, I3.2, 1X, A,I9,1X,8F15.8)') yy,mm,dd, day,x(i)
                                        end if
                                        write(*,*) "Do you want to include this data?"
                                        write(*,*) " y=include for this time, "
                                        write(*,*) " Y=include all such data, "
                                        write(*,*) " n=no for this time,"
                                        write(*,*) " N=no for all such data."
                                        read(*,*) ans
                                        if (ans.eq. 'n') i=i-1
                                        end if
                		end if
			end if
                end if
                nl=nl+1
        end do
        nrec=i
        close(1)
        !write(*,*) "Data_input finished." 
end subroutine data_input



	!synchronization two stock data with respect to the first one
   	!size(z,1)=size(x)
      
      	subroutine syn_2stock(x,y,z)
        implicit none

        type(daily), intent(in) :: x(:)
        type(daily), intent(in) :: y(:)
        type(daily), intent(out) :: z(size(x),2)

        integer :: i,nx,ny,nz,j

        nx=size(x)
        ny=size(y)
        nz=nx
        
	!check if there zero-volume in either stock
        !if(any(x(:)%vol.eq.0.0d0).or.any(y(:)%vol .eq. 0.0d0)) then
        !       write(*,*) "Zero trading in orignal series. Check data first."
        !       z(:,:)%vol=0.0d0
        !else
                !check if order is increasing
                if( check_order(x(:)%julian).ne.1 .or.&
                check_order(x(:)%julian).ne.1) then
                  write(*, *) "The date is not in increasing order."
                end if
                !endif
                if( x(nx)%julian .ge. y(ny)%julian) then
                        print *, "The second data is not up to date with the first one &
                                & day."
                end if

                j=ny
                z(:,1)=x
                if(all(y(ny-nz+1:ny)%julian.eq.x(:)%julian)) then
                	z(:,2)=x(ny-nz+1:ny)
                else
                	do i=nz,1,-1
                	!write(*,*) i
                		do
                			if (y(j)%julian.gt.z(i,1)%julian) then
                				j=j-1
                			else
                				exit
                			end if
                		end do
                		if(y(j)%julian.eq.z(i,1)%julian) then
                			z(i,2)=y(j)
                			j=j-1
                		else
                			write(*,*) "data missing for 2nd data at date", z(i,1)%julian
                			z(i,2)%vol=0.0d0
                		end if
                		if( j.eq. 0) then
                			write(*,*) "Error, Exit."
                			return
                		end if
                      	end do
                end if
        !end if
      end subroutine syn_2stock


!ip= the garch order in the volatility
!iq= the arch order in the volattility
!theta= the vector of parameters
!se = standard errors
!sc= scores

subroutine elfmts(nobs,nvar,r,nf,flm)
use gnuplot_int
        implicit none
        integer :: nobs,nvar
        real*8 :: r(nobs,nvar),flm(nvar,nvar)
        real*8 :: f(nobs,nf)

        real*8 ::v(nvar,nvar),w(nvar)
        integer :: nf
        integer,parameter :: mlag=5
        real*8 :: xcrr(0:mlag,nvar,nvar),a(nvar,nvar),rw(nvar-1)
        integer :: lwork
        real*8 :: work(3*nvar-1)
        integer :: info,i
        real*8 :: minr

        call xacv(r,mlag,xcrr)
        a=0.0d0
        do i=1,mlag
         a=a+matmul(xcrr(i,:,:),transpose(xcrr(i,:,:)))
        end do
        v=a
        lwork=3*nvar-1
        call dsyev('V','U',nvar,v,nvar,w,work,lwork,info)
        rw=w(1:nvar-1)/w(2:nvar)
        write(*, *) "Eigenvalues,    consecutive ratio"
        write(*, '(2F13.3)') w(1)
        write(*, '(2F13.3)') (w(i),rw(i-1),i=2,nvar)

        call gplot1(w)
        !call gphisto(w)
        minr=rw(nvar-1)
        nf=nvar-1
        do i=nvar-2,1,-1
        if ( rw(i) .lt. minr) then
                minr=rw(i)
                nf=i
        end if
        end do
        nf=nvar-nf
        flm=0.0d0
        flm(:,1:nf)=v(:,nvar:nvar-nf+1:-1)
        f=matmul(r,flm(:,1:nf))
        

        print *, "number of factors:", nf
end subroutine elfmts

function outer_product(x,y)
	implicit none
	real*8 :: x(:),y(:)
	real*8 :: outer_product(size(x),size(y))

	integer :: i

	do i=1,size(x)
		outer_product(i,:)=x(i)*y
	end do
	end function outer_product

end module nvst_int
