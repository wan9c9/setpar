
      subroutine drift(lambda,r,para,k,v)
      implicit none
      real*8,intent(in) :: lambda
      integer,intent(in) :: r
      real*8,intent(in) :: para(6)
      real*8,intent(in) :: k
      real*8,intent(out) :: v
      
      integer :: i

      real*8 :: prob
      
      prob=exp(-lambda)
      v=1.0d0+(para(1)+para(2)*lambda)**k*prob
      do i=1,1000
        prob=prob*lambda/i
        if( i.le.r) then
                 v=v+(para(1)+para(2)*lambda+para(3)*i)*prob
        else
                 v=v+(para(4)+para(5)*lambda+para(6)*i)*prob
        end if
      end do
        
      end subroutine drift


      program test_drift
      implicit none
      real*8 :: lambda
      integer :: r=5
      real*8 :: para(6)
      real*8 :: k=2.0d0
      real*8 :: v

      integer :: i

      para(1:3)=(/0.5d0,0.4d0,0.3d0/)
      para(4:6)=(/0.5d0,0.7d0,0.4d0/)


      do i=1,100
      lambda=i*0.1d0
      call drift(lambda,r,para,k,v)
      if( v .lt. 1.0d0+lambda**k) then
              write(*,'(10F15.5)') lambda,v,1.0d0+lambda**k
      end if

      end do

      end program





