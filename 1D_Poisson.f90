! Written: DANG Truong
! Purpose: Solving Poisson equation in 1D.
! Reference: johansen1998.pdf

program johansen1998_1D
implicit none

interface

  real(8) function exact(x)
  implicit none
  real(8),intent(in)::x
  end function exact

  real(8) function second_derivative(x)
  implicit none
  real(8),intent(in)::x
  end function second_derivative


end interface

real(8),parameter :: lambda = 0.5d0
real(8),parameter :: l = 1.0d0

integer(4),parameter:: N = 80
integer(4):: i
integer(4):: i1,i2
real(8)::deltax
real(8),dimension(N+1)::x ! mesh points
real(8),dimension(N+1)::xcp ! control points

real(8),dimension(N,N)::A,A_tempo,A_inv
real(8),dimension(N)::F,u

real(8),dimension(N+2)::u_app,u_exact

real(8)::d2,d1

real(8):: L2norm
logical:: exist

deltax = l/(N-1+lambda)
!print*,'deltax =',deltax

do i=1,N
  x(i) = (i-1)*deltax
end do
x(N+1) = l

d1 = x(N+1)-x(N)
d2 = deltax+d1

!print*,'x = ',x

xcp(1)=0.0d0

do i=1,N
  if (i==1) then
    xcp(i+1)=xcp(i) + deltax/2.0d0
  else
    xcp(i+1)=xcp(i) + deltax
  end if
end do

!print*,'xcp = ',xcp

A = 0.0d0

do i1=1,N
  do i2=1,N
    if(i1.ne.1 .and. i1.ne.N) then
      if(i1==i2) then
        A(i1,i2) = -2.0d0/deltax**2
        A(i1,i2-1) = 1.0d0/deltax**2
        A(i1,i2+1) = 1.0d0/deltax**2
      end if
    end if
    if (i1==1) then
       A(i1,i1) = -3.0d0/deltax**2
       A(i1,i1+1) = 1.0d0/(3.0d0*deltax**2)
    end if
    if (i1==N) then
       A(i1,i1) = -1.0d0/(deltax**2 * lambda)
       A(i1,i1-1) = 1.0d0/(deltax*lambda)* &
                    (1.0d0/deltax -d2/(d1*(d2-d1)))
       A(i1,i1-2) = 1.0d0/(deltax*lambda)* &
                    d1/(d2*(d2-d1))
    end if
  end do
end do

A_tempo = A
call matinv(A_tempo,N,A_inv)

u = 0.0d0
F=0.0d0
do i1=1,N
   F(i1) = second_derivative(xcp(i1+1))
end do


u = matmul(A_inv,F)
!print*,'u =',u


u_app = 0.0d0
u_exact = 0.0d0

do i1=1,N
   u_app(i1+1) = u(i1)
end do


u_exact(1) = exact(0.0d0)

do i1=1,N+1
   u_exact(i1) = exact(xcp(i1))
end do


L2norm = 0.0d0
do i1= 1,N+1
  L2norm = L2norm + (u_exact(i)-u_app(i))**2*deltax
end do

print*,'L2 norm =',L2norm
!open(unit=3,file='L2norm')
!write(3,*) N, L2norm
!close(3)

inquire(file="L2norm", exist=exist)
  if (exist) then
    open(3, file="L2norm", status="old", position="append", action="write")
  else
    open(3, file="L2norm", status="new", action="write")
  end if
  write(3, *) log(N*1.0d0), log(L2norm)
  close(3)


open(unit=1,file='u_app')
do i1=1,N+1
  write(1,*) xcp(i1),u_app(i1)
end do
close(1)

open(unit=2,file='u_exact')
do i1=1,N+1
  write(2,*) xcp(i1),u_exact(i1)
end do
close(2)


end program

real(8) function exact(x)
implicit none
real(8),intent(in)::x

exact = x**2*(x-1.0d0)
end function exact

real(8) function second_derivative(x)
implicit none
real(8),intent(in)::x

second_derivative = 6.0d0*x - 2.0d0
end function second_derivative

