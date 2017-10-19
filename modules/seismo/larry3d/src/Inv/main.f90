program main
real, dimension(:), allocatable :: x, t, tt
integer, dimension(:),allocatable :: i, mpoin
!integer :: j
print*,"dimension of right-hand-side vector?"
read*,m
print*,"dimension of sparse matrix arrays?"
read*,n
allocate(x(n),i(n),t(m),tt(m),mpoin(0:m))
call lsqrinv(x,i,n,t,tt,mpoin,m)
end program

