program main
  implicit none

  integer ::STATUS=1
  character(20) :: prog_name, infile, outfile
  real*8, dimension(:,:,:), allocatable :: conc
  integer nt, n, l
  integer m, i, j, leng
  logical left, right, bottom, top
  real*8 t, delta, delta_t
  real*8 u, v, courant_stab, temp1, temp2
 
  call getarg(0, prog_name)
  call getarg(1, infile)

  open(9, file=trim(infile))
  
  outfile = trim(prog_name(3:))//".out"
  open(99, file=outfile)

  !l=1
  !n=400
  !t=1000000
  !nt=20000

  !!Set intitial conditions Concentration at time 0
  !u=0.0000005
  read(9,*)
  read(9,*) l, n, t, nt, u
  v=u*0.57
  close(9)
  print *, l,n,t,nt,u,v

  delta=l/(n*1.0)
  delta_t=t/nt
  courant_stab = delta/(sqrt(2.)*sqrt(u**2+v**2))

  if(delta_t>courant_stab) then
    write(*,*) "Courant Stability criteria not met"
    write(*,*) "time resolution is greater than discretized mesh time step"
    write(*,*) "delta t > (L/N)/[sqrt(2) velocity_magnitude]"
    write(*,*) delta_t, ">",  courant_stab
    call EXIT(STATUS)
  end if
  
  !!allocate 3D array. The first dimension holds 2 values, 
  !!current time and previous time
  !!For each time step, the x and y values are stored respectively
  !!Alter default parameters of Fortran array starting at 1, and 
  !!start it at 0 and then extend it 1 beyond the required dimensions
  !!This is in order to store a buffer of data for wrap around
  !!and neighbor data
  allocate(conc(2,0:n+1,0:n+1))
 
  call gaussian_dist(dble(1),conc(1,:,:))
  !print *, conc(1,:,:)
  print *, "Gaussian Complete"
  leng=size(conc(1,1,:))-2
  do m=2, nt
  !  conc(m,:,:,:,:)=conc(1,:,:,:,:)
  !  do k=1, delta
  !    do kk=1, delta
        do i=1, n
          if(i==1) then
            conc(1,i-1,:)=conc(1,n,:)
          end if

          if(i==n) then
            conc(1,i+1,:)=conc(1,1,:)
          end if

          do j=1, n
            if(j==1) then
              conc(1,:,j-1)=conc(1,:,n)
            end if

            if(j==n) then
              conc(1,:,j+1)=conc(1,:,1)
            end if 
            !print *, conc(1,i,j)
            
            temp1= conc(1,i+1,j)+conc(1,i-1,j)+conc(1,i,j+1)+conc(1,i,j-1)
            temp2=u*(conc(1,i+1,j)-conc(1,i-1,j))+v*(conc(1,i,j+1)-conc(1,i,j-1))
            !print *, temp1, temp2
            conc(2,i,j) = .25*temp1-((delta_t/(2*delta))*temp2)
          end do
          if(m==10000) then
            write(99,*) conc(2,i,:)
          end if
        end do
        !print *, "****", m, "*****"
        !print *, conc(2,1:leng,1:leng)
        conc(1,:,:) = conc(2,:,:)
  !    end do
  !  end do
  end do
  !print *, conc(2,1:leng,1:leng)
  !print*,leng
  !do m=1, nt
  !do k=1,delta
  !  do kk=1, delta
  !    print*,k, "-", kk
  !    do i=1,n
  !      print *, conc(4,i,1:leng)
  !    end do
  !    print*, " "
  !  end do
  !end do
  !end do
  
  deallocate(conc)
  close(99)
contains

subroutine gaussian_dist(ampl, y)
  implicit none
  real*8, intent(in) :: ampl
  real*8, dimension(:,:), intent(inout) :: y
  
  real*8 p, xnot, ynot, sigmax, sigmay
  real*8 xterm, yterm
  integer i, j
 
  sigmax=.25
  sigmay=.25
  !Get midpoint of global matrix
  xnot=size(y(1,:))/2
  ynot=size(y(:,1))/2
  
  !column-major access for efficiency 
  do j=1, size(y(1,:))
    do i=1, size(y(:,1))
      xterm = (i-xnot)**2/(2*sigmax**2)
      yterm = (j-ynot)**2/(2*sigmay**2)
      p=-(xterm+yterm)
      y(i,j)=ampl*exp(p)
    end do
  end do
end subroutine gaussian_dist

!used http://www.design.caltech.edu/erik/Misc/Gaussian.html
!http://en.wikipedia.org/wiki/Box%E2%80%93Muller_transform
!to define an algorithm that transformed the uniform distribution
!to a gausian distribution. Standard deviation of 1 is assumed
subroutine gaussian_dist_rand(mean, y)
  implicit none
  real*8, intent(in) :: mean
  real*8, dimension(:,:), intent(inout) :: y
  
  real*8 x1, x2, w
  integer i, j
  !call random_seed()
  !print*,size(y)
  do i=1, size(y(1,:))
    do j=1, size(y(:,1)),2
      !print*,j 
      call box_muller(x1,x2, w)
      y(i,j)=x1*w+mean
      if(j+1 <= size(y(:,1))) then
        y(i,j+1)=x2*w+mean
      end if
    end do
  end do
end subroutine gaussian_dist_rand

subroutine box_muller(z0,z1,s)
  implicit none
  real*8 r1, r2, z0,z1,s 
  do
    call random_number(r1)
    call random_number(r2)
    z0 = (2. *r1) - 1.
    z1 = (2. *r2) - 1.
    s = z0**2 + z1**2
    if(s<1.0) then
      exit
    end if
  end do
  !log is natural log, not log base 10
  s=sqrt( (-2. * log(s))/s)

end subroutine box_muller

subroutine set_c_ij(x)
  real*8, dimension(:,:), intent(inout) :: x
  
  real*8, dimension(:,:), allocatable :: a
  real*8 mean
  !print*, size(x)
  allocate(a(size(x(1,:)),size(x(:,1))))
  mean=x(1,1)
  call gaussian_dist_rand(mean,a)
  x=a
  
end subroutine set_c_ij

end program main
