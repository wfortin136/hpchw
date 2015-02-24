program main
  implicit none

  integer ::STATUS=1, write_tstep=5000
  character(20) :: prog_name, infile, outfile, m_count
  real*8, dimension(:,:,:), allocatable :: conc
  integer nt, n, l
  integer m, i, j, leng
  logical left, right, bottom, top
  real*8 t, delta, delta_t, time1, time2
  real*8 u, v, courant_stab, temp1, temp2
 
  call getarg(0, prog_name)
  call getarg(1, infile)

  open(9, file=trim(infile))

  !!Set intitial conditions Concentration at time 0
  read(9,*)
  read(9,*) l, n, t, nt, u
  v=u*0.57
  close(9)
  print *, "LxL=",l,"x",l
  print *, "NxN=",n,"x",n
  print*, "T=",t
  print*, "NT=",nt
  print*, "u=",u
  print*, "v=",v

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
  
  do m=2, nt
    !if(mod(m,write_tstep) == 0) then
    !  write(m_count,*) m
    !  outfile = trim(prog_name(3:))//"_"//trim(adjustl(m_count))//".out"
    !  open(m, file=outfile)
    !end if
    
    do j=1, n
          if(j==1) then
              conc(1,:,j-1)=conc(1,:,n)
          end if

          if(j==n) then
              conc(1,:,j+1)=conc(1,:,1)
          end if

          do i=1, n
            if(i==1) then
              conc(1,i-1,:)=conc(1,n,:)
            end if

            if(i==n) then
              conc(1,i+1,:)=conc(1,1,:)
            end if 
            !print *, conc(1,i,j)
            
            temp1= conc(1,i+1,j)+conc(1,i-1,j)+conc(1,i,j+1)+conc(1,i,j-1)
            temp2=u*(conc(1,i+1,j)-conc(1,i-1,j))+v*(conc(1,i,j+1)-conc(1,i,j-1))
            !print *, temp1, temp2
            conc(2,i,j) = .25*temp1-((delta_t/(2*delta))*temp2)
          end do
        end do
        !if(mod(m,write_tstep)==0) then
        !  do i=1, n
        !    write(m,*) conc(2,i,:)
        !  end do
        !  close(m)
        !end if
        conc(1,:,:) = conc(2,:,:)
  end do
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

end program main
