program main
  use mpi
  implicit none

  integer ::STATUS=1, write_tstep=5000
  integer stat(MPI_STATUS_SIZE), ierr, rank, numrank,r
  character(30) :: prog_name, infile, outfile, m_count 
  real*8, dimension(:,:,:), allocatable :: conc
  real*8, dimension(:), allocatable :: leftbuf, rightbuf, upbuf, downbuf
  real*8, dimension(:,:), allocatable :: recvbuf, sendbuf
  integer dims(2), coords(2), borders(4), cartcom
  integer reqs(8), stats(MPI_STATUS_SIZE,8)
  logical periods(2), reorder
  integer entry_count
  integer nt, n, l, te, te2
  integer m, x, i, j, icount, jcount, leng
  integer l_iwidth, istart, iend, l_jwidth, jstart, jend
  integer istart2, iend2, jstart2, jend2
  real*8 t, delta, delta_t, time1, time2
  real*8 u, v, courant_stab, temp1, temp2

  call MPI_INIT(ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD, numrank, ierr)
  
  dims(:) = sqrt(numrank*1.)
  !print * , dims(:)
  periods(:) = .TRUE.
  reorder=.FALSE.
  call MPI_CART_CREATE(MPI_COMM_WORLD, 2,dims, periods, reorder, cartcom, ierr)
  call MPI_COMM_RANK(cartcom, rank, ierr)
  call MPI_CART_COORDS(cartcom, rank, 2, coords, ierr)
  !1:up, 2:down, 3:left, 4:right
  call MPI_CART_SHIFT(cartcom, 0,1, borders(1), borders(2), ierr)
  call MPI_CART_SHIFT(cartcom, 1,1, borders(3), borders(4), ierr)
  if(rank ==0) then 
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
  end if
  !broadcast parameters to all ranks
  call MPI_BCAST(l, 1, MPI_INTEGER, 0, cartcom, ierr)
  call MPI_BCAST(n, 1, MPI_INTEGER, 0, cartcom, ierr)
  call MPI_BCAST(t, 1, MPI_REAL8, 0, cartcom, ierr)
  call MPI_BCAST(nt, 1, MPI_INTEGER, 0, cartcom, ierr)
  call MPI_BCAST(u, 1, MPI_REAL8, 0, cartcom, ierr)
  call MPI_BCAST(v, 1, MPI_REAL8, 0, cartcom, ierr)
  
  !!allocate 3D array. The first dimension holds 2 values, 
  !!current time and previous time
  !!For each time step, the x and y values are stored respectively
  !!Alter default parameters of Fortran array starting at 1, and 
  !!start it at 0 and then extend it 1 beyond the required dimensions
  !!This is in order to store a buffer of data for wrap around
  !!and neighbor data
  allocate(conc(2,0:n+1,0:n+1))
  entry_count=(n+2)*(n+2)
  
  delta=l/(n*1.0)
  delta_t=t/nt
  !do check and initialization once
  if(rank==0) then
    courant_stab = delta/(sqrt(2.)*sqrt(u**2+v**2))

    !print *, delta, delta_t, courant_stab
  
    if(delta_t>courant_stab) then
      write(*,*) "Courant Stability criteria not met"
      write(*,*) "time resolution is greater than discretized mesh time step"
      write(*,*) "delta t > (L/N)/[sqrt(2) velocity_magnitude]"
      write(*,*) delta_t, ">",  courant_stab
      !call MPI_FINALIZE(ierr)
      call EXIT(STATUS)
    end if
  
    call gaussian_dist(dble(1),conc(1,:,:))
    !set buffers for initial timestep to wrap around conditions
    conc(1,0,:)=conc(1,n,:)
    conc(1,n+1,:)=conc(1,1,:)
    conc(1,:,0)=conc(1,:,n)
    conc(1,:,n+1)=conc(1,:,1)
    print *, "Gaussian Complete"
  end if 
  !broadcast initialization matrix to all ranks
  call MPI_BCAST(conc(1,:,:), entry_count, MPI_REAL8, 0, cartcom, ierr)

  !Set parameters for each rank in loop
  l_iwidth=n/dims(1)
  istart=coords(2)*l_iwidth+1
  iend = istart+l_iwidth-1
  l_jwidth=n/dims(2)
  jstart=coords(1)*l_jwidth+1
  jend = jstart+l_jwidth-1

  !Must allocate in row major order in order to keep data contiguous for 
  !Irecv and Isend
  !http://www.mpi-forum.org/docs/mpi-2.2/mpi22-report/node343.htm
  !Problems Due to Data Copying and Sequence Association
  allocate(recvbuf(l_jwidth,4))
  allocate(sendbuf(l_jwidth,4))
  !print *, rank, istart, iend, jstart, jend

  time1=MPI_Wtime()
  do m=2, nt
    !if(mod(m,write_tstep) == 0 .AND. rank==0 ) then
    !  write(m_count,*) m
    !  outfile = trim(prog_name(3:))//"_"//trim(adjustl(m_count))//".out"
    !  open(m, file=outfile)
    !end if
    !if(rank ==0) then
    !  print*, m
    !end if
        do i=istart, iend

          if(i==istart .AND. m/=2) then
            !Test for receive from left
            !print *, rank, reqs(3)
            call MPI_WAIT(reqs(3),stat, ierr)
            conc(1,istart-1, jstart:jend)=recvbuf(:,3)
        
            !print *, rank, reqs(3)
          end if

          if(i==iend .AND. m/=2) then
            !Test for receive from right
            call MPI_WAIT(reqs(4),stat, ierr)
            conc(1,iend+1, jstart:jend)=recvbuf(:,4) 
          end if

          do j=jstart, jend
            if(j==jstart .AND. m/=2) then
              !Test for receive from up
              call MPI_WAIT(reqs(1),stat, ierr)
              conc(1,istart:iend, jstart-1)=recvbuf(:,1)
            end if

            if(j==jend .AND. m/=2) then
              !Test for receive from down
              call MPI_WAIT(reqs(2),stat, ierr)
              conc(1,istart:iend, jend+1)=recvbuf(:,2)
            end if 
            !print *, conc(1,i,j)
            
            temp1= conc(1,i+1,j)+conc(1,i-1,j)+conc(1,i,j+1)+conc(1,i,j-1)
            temp2=u*(conc(1,i+1,j)-conc(1,i-1,j))+v*(conc(1,i,j+1)-conc(1,i,j-1))
            !print *, rank, temp1, temp2
            conc(2,i,j) = .25*temp1-((delta_t/(2.*delta))*temp2)
            !print *, rank, conc(2,i,j)
          end do
        end do
        !print *, rank, "send calls"
        !if(m/=2) then
        !  call MPI_WAITALL(8, reqs, stats , ierr)
        !  print*, rank
          !call MPI_WAIT(reqs(5),stat, ierr)
          !call MPI_WAIT(reqs(6),stat, ierr)
          !call MPI_WAIT(reqs(7),stat, ierr)
          !call MPI_WAIT(reqs(8),stat, ierr)
        !end if
        !FROM LEFT
        call MPI_IRECV(recvbuf(:,3),l_jwidth, & 
              MPI_REAL8, borders(3), 4, cartcom, reqs(3), ierr)
        !FROM RIGHT
        call MPI_IRECV(recvbuf(:,4),l_jwidth, & 
              MPI_REAL8, borders(4), 3, cartcom, reqs(4), ierr)
        !FROM UP
        call MPI_IRECV(recvbuf(:,1),l_iwidth, & 
              MPI_REAL8, borders(1), 2, cartcom, reqs(1), ierr)
        !FROM DOWN
        call MPI_IRECV(recvbuf(:,2),l_iwidth, & 
              MPI_REAL8, borders(2), 1, cartcom, reqs(2), ierr)
        
        !TO LEFT
        sendbuf(:,3) = conc(2,istart, jstart:jend)
        call MPI_ISEND(sendbuf(:,3),l_jwidth, & 
              MPI_REAL8, borders(3), 3, cartcom, reqs(4+3), ierr)
        !TO RIGHT
        sendbuf(:,4) = conc(2,iend, jstart:jend)
        call MPI_ISEND(sendbuf(:,4),l_jwidth, & 
              MPI_REAL8, borders(4), 4, cartcom, reqs(4+4), ierr)
        !TO UP
        sendbuf(:,1) = conc(2,istart:iend, jstart)
        call MPI_ISEND(sendbuf(:,1),l_iwidth, & 
              MPI_REAL8, borders(1), 1, cartcom, reqs(4+1), ierr)
        !TO DOWN
        sendbuf(:,2) = conc(2,istart:iend, jend)
        call MPI_ISEND(sendbuf(:,2),l_iwidth, & 
              MPI_REAL8, borders(2), 2, cartcom, reqs(4+2), ierr)

        !if(mod(m,write_tstep)==0) then
        !  if(rank/=0) then
        !    call MPI_SEND(conc(2,istart:iend, jstart:jend),l_iwidth*l_jwidth, & 
        !        MPI_REAL8, 0, rank, cartcom, ierr) 
        !  else
        !    do r=1, numrank-1
        !      call MPI_CART_COORDS(cartcom, r, 2, coords, ierr)
        !      istart2=coords(2)*l_iwidth+1
        !      iend2 = istart2+l_iwidth-1
        !      jstart2=coords(1)*l_jwidth+1
        !      jend2 = jstart2+l_jwidth-1
        !      call MPI_RECV(conc(2,istart2:iend2, jstart2:jend2),l_iwidth*l_jwidth, &
        !        MPI_REAL8, r, r, cartcom, stat, ierr)
        !    end do
        !    do x=1, n
        !      write(m,*) conc(2,x,:)
        !    end do
        !  end if
        !close(m)
        !end if
        !set current time step to previous and increment forward in loop
        
        conc(1,istart:iend,jstart:jend) = conc(2,istart:iend,jstart:jend)
  end do
  time2=MPI_Wtime()
  print*, time2-time1
  deallocate(conc)
  deallocate(recvbuf)
  CLOSE(99)
  call MPI_FINALIZE(ierr)
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
