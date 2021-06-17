subroutine z_root2(output, func, x1, x2, xacc, iprint, ierr)

!real function z_root2(func,x1,x2,xacc,iprint,ierr)

!-----------------------------------------------------------------------------!
!
!   +-------+    ALKYON Hydraulic Consultancy & Research
!   |       |
!   |   +---+
!   |   | +---+
!   +---+ |   |
!         +---+
!
!  0. Update history
!
!     Version Date       Modification
!
!     0.01    29/11/1999 Initial version
!     0.02    07/11/1999 Test added to check boundaries, and reverse if necessary
!                        Bug fixed in assigning answer
!     0.03    02/09/2002 Maximum number of iterations set to 20, instead of 10
!
!  1. Purpose
!
!     Find zero crossing point of function FUNC between the
!     initial values on either side of zero crossing
!
!  2. Method
!
!     Ridders method of root finding
!
!     adapted from routine zridddr
!     Numerical Recipes
!     The art if scientific computing, second edition, 1992
!     W.H. Press, S.A. Teukolsky, W.T. Vetterling and B.P. Flannery
!
!  3. Parameter list
!
!     Name    I/O  Type  Description
!
!     func     i    r    real function
!     x1       i    r    initial x-value on left/right side of zero-crossing
!     x2       i    r    initial x-value on right/left side of zero-crossing
!     xacc     i    r    accuracy, used as |x1(i)-x2(i)|< xacc
!     iprint   i    i    Output channel and test level
!     ierr     o    i    Error indicator
!
!  4. Subroutines used
!
!     Func      user supplied real function
!
!  5. Error messages
!
!     ierr = 0   No errors occured during iteration process
!            1   Iteration halted in dead end, this combination may NEVER occur
!            2   Maximum number of iterations exceeded
!            3   Solution jumped outside interval
!
!  6. Remarks
!
!     It is assumed that the x1- and x2-coordinate lie
!     on different sides of the actual zero crossing
!
!     The input parameter IPRINT is used to generate test output.
!     If IPRINT==0, no test output is created
!              > 0, test output is directed to the file connected to unit LUPRINT=IPRINT
!                   if no file is connected to this unit, no output is written
!
!
!implicit none
!
real, intent(out) :: output

real func                         ! external function
real, intent (in) :: x1           ! x-value at one side of interval
real, intent (in) :: x2           ! x-value at other side of interval
real, intent (in) :: xacc         ! requested accuracy
integer, intent (in) :: iprint    ! number of output channel, only used when
!!!!!!!!!integer, intent (out) :: ierr     ! error indicator
!
real unused                       ! default value
real zriddr                       ! intermediate function value
real xx1,xx2,xx                   ! local boundaries during iteration
integer maxit                     ! maximum number of iteration
integer luprint                   ! unit of test output
logical lopen                     ! check if a file is opened

parameter (maxit = 20)
external func
!
integer iter      ! counter for number of iterations
real fh           ! function value FUNC(xh)
real fl           ! function value FUNC(xl)
real fm           ! function value FUNC(xm)
real fnew         ! function value FUNC(xnew)
real s            ! temp. function value, used for inverse quadratic interpolation
real xh           ! upper (high) boundary of interval
real xl           ! lower boundary of interval
real xm           ! middle point of interval
real xnew         ! new estimate according to Ridders method
!
!!!!!!!!!ierr   = 0        ! set error level
unused =-1.11e30  ! set start value
!
xx1 = x1          ! copy boundaries of interval to local variables
xx2 = x2
!
luprint = iprint
!
if(luprint > 0) then
  inquire(unit=luprint,opened=lopen)
  if(.not.lopen) then
    luprint = 0
    write(*,'(a,i4)') 'Z_ROOT2: invalid unit number:',iprint
  end if
end if
!
! check boundaries on requirement x2 > x1
!
if(xx1 > xx2) then
  xx  = xx1
  xx1 = xx2
  xx2 = xx
end if
!
fl = func(xx1)
fh = func(xx2)
!
!if(luprint > 0) write(luprint,'(a,4e13.5)') &
!&  'Z_ROOT2: xx1 xx2 fl fh:',xx1,xx2,fl,fh
!
if((fl > 0. .and. fh < 0.) .or. (fl < 0. .and. fh > 0.))then
   xl = xx1
   xh = xx2
   zriddr = unused
!
   do iter=1,maxit
      xm = 0.5*(xl+xh)
      fm = func(xm)
      s = sqrt(fm**2-fl*fh)
      if(s == 0.) goto 9000
      xnew = xm+(xm-xl)*(sign(1.,fl-fh)*fm/s)
!
!      if(luprint>0) write(luprint,'(a,4e13.5)') &
!&       'Z_ROOT2: xm,fm,s,xnew:',xm,fm,s,xnew
!
      if (abs(xnew-zriddr) <= xacc) then
!        if(luprint>0) write(luprint,'(a)') 'Z_ROOT2: xnew=zriddr'
        goto 9000
      end if
!
      zriddr = xnew
      fnew = func(zriddr)
      if (fnew == 0.) goto 9000
!
      if(sign(fm,fnew) /= fm) then
         xl = xm
         fl = fm
         xh = zriddr
         fh = fnew
      elseif(sign(fl,fnew) /= fl) then
         xh = zriddr
         fh = fnew
      elseif(sign(fh,fnew) /= fh) then
         xl = zriddr
         fl = fnew
      else
         !!!!!!!!ierr = 1
         goto 9000
      endif
!
      if(abs(xh-xl) <= xacc) goto 9000
!
      if(luprint > 0) write(luprint,'(a,i4,5e14.6)') &
&     'Z_ROOT2: iter,x1,x2,|x1-x2|,xacc,z:', iter,xl,xh,abs(xl-xh),xacc,fnew
!
   end do
!!!!!!!   ierr = 2
   if(luprint > 0) write(luprint,'(a)') 'Z_ROOT2: -> ierr=2'
   goto 9000
else if (fl == 0.) then
  zriddr = xx1
else if (fh == 0.) then
  zriddr = xx2
else
!!!!!!!!!!!  ierr = 3
  goto 9999
! 'root must be bracketed in zriddr'
endif
!
9000 continue
!
output = zriddr
!z_root2 = zriddr
!
!!!IGNORE THEE ERROR :P
!!!if(luprint > 0) write(luprint,'(a,2i3,5e13.5)') &
!!!&     'Z_ROOT2: ierr,iter,xl,xh,acc,x0,z0:', ierr,iter,xl,xh,xacc,z_root2,func(z_root2)
!!!!
9999 continue
!
!!!return
end subroutine z_root2
