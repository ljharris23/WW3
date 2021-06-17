!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!--------------------------------!!!!!
!!-------------Q-XNL4v4-----------------!!
!!!!!--------------------------------!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine q_xnl4v4(aspec,sigma,angle,nsig,nang,depth,xnl,diag,ierr)
!------------------------------------------------------------------------------
!
!   +-------+    ALKYON Hydraulic Consultancy & Research
!   |       |    Gerbrant van Vledder
!   |   +---+
!   |   | +---+  Last update: 25 June 2003
!   +---+ |   |  Release: 5.0
!         +---+
!
! do not use m_xnldata
use m_constants
use serv_xnl4v5
implicit none
!------------------------------------------------------------------------------
!  0. Update history
!
!     08/01/2000 Initial version
!     12/01/2001 Updated interface
!     13/01/2001 Inclusion of diagonal term
!     14/02/2002 Upgrade to release 4.0, depth added to input
!     20/08/2002 quad depth adapted in the case of WAM-depth scaling
!                then deep water is assumed for conversion of A(sig,theta) -> N(kx,ky)
!                Search option for nearest grid included
!     23/08/2002 Allocation of work arrays set to fixed size
!     11/09/2002 Filtering of energy densities introduced and restructure
!     14/04/2003 Format of test write statement corrected
!     03/05/2003 Computation and output of triplets enabled
!     12/06/2003 Export spectral grid in case of Q_INTEG>1
!     16/06/2003 Switch IQ_SYM included
!                Allocation of dynamic data array's moved to Q_ALLOCATE
!     24/06/2003 Range of loop for IK3 made dependent on value of IQ_SYM
!     25/06/2003 Bug fixed in assigment of contribution of diagonal term
!
!  1. Purpose:
!
!     Compute nonlinear transfer for a given action density spectrum
!     on a given wave number and direction grid
!
!  2. Method
!
!     Compute nonlinear transfer in a surface gravity wave spectrum
!     due to resonant four wave-wave interactions
!
!     Methods: Webb/Resio/Tracy/VanVledder
!
!
!  3. Parameter list:
!
! Type    I/O          Name               Description
!------------------------------------------------------------------------------
integer,intent(in)  :: nsig             ! number of radian frequencies
integer,intent(in)  :: nang             ! number of directions
real,   intent(in)  :: aspec(nsig,nang) ! Action density spectrum as a function of (sigma,theta)
real,   intent(in)  :: sigma(nsig)      ! radian frequencies
real,   intent(in)  :: angle(nang)      ! directions in radians (sector or full circle)
real,   intent(in)  :: depth            ! water depth in m
real,   intent(out) :: xnl(nsig,nang)   ! nonlinear quadruplet interaction computed with
!                                         a certain exact method (k,theta)
real,   intent(out) :: diag(nsig,nang)  ! Diagonal term for WAM based implicit integration scheme
integer, intent(out) :: ierr          ! error indicator
!
!  4. Error messages
!
!  5. Called by:
!
!     XNL_MAIN
!
!  6. Subroutines used
!
!     Q_STACK --- check
!     Q_INIT --- check
!     Q_CTRGRID --- check
!     Q_T13V4 --- check
!     Q_SEARCHGRID --- check
!
!  7. Remarks
!
!     The external action density spectrum is given as N(sigma,dir)
!     The internal action density spectrum is given as N(kx,ky)
!
!     These 2 spectra are conected via the Jacobian transformation
!
!                cg
!     N(kx,ky) = -- N(sig,theta)
!                 k
!  8. Structure
!
!  9. Switches
!
! 10. Source code
!------------------------------------------------------------------------------
! local variables
!---------------------------------------------------------------------------------------
integer iaq      ! counter for directions
integer jaq      ! counter for directions
integer ikq      ! counter for wave numbers
integer iang     ! counter for directions
integer ia       ! counter for directions
integer ik       ! counter for wave numbers
integer idir1    ! direction in degrees of k1 (for integration test)
integer idir3    ! direction in degrees of k3 (for integration test)
real period      ! periodicity for direction, used in conversion of 2-spectra
real diagk1      ! diagonal term for k1
real diagk3      ! diagonal term for k3
!
real qn_max      ! maximum action density
real qn_min      ! minimum action density
!
real cg(nsig)         ! group velocity for conversion of spectrum and transfer
!
integer ia1,ia3,ja3   ! limits for directional loops
integer jk3           ! start of k3 loop
integer ik1,ik3       ! counters for wave number loop
integer nloc          ! number of points on locus
!
integer igrid         ! status of grid file
real t13              ! value of sub-integral
real k_rat            ! local ratio of wave numbers
real a_dif            ! directional difference
real jacobian         ! Jacobian
real qn1,qn3          ! action densities in k1 and k3
!
!  testing of diagonal term on a low level
!
real diagk1_0         ! saved value of diagk1
real diagk3_0         ! saved value of diagk3
real dq1              ! small change in action density of n1
real dq3              ! small change in action density of n3
real t13_0            ! Original estimated of diagonal term
real t13_1,t13_3      ! perturbed estimated of diagonal term
!
integer ifil_dir      ! indicator for filtering of directional criterion
integer ifil_krat     ! indicator for filtering of wave number ratio criterion
integer ifil_dens     ! indicator for filtering of action density criterion
integer ifil_tot      ! indicator for filtering due to any criterion
integer nfil_dir      ! counter to indicate filtering of directional criterion
integer nfil_krat     ! counter to indicate filtering of wave number criterion
integer nfil_dens     ! counter to indicate filtering of action density criterion
!
integer ntot_conf     ! total number of configurations
integer ntot_filt     ! total number of filtered configurations
!
!
!------------------------------------------------------------------------------
!!!!!!!!!!call q_stack('+q_xnl4v4')
!
! initialisations
!------------------------------------------------------------------------------
ierr = 0              ! error status
diag = 0                ! initialize output diagonal term
!
!
if(iq_type==3) then
  q_depth = depth         ! water depth to be used in computation
else
  q_depth = q_maxdepth
end if
!--------------------------------------------------------------------------
!  generate basic grid of loci and store loci in memory and to datafile
!--------------------------------------------------------------------------
if(iq_screen >= 1) write(iscreen,'(a)') 'Q_XNL4V4: Checking interaction grid '
!
if(iq_search==0 .or. iq_type/=3) then
  call q_init
  call q_ctrgrid(2,igrid)
  if(iq_err /= 0) goto 9999
!
  if(igrid/=0) then
    !!!!!!!!call q_error('e','NOGRID','No proper grid exists')
    print *, "Error 23568"
    goto 9999
  end if
!
  if(iq_make ==3) then
    !!!!!!!call q_error('e','MAKEGRID','Only computation of grid')
    print *, "Error 2352678"
    goto 9999
  end if
!------------------------------------------------------------------------------
!  set overall scale factor resulting from optional SEARCH for nearest grid
!------------------------------------------------------------------------------
!
  q_scale = 1.
!------------------------------------------------------------------------
else
!
!  search nearest valid grid and compute additional WAM scale factor
!  only active when IQ_SEARCH==1 .AND. IQ_TYPE==3
!
  call q_searchgrid(depth,igrid)


  if(igrid/=0) then
    !!!!!!!!!call q_error('e','NOGRID','No proper grid exists')
    print *, "Error 6782"
    goto 9999
  end if
!
  if(iq_err /=0) goto 9999
end if
!
    goto 9999
  end if
!
  if(iq_err /=0) goto 9999
end if
!
!------------------------------------------------------------------------------
!  convert input action density spectrum from A(sigma,theta) -> N(kx,ky)
!
do ikq=1,nkq
  call z_cmpcg(sigma(ikq),q_depth,q_grav,cg(ikq))
  do iaq=1,naq
    nspec(ikq,iaq) = aspec(ikq,iaq)/q_k(ikq)*cg(ikq)
  end do
end do
!
!------------------------------------------------------------------------------
!
!--------------------------------------------------------------------------------------
!  integration over all possible configurations
!--------------------------------------------------------------------------------------
xnl = 0.
qn_max = maxval(nspec)
!
!--------------------------------------------------------------------------------------
do ik1 = 1,nkq
  if(iq_screen >= 1) write(iscreen,'(a,2i4,e12.3)') 'Q_XNL4V4: k1 nk d:',ik1,nkq,q_depth
  jk3 = ik1
  if(iq_sym==0) jk3 = 1
!
  do ia1 = iaq1,iaq2          ! loop over selected part of grid, set in q_init
!
    qn1 = nspec(ik1,ia1)
!
    do ik3 = jk3,nkq           ! compute only half-plane
      do ia3 = 1,naq           ! loop over all possible wave directions
        qn3 = nspec(ik3,ia3)
!
        if(iq_screen>=3) write(iscreen,'(a,4i4)') 'Q_XNL4V4: ik1 ia1 ik3 ia3:',ik1,ia1,ik3,ia3
!
!  computes distances in wave number space
!
        a_dif = 180. - abs(180. - abs(q_ad(ia1) - q_ad(ia3)))
        k_rat = max(q_k(ik1)/q_k(ik3), q_k(ik3)/q_k(ik1))
        qn_min = qf_frac*qn_max/(q_k(ik3)/q_k(1))**7.5
        qn_min = qf_frac*qn_max*q_kpow(ik3)
!
        ifil_dir  = 0
        ifil_krat = 0
        ifil_dens = 0
        ifil_tot  = 0
!
!  perform filtering
!
!  directional difference
!
        if(a_dif > qf_dmax) then
          ifil_dir = 1
        end if
!
!  wave number ratio
!
        if(k_rat > qf_krat) then
          ifil_krat = 1
        end if
!
!  energy density filtering
!
        if(qn1 < qn_min .and. qn3 < qn_min) then
          ifil_dens = 1
        end if
!
!
        if(ifil_dir==0 .and. ifil_krat==0 .and. ifil_dens==0 .or. iq_filt==0) then
!?        if(a_dif < qf_dmax .and. k_rat < qf_krat .or. iq_filt==0) then
!
!  perform integration along locus
!
          call q_t13v4(ik1,ia1,ik3,ia3,t13,diagk1,diagk3)
!
!
          if(iq_err /= 0) goto 9999
!
!  check contribution T13 with the computed with triplet method
!
!!/R           qt13 = 0.
!!/R           do iqtr = 1,ktriplets
!!/R             qt13 = qt13 + w_qtr(iqtr)*nspec(i_qtr(iqtr,1),i_qtr(iqtr,2))* &
!!/R &                   nspec(i_qtr(iqtr,3),i_qtr(iqtr,4))*nspec(i_qtr(iqtr,5),i_qtr(iqtr,6))
!!/R           end do
!!/R           write(iscreen,*) 'CHECK T13 QT13:',t13,qt13
!
!
!  take care of additional scale factor aring from search of nearest grid
!
          t13    = t13*q_scale
          diagk1 = diagk1*q_scale
          diagk3 = diagk3*q_scale
!
!  take care of symmetric storing of interactions
!  and factor 2 due to symmetry (if activated)
!
          if(iq_sym==1) then
            t13 = 2.*t13
            diagk1 = 2.*diagk1
            diagk3 = 2.*diagk3
          end if
!
          ja3 = ia3
          if(iq_grid==1 .and. ia3 < iaref) ja3 = naq-ia1+1
          xnl(ik1,ia1)  = xnl(ik1,ia1) + t13*q_k(ik3)*q_delta*q_dk(ik3)
          if(iq_sym==1) xnl(ik3,ja3)  = xnl(ik3,ja3) - t13*q_k(ik1)*q_delta*q_dk(ik1)
!
!  add diagonal term
!
          diag(ik1,ia1) = diag(ik1,ia1) + diagk1*q_k(ik3)*q_delta*q_dk(ik3)
!!/F &       ik1,ia1,ik3,ia3,qn1,qn3,t13,a_dif,k_rat,&
!!/F &       ifil_dir,ifil_krat,ifil_dens,ifil_tot
      end do
    end do
  end do
end do
!
!
!
!
!
! write number of triplets that have been written
!
!
!------------------------------------------------------------------------------
! in the case of a symmetric sector, copy results to other part
!
! Examples: naq=5, iaref=3: 1,2,3,4,5 ->   Q(1)=Q(5)
!                                          Q(2)=Q(4)
!                                          Q(3)=Q(3)
!                                                     iaq+jaq=naq+1
!           naq=6, iaref=4: 1,2,3,4,5,6 -> Q(1)=Q(6)
!                                          Q(2)=Q(5)
!                                          Q(3)=Q(4)
!
if(iq_grid==1) then
  do ikq = 1,nkq
    do iaq=iaref,naq
      jaq = naq+1-iaq
      xnl(ikq,jaq) = xnl(ikq,iaq)
    end do
  end do
end if
!
!------------------------------------------------------------------------------
if(iq_screen>=2) write(iscreen,'(a)') 'Q_XNL4V4: Main computation ended'
!
!  Convert transfer from (kx,ky) grid to (sigma,theta) grid
!
do ikq=1,nkq
  jacobian = q_k(ikq)/cg(ikq)
  do iaq=1,naq
    xnl(ikq,iaq) = xnl(ikq,iaq)*jacobian
  end do
end do
!                                                   !
9999 continue
!
!!!!!!!!!!call q_stack('-q_xnl4v4')
!
return
end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!--------------------------------!!!!!
!!-------------Q-T13V4------------------!!
!!!!!--------------------------------!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine q_t13v4(ik1,ia1,ik3,ia3,t13,diagk1,diagk3)
!------------------------------------------------------------------------------
!
!   +-------+    ALKYON Hydraulic Consultancy & Research
!   |       |    Gerbrant van Vledder
!   |   +---+
!   |   | +---+  Last update: 5 September 2003
!   +---+ |   |  Release: 5.0
!         +---+
!
! do not use m_xnldata
use m_constants
implicit none
!
!  0. Update history
!
!     25/02/1999  Initial version
!     14/04/1999  Extra check in GET_LOC if locus exists in internal database
!     12/10/1999  Error handling improved
!     15/01/2001  Interface extended with diagonal term
!     06/05/2002  Criterion f34_mod added to computational procedure
!     14/08/2002  Integration simplified
!     22/08/2002  Integration modified depending on actual number of non-zero points
!     26/09/2002  Boundary check for sector grid activated
!     15/04/2003  Bug fixed in handling of periodicity
!                 Nearest bin integration enabled, including diagonal term
!     25/04/2003  Output to triplet arrays for nearest bin
!     03/05/2003  Output of triplets for bi-linear interpolation enabled
!     04/06/2003  Parameter IQ_INT renamed IQ_INTEG
!     13/06/2003  Test of integration for case of nearest bin interpolation
!     25/06/2003  Bug fixed in computation of partial derivatives for contribution to
!                 diagonal term
!     27/08/2003  Short-cut when number of non-zero points on locus is ZERO
!     05/09/2003  Switches for test output in nearest bin approach modified
!
!  1. Purpose:
!
!     Compute the function T13, defined as a line integral around a locus
!
!  2. Method
!
!     See Tracy and Resio (1982) and Van Vledder (1999)
!
!  3. Parameter list:
!
! Type    I/O          Name             Description
!------------------------------------------------------------------------------
integer, intent(in) ::  ik1     !    Index of k-component of wave number k1
integer, intent(in) ::  ia1     !    Index of a-component of wave number k1
integer, intent(in) ::  ik3     !    Index of k-component of wave number k3
integer, intent(in) ::  ia3     !    Index of a-component of wave number k3
real, intent(out)   ::  t13     !    Value of line integral over a specific locus
real, intent(out)   ::  diagk1  !    Contribution to diagonal term of k1
real, intent(out)   ::  diagk3  !    Contribution to diagonal term of k3
!
!  4. Error messages
!
!  5. Called by:
!
!     Q_XNL4V4
!
!  6. Subroutines used
!
!     Q_GETLOCUS
!     Q_PUT_BTRIPLETS
!     Q_PUT_NTRIPLETS
!
!  7. Remarks
!
!     The action density product term is given by:
!     P = n1.n2.(n3+n4)-(n1+n2).n3.n4
!
!     This term is rewritten as:
!
!     P = n1.n2.n3 + n1.n2.n4 - n1.n3.n4 - n2.n3.n4
!       = n1.n3.(n2-n4) + n2.n4.(n1-n3)

!
!  8. Structure
!
!  9. Switches
!
!     /S  enable subroutine tracing
!     /T  enable test output
!     /N  enable interpolation using nearest point
!
! 10. Source code:
!------------------------------------------------------------------------------
!     Local variables
!
integer iloc          ! counter along locus
integer ifnd           ! indicator if correct locus is found
integer ja2,ja2p      ! direction indices for interpolation of k2
integer jk2,jk2p      ! wave number indices for interpolation of k2
integer ja4,ja4p      ! direction indices for interpolation of k4
integer jk4,jk4p      ! wave number indices for interpolation of k4
integer ikq,iaq       ! counters
!
real sumt13           ! sum along locus
real qn1,qn2,qn3,qn4  ! action densities at wave numbers k1, k2, k3 and k4
real nprod            ! wave number product
real t2,t4            ! tail factors for k2 and k4
real qd1,qd3          ! contribution to diagonal term
real rterm            ! product term along locus
!
real qn13p            ! product of N1 and N3
real qn13d            ! difference of N1 and N3
!
!
!------------------------------------------------------------------------------
!!!!!!!!!!!!call q_stack('+q_t13v4')
!
t13    = 0.
diagk1 = 0.
diagk3 = 0.
!
!
if(ik1==ik3 .and. ia1==ia3) goto 9999  ! skip routine if k1=k3
!
!  obtain information requested locus based on a information
!  about a precomputed locus, as stored in the database file
!
call q_getlocus(ik1,ia1,ik3,ia3,ifnd)
!
if(ifnd==0 .or. nlocusx==0) then
  t13 = 0.
  goto 9999
end if
!---------------------------------------------------------------------------------------
qn1 = nspec(ik1,ia1)
qn3 = nspec(ik3,ia3)
!
qn13p = qn1*qn3      ! compute product
qn13d = qn3-qn1      ! compute difference
!
sumt13 = 0
!
!    3-----------4 ja2p         w1 = (1-wk)*(1-wa)
!    |    .      |              w2 = wk*(1-wa)
!    |. . + . . .| wa2   A      w3 = (1-wk)*wa
!    |    .      |       |      w4 = wk*wa
!    |    .      |       wa
!    |    .      |       |
!    1-----------2 ja2   V
!   jk2  wk2  jk2p
!
!    <-wk->
!
!
t2 = 1.
t4 = 1.
!
!-----------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------
!  Main loop over the locus
!
do iloc=1,nlocusx
!
  jk2  = t_ik2(iloc)
  jk2p = min(jk2+1,nkq)
  ja2  = mod(t_ia2(iloc)-1+naq,naq)+1
  ja2p = mod(t_ia2(iloc)+naq,naq)+1
!
! compute tail parameters
!
!!  if(iq_geom==1) then
!!    jk2 = max(1,jk2)
!!    jk4 = max(1,jk4)
!!    t2 = max(1.,q_kfac**real(t_ik2(iloc)-nkq))
!!    t2 = t2**qk_tail
!!    t4 = max(1.,q_kfac**real(t_ik4(iloc)-nkq))
!!    t4 = t4**qk_tail
!!  end if
!---------------------------------------------------------------------------------------
!  check boundaries of sector grid
!
  if(iq_grid < 3) then
    ja2  = max(ja2,1)
    ja2  = min(ja2,naq)
    ja2p = max(ja2p,1)
    ja2p = min(ja2p,naq)
  end if
!
  qn2  = (t_w1k2(iloc)*nspec(jk2,ja2)  + t_w2k2(iloc)*nspec(jk2p,ja2) + &
&         t_w3k2(iloc)*nspec(jk2,ja2p) + t_w4k2(iloc)*nspec(jk2p,ja2p))*t2
!
  jk4  = t_ik4(iloc)
  jk4p = min(jk4+1,nkq)
  ja4  = mod(t_ia4(iloc)-1+naq,naq)+1
  ja4p = mod(t_ia4(iloc)+naq,naq)+1
!
!  special treatment for sector grids
!  limit range of indices
!  QQQ: in fact energy density should be set to ZERO
!
  if(iq_grid < 3) then
    ja4  = max(ja4,1)
    ja4  = min(ja4,naq)
    ja4p = max(ja4p,1)
    ja4p = min(ja4p,naq)
  end if
!
  qn4  = (t_w1k4(iloc)*nspec(jk4,ja4)  + t_w2k4(iloc)*nspec(jk4p,ja4) + &
&         t_w3k4(iloc)*nspec(jk4,ja4p) + t_w4k4(iloc)*nspec(jk4p,ja4p))*t4
!
!-------------------------------------------------------------------------------
!
  nprod     = qn13p*(qn4-qn2) + qn2*qn4*qn13d
  rterm     = t_zz(iloc)
  t13        = t13 + rterm*nprod
!
! output to triplets
!
!
!
!  add diagonal terms
!
!!  qd1    = qn3*(qn4-qn2) + qn2*qn4*qn3
!!  qd3    = qn1*(qn4-qn2) - qn2*qn4*qn1
!
  qd1    = qn3*(qn4-qn2) - qn2*qn4
  qd3    = qn1*(qn4-qn2) + qn2*qn4
  diagk1 = diagk1 + qd1*rterm
  diagk3 = diagk3 + qd3*rterm
!-----------------------------------------------------------------------------------
end do
!
!!/T if(iq_test>=4) then
!!/T   write(luq_tst,'(a)') 'Q_T13V4: NSPEC'
!!/T  do ikq=1,nkq
!!/T    write(luq_tst,'(100e12.4)') (nspec(ikq,iaq),iaq=1,naq)
!!T  end do
!!T end if
!!if(iq_integ==3) write(luq_int,'(4i3,i5,1000e13.5)') ik1,ia1,ik3,ia3,nloc, &
!!& t_s(nloc),t13,(dt13(iloc),iloc=1,nloc)
!
9999 continue
!
!!!!!!!!!!!!!!call q_stack0('-q_t13v4')
!
return
end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!--------------------------------!!!!!
!!-------------Q-GETLOCUS---------------!!
!!!!!--------------------------------!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine q_getlocus(ik1,ia1,ik3,ia3,ifnd)
!------------------------------------------------------------------------------

!   +-------+    ALKYON Hydraulic Consultancy & Research
!   |       |    Gerbrant van Vledder
!   |   +---+
!   |   | +---+  Last update: 27 August 2003
!   +---+ |   |  Release: 5.0
!         +---+
!
! do not use m_xnldata
use m_constants
use serv_xnl4v5
!----------------------------------------------------------------------------------
implicit none
!
!  0. Update history
!
!    25/02/1999  Initial version
!    15/04/1999  Extra parameter IFND to indicate if a
!                reference locus exists in the data base
!    19/07/1999  Restructured and some bugs removed
!    20/07/1999  Option added to compute locus directly, no database
!                or scaling involved. IFND < 0
!    15/10/1999  Information to transformation file updated
!    16/10/1999  Equation for computing address for storage updated
!    25/10/1999  Transformations updated
!    28/10/1999  Local variables ia1 and ia3 may not be changed, temp. variables
!                it1 and it3 included
!    29/10/1999  Use of IQ_TRF modified
!                EPSK introduced to check equality of loci
!    28/12/1999  A_CMPLOC renamed to Q_CMPLOC
!    03/01/2000  IQ_START replaced by IQ_LOCUS
!    05/01/2000  Interface with Q_CMPLOC modified
!    08/08/2002  Upgrade to release 4.0
!    13/08/2002  Indexing in terms of integers and real weights
!                upgrade to release 4.0
!    19/08/2002  Bug fixed in transforming CASE 6
!                Interpolation option added
!    12/06/2003  Parameter t_ws set equal to r_ws
!    13/06/2003  Parameter t_cple, t_jac and t_sym assigned
!                Bug fixed in nearest bin approach, symmetry regained
!    27/08/2003  Short-cut when number of points on locus is ZERO
!
!  1. Purpose:
!
!     Retrieve locus from basic locus as stored in the database
!
!  2. Method
!
!     In the case of geometric scaling, k-scaling is used using scale laws
!     described by Tracy
!
!     Directional transformation using linear transformations, shifting and mirror
!     imaging.
!
!
!  3. Parameter list:
!
!Type      I/O           name       Description
!------------------------------------------------------------------------------
integer, intent(in)  ::  ik1    !  k-index of wave number k1
integer, intent(in)  ::  ia1    !  theta-index of wave number k1
integer, intent(in)  ::  ik3    !  k-index of wave number k3
integer, intent(in)  ::  ia3    !  theta-index of wave number k3
integer, intent(out) :: ifnd    !  indicator if reference locus exists in database
!
!  4. Error messages
!
!  5. Called by
!
!     Q_T13V4
!
!  6. Subroutines used
!
!  7. Remarks
!
!  8. Structure
!
!  9. Switches
!
! 10. Source code
!------------------------------------------------------------------------------
!     Local variables
!
integer it1,it3           ! work indices for directions, copy of ia1 and ia3
integer idir              ! switch to indicate if locus should be inverted
integer itrans            ! type of transformation
integer iloc,jloc         ! counters along locus
integer iadif,ikdif       ! difference in angular and k-index
integer ja1r,ja3r
integer imirror           ! extra step when locus is mirrorred
!
integer ibeta,kdif
integer nloc             ! number of points on locus
!
integer ierr
integer kmem                 ! index for storing 2-d matrix in 1-d array
integer amem                 ! index for storing direction of reference wave number k3
!
real lambda                  ! geometric scale factor
real j_lambda                ! scale factor for Jacobian term
real c_lambda                ! scale factor for coupling coefficient
real zz_lambda               ! combined scale factor
!
real xt2(nlocus),yt2(nlocus) ! xy-components of test k2-locus
real xt4(nlocus),yt4(nlocus) ! xy-components of test k4-locus
real wk,wa,vk,va
!! \A
!! real x_kfunc                 ! real function to compute wave number
!! \Z
!
integer ikmin,ja1,ja3,jk1,jk3,itmin
integer ibdif,nhalf
integer iaq,ikq              ! counters for loop over direction and wave numbers
!------------------------------------------------------------------------------
!
!! data i_getloc /0/  ! Initialise counter
!-------------------------------------------------------------------------------------
!!!!!!!!!!!call q_stack('+q_getlocus')
!
!------------------------------------------------------------------------------
! initialisations
!------------------------------------------------------------------------------
!
it1 = ia1
it3 = ia3
!
imirror = 1
!
ikmin = min(ik1,ik3)    ! compute minimum of wave number index
ikdif = abs(ik1-ik3)    ! compute difference between wave number indices
!
if (iq_geom ==0) then
  jk1 = min(ik1,ik3)
  jk3 = max(ik1,ik3)
else
  jk1 = 1
  jk3 = ikdif + 1       ! compute k-index of wave number k3 relative to reference wave number
end if
!
itmin = min(it1,it3)    ! compute minimum angle of k1 and k3
iadif = abs(it1-it3)    ! difference index
ja1   = 1               ! index of direction of reference wave number k1
ja3   = iadif+iaref     ! compute theta-index of direction of wave number k3
!
!------------------------------------------------------------------------------
!  circle grid, modify ranges and transformation variables
!------------------------------------------------------------------------------
!
if (iq_grid==3) then
  nhalf = naq/2
  if (iadif > nhalf) then
    if(it1 > nhalf) it1 = it1 - naq
    if(it3 > nhalf) it3 = it3 - naq
  end if
  itmin = min(it1,it3)
  ibdif = (naq - abs(naq-2*abs(it1-it3)))/2   ! compute shortest difference in indices
                                              ! while taking care of periodicity
  ja3   = ibdif + 1
  iadif = ibdif
end if
!
ja1r = 1
ja3r = iadif + 1
amem = iadif + 1   ! compute index of reference wave number k3 in interaction grid
!
!------------------------------------------------------------------------------
! obtain k-index of reference wave number
!------------------------------------------------------------------------------
!
if(iq_geom==0) then
  kmem = (jk3-jk1+1) - (jk1-2*nkq-2)*(jk1-1)/2
else
  kmem = ikdif+1
end if
!
!
!------------------------------------------------------------------------------
!  check memory indexing
!------------------------------------------------------------------------------
!
if (amem > iamax) then
  ifnd = 0
  !!!!call q_error('e','MEMORY','Incorrect addres')
  print *, "Error 13879"
  goto 9999
end if
!
!-----------------------------------------------------------------------------
! retrieve info from reference locus in
! get actual number of valid points along locus (NLOCUSZ)
! depending on value of switch IQ_COMPACT
!------------------------------------------------------------------------------
!
nloc    = quad_nloc(kmem,amem)
nlocusx = nloc
!
!  short-cut when number of NON-ZERO points on locus is ZERO [27/8/2003]
!
if(nlocusx==0) goto 9999
!
r_ik2(1:nloc)  = quad_ik2(kmem,amem,1:nloc)
r_ia2(1:nloc)  = quad_ia2(kmem,amem,1:nloc)
r_ik4(1:nloc)  = quad_ik4(kmem,amem,1:nloc)
r_ia4(1:nloc)  = quad_ia4(kmem,amem,1:nloc)
!
r_w1k2(1:nloc) = quad_w1k2(kmem,amem,1:nloc)
r_w2k2(1:nloc) = quad_w2k2(kmem,amem,1:nloc)
r_w3k2(1:nloc) = quad_w3k2(kmem,amem,1:nloc)
r_w4k2(1:nloc) = quad_w4k2(kmem,amem,1:nloc)
!
r_w1k4(1:nloc) = quad_w1k4(kmem,amem,1:nloc)
r_w2k4(1:nloc) = quad_w2k4(kmem,amem,1:nloc)
r_w3k4(1:nloc) = quad_w3k4(kmem,amem,1:nloc)
r_w4k4(1:nloc) = quad_w4k4(kmem,amem,1:nloc)
!
r_zz(1:nloc)   = quad_zz(kmem,amem,1:nloc)
!
!------------------------------------------------------------------------------
kdif = ikmin - 1
if(iq_geom==0) then
  lambda = 1.
  kdif  = 0.
else
  lambda = q_kfac**(ikmin-1.)
end if
!
j_lambda = 1./sqrt(lambda)
c_lambda = lambda**6
!
!  compute combined scale factor
!
zz_lambda = lambda*c_lambda/j_lambda
!
!------------------------------------------------------------------------------
! select case to transform reference locus
!
! Transform of weigths reduces to an addition or subtraction
! because of log-spacing of wave numbers in the case of deep water
!
if(ik3 > ik1 .and. it3 >= it1) then      ! Case 1
  itrans = 1
  t_ik2(1:nloc)  = kdif + r_ik2(1:nloc)
  t_ik4(1:nloc)  = kdif + r_ik4(1:nloc)
  ibeta          = itmin-iaref
  t_ia2(1:nloc)  = r_ia2(1:nloc) + ibeta
  t_ia4(1:nloc)  = r_ia4(1:nloc) + ibeta
  idir   = 1
  t_w1k2(1:nloc)  = r_w1k2(1:nloc)
  t_w2k2(1:nloc)  = r_w2k2(1:nloc)
  t_w3k2(1:nloc)  = r_w3k2(1:nloc)
  t_w4k2(1:nloc)  = r_w4k2(1:nloc)
  t_w1k4(1:nloc)  = r_w1k4(1:nloc)
  t_w2k4(1:nloc)  = r_w2k4(1:nloc)
  t_w3k4(1:nloc)  = r_w3k4(1:nloc)
  t_w4k4(1:nloc)  = r_w4k4(1:nloc)
!
elseif(ik3 > ik1 .and. it3 < it1)  then  ! Case 2
  itrans = 2
  t_ik2(1:nloc)  = kdif + r_ik2(1:nloc)
  t_ik4(1:nloc)  = kdif + r_ik4(1:nloc)
  ibeta          = int(q_ad(ia1)/q_deltad+0.01)
  t_ia2(1:nloc)  = ibeta + 2.*iaref - r_ia2(1:nloc) -imirror
  t_ia4(1:nloc)  = ibeta + 2.*iaref - r_ia4(1:nloc) -imirror
  t_w1k2(1:nloc)  = r_w3k2(1:nloc)
  t_w2k2(1:nloc)  = r_w4k2(1:nloc)
  t_w3k2(1:nloc)  = r_w1k2(1:nloc)
  t_w4k2(1:nloc)  = r_w2k2(1:nloc)
  t_w1k4(1:nloc)  = r_w3k4(1:nloc)
  t_w2k4(1:nloc)  = r_w4k4(1:nloc)
  t_w3k4(1:nloc)  = r_w1k4(1:nloc)
  t_w4k4(1:nloc)  = r_w2k4(1:nloc)
  idir   = -1   ! according to theory
!  idir   = 1    ! as it should be to get symmetry
!
elseif(ik1 > ik3 .and. it3 >= it1) then      ! Case 3
  itrans = 3
  t_ik2(1:nloc)  = kdif + r_ik4(1:nloc)
  t_ik4(1:nloc)  = kdif + r_ik2(1:nloc)
  ibeta          = int(q_ad(ia3)/q_deltad+0.01)
  t_ia2(1:nloc)  = ibeta + 2.*iaref - r_ia4(1:nloc) -imirror
  t_ia4(1:nloc)  = ibeta + 2.*iaref - r_ia2(1:nloc) -imirror
  t_w1k2(1:nloc)  = r_w3k2(1:nloc)
  t_w2k2(1:nloc)  = r_w4k2(1:nloc)
  t_w3k2(1:nloc)  = r_w1k2(1:nloc)
  t_w4k2(1:nloc)  = r_w2k2(1:nloc)
  t_w1k4(1:nloc)  = r_w3k4(1:nloc)
  t_w2k4(1:nloc)  = r_w4k4(1:nloc)
  t_w3k4(1:nloc)  = r_w1k4(1:nloc)
  t_w4k4(1:nloc)  = r_w2k4(1:nloc)
  idir   = 1
!
elseif(ik1 > ik3 .and. it1 > it3) then   ! Case 4
  itrans = 4
  t_ik2(1:nloc)  = kdif + r_ik4(1:nloc)
  t_ik4(1:nloc)  = kdif + r_ik2(1:nloc)
  ibeta          = itmin-iaref
  t_ia2(1:nloc)  = ibeta + r_ia4(1:nloc)
  t_ia4(1:nloc)  = ibeta + r_ia2(1:nloc)
  idir   = -1
  t_w1k2(1:nloc)  = r_w1k2(1:nloc)
  t_w2k2(1:nloc)  = r_w2k2(1:nloc)
  t_w3k2(1:nloc)  = r_w3k2(1:nloc)
  t_w4k2(1:nloc)  = r_w4k2(1:nloc)
  t_w1k4(1:nloc)  = r_w1k4(1:nloc)
  t_w2k4(1:nloc)  = r_w2k4(1:nloc)
  t_w3k4(1:nloc)  = r_w3k4(1:nloc)
  t_w4k4(1:nloc)  = r_w4k4(1:nloc)
!
elseif(ik1==ik3 .and. it3 > it1) then  ! Case 5
  itrans = 5
  t_ik2(1:nloc)  = kdif + r_ik2(1:nloc)
  t_ik4(1:nloc)  = kdif + r_ik4(1:nloc)
  ibeta          = itmin-iaref
  t_ia2(1:nloc)  = r_ia2(1:nloc) + ibeta
  t_ia4(1:nloc)  = r_ia4(1:nloc) + ibeta
  idir   = 1
  t_w1k2(1:nloc)  = r_w1k2(1:nloc)
  t_w2k2(1:nloc)  = r_w2k2(1:nloc)
  t_w3k2(1:nloc)  = r_w3k2(1:nloc)
  t_w4k2(1:nloc)  = r_w4k2(1:nloc)
  t_w1k4(1:nloc)  = r_w1k4(1:nloc)
  t_w2k4(1:nloc)  = r_w2k4(1:nloc)
  t_w3k4(1:nloc)  = r_w3k4(1:nloc)
  t_w4k4(1:nloc)  = r_w4k4(1:nloc)
!
elseif(ik1==ik3 .and. it1 > it3) then  ! Case 6
  itrans = 6
  t_ik2(1:nloc)  = kdif + r_ik4(1:nloc)
  t_ik4(1:nloc)  = kdif + r_ik2(1:nloc)
  ibeta          = int(q_ad(ia1)/q_deltad+0.01)
  t_ia2(1:nloc)  = ibeta + 2.*iaref - r_ia2(1:nloc) -imirror
  t_ia4(1:nloc)  = ibeta + 2.*iaref - r_ia4(1:nloc) -imirror
!!  ibeta          = itmin-iaref
!!  t_ia2(1:nloc)  = r_ia4(1:nloc) + ibeta
!!  t_ia4(1:nloc)  = r_ia2(1:nloc) + ibeta
  idir   = -1
  t_w1k2(1:nloc)  = r_w3k2(1:nloc)
  t_w2k2(1:nloc)  = r_w4k2(1:nloc)
  t_w3k2(1:nloc)  = r_w1k2(1:nloc)
  t_w4k2(1:nloc)  = r_w2k2(1:nloc)
  t_w1k4(1:nloc)  = r_w3k4(1:nloc)
  t_w2k4(1:nloc)  = r_w4k4(1:nloc)
  t_w3k4(1:nloc)  = r_w1k4(1:nloc)
  t_w4k4(1:nloc)  = r_w2k4(1:nloc)
end if
!
t_zz(1:nloc)   = lambda*c_lambda/j_lambda * r_zz(1:nloc)
!
ifnd = 1
!
!------------------------------------------------------------------------------
!
9999 continue
!
!!!!!!!!!!!!call q_stack('-q_getlocus')
!
return
end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!--------------------------------!!!!!
!!-------------Q-INIT-------------------!!
!!!!!--------------------------------!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


subroutine q_init
!------------------------------------------------------------------------------
!
!   +-------+    ALKYON Hydraulic Consultancy & Research
!   |       |    Gerbrant van Vledder
!   |   +---+
!   |   | +---+  Last update: 25 Sep. 2002
!   +---+ |   |  Release: 5.0
!         +---+
!
! do not use m_xnldata
use m_fileio
use m_constants
use serv_xnl4v5
implicit none
!--------------------------------------------------------------------------------
!  0. Update history
!
!     25/02/1999 Initial version
!     13/10/1999 Error handling improved
!     18/10/1999 Test output to MATCHK.GRD added
!     21/10/1999 Extra output to MATCHK.GRD, iaref=1 for circle grids
!     01/11/1999 Allocatable arrays Q_XK and Q_SK added
!     14/02/2001 Version ready for WAVEWATCH III
!      8/08/2002 Release 4.
!     16/08/2002 Group velocity computed
!     22/08/2002 First and last used defined direction accounted for
!     11/09/2002 Call of Q_ALOC moved to higher level, viz. XNL_INIT
!                q_kpow initialized
!     25/09/2002 User defined directions used in the case of a sector grid
!
!  1. Purpose:
!
!     Initializing module for quadruplets
!     and setting default settings
!
!  2. Method
!
!     Conversion of power of spectral tail from E(f) to N(k) using the following
!     relations:
!
!       E(f) ~ f^qf_tail
!
!       N(k) ~ k^qk_tail
!
!       qk_tail = qf_tail/2 -1
!
!     See also Note 13 of G.Ph. van Vledder
!
!  3. Parameter list:
!
!     Name    I/O  Type  Description
!
!
!  4. error meassage
!
!  5. Called by
!
!     XNL_INIT
!
!  6. Subroutines used
!
!     Q_STACK
!     Z_CMPCG
!     Z_STEPS
!     Z_WNUMB
!
!  7. Remarks
!
!  8. Structure
!
!  9. Switches
!
!     /S Enable subroutine tracing
!
! 10. Source code
!------------------------------------------------------------------------------------------
!     Local variables
!
integer iaq,ikq    ! counters for loops over directions and wave numbers
real ff            ! frequency
!!!real z_wnumb       ! service function to compute wave number
!
!------------------------------------------------------------------------------
!
!!!!!!!!!!!!!!call q_stack('+q_init')
!
! set general settings
!
! convert power of E(f) f^qf_tail to power of N(k) k^qk_tail
! See Note 13 of G.Ph. van Vledder
!
qk_tail = (qf_tail-2.)/2. ! power of spectral tail, of N(k)
!
if(iq_prt >=2) then
  write(luq_prt,*)
  write(luq_prt,'(a,f6.1)') 'Q_INIT:  E(f)_tail: ',qf_tail
  write(luq_prt,'(a,f6.1)') 'Q_INIT:  N(k)_tail: ',qk_tail
end if
!
! set absolute and relative accuracies
!
eps_q   = 0.001           ! absolute accuracy for check of q==0
eps_k   = 1.e-5           ! absolute accuracy for equality check of k
rel_k   = 0.001           ! relative accuracy for equality check of k
!
sk_max = 50.              ! set maximum waver number
wk_max = real(nkq+0.9999) ! set maximum wave number index
!
! compute frequency and wave number grid
! assume that frequencies are always geometrically spaced,
! in the case of deep water this also holds for the wave numbers
!
q_ffac = (fqmax/fqmin)**real(1./(nkq-1.))         ! geometric spacing factor of frequencies
!
ff = fqmin                                        ! set minimum frequency
!
if(iq_prt>=2) then
  write(luq_prt,*)
  write(luq_prt,'(a)') 'Basic wave numbers, frequencies'
end if
!
do ikq=1,nkq                                       ! Generate wave number dependent variables
  q_f(ikq)    = ff                                 ! Frequency
  q_sig(ikq)  = ff*2.*pi                           ! Radian frequency
  q_k(ikq)    = z_wnumb(q_sig(ikq),q_depth,q_grav) ! compute wave number
  q_kpow(ikq) = (q_k(1)/q_k(ikq))**7.5             ! used in filtering
  ff          = ff*q_ffac                          ! Increase frequency
!
  call z_cmpcg(q_sig(ikq),q_depth,q_grav,q_cg(ikq))
  if(iq_prt >= 2) then
    write(luq_prt,'(a,i4,3f10.5,e12.4)') 'Q_INIT: ikq f sigma k k^p:', &
&    ikq,q_f(ikq),q_sig(ikq),q_k(ikq),q_kpow(ikq)
  end if
end do
!
! compute characteristics of extended k-array
!
if(iq_prt>=2) then
  write(luq_prt,*)
  write(luq_prt,'(a)') 'Extended wave numbers and spacing'
end if
!
do ikq=0,nkq
  if(ikq==0) then
    q_xk(ikq) = 0.
    q_sk(ikq) = q_k(1)
  elseif(ikq==nkq) then
    q_xk(ikq) = q_k(ikq)
    q_sk(ikq) = sk_max
  else
    q_xk(ikq) = q_k(ikq)
    q_sk(ikq) = q_k(ikq+1) - q_k(ikq)
  end if
!
end do
!
!
kqmin = q_k(1)
kqmax = q_k(nkq)
q_kfac = (kqmax/kqmin)**real(1./(nkq-1))  ! this value makes only sense in the
                                          ! case of deep water, IQ_DISP==1
!
! compute step size of frequency grids and wave number grid
!
call z_steps(q_f,  q_df,  nkq)           ! step size of frequencies
call z_steps(q_sig,q_dsig,nkq)           ! step size of radian frequencies
call z_steps(q_k,  q_dk,  nkq)           ! step size of wave numbers
!
if(iq_prt >= 2) then
  write(luq_prt,*)
  write(luq_prt,'(a)') 'Q_INIT: Additional information'
  write(luq_prt,'(a,f8.1)')  'Q_depth (m):',q_depth
  write(luq_prt,'(a,i3)')    'Number of frequencies:',nkq
  write(luq_prt,'(a,f8.4)')  'Geometric f-spacing factor:',q_ffac
  write(luq_prt,'(a,f8.4)')  'Geometric k-spacing factor:',q_kfac
  write(luq_prt,'(a,2f8.3)') 'fmin fmax (Hz):',fqmin,fqmax
  write(luq_prt,'(a,2f8.3)') 'kmin kmax (Hz):',kqmin,kqmax
  write(luq_prt,*)
!
  write(luq_prt,*) '     i      f         df       sig      dsig       k         dk         cg'
!
  do ikq=1,nkq
    write(luq_prt,'(1x,i4,7f10.4)') &
 &  ikq,q_f(ikq),q_df(ikq),q_sig(ikq),q_dsig(ikq),q_k(ikq),q_dk(ikq),q_cg(ikq)
  end do
end if
!
! =============== D I R E C T I O N S ===============================================
!
! the directions in the array ANGLE are running from 1 to NAQ
! for a sector definition the middle direction has index IAREF
!
!  compute index IAREF of middle wave direction for sector grids
!
if(iq_grid ==1 .or. iq_grid==2) then
  iaref = (naq/2)+1
elseif(iq_grid==3) then
  iaref = 1
end if
!
if(iq_prt >= 2) write(luq_prt,'(a,i4)') &
&  'Q_INIT: Index of first direction for reference:',iaref
!
!  set loops indices
!
if(iq_grid==1) then    ! symmetric sector
  iaq1 = iaref
  iaq2 = naq
!
! non-symmetric sector and full circle
!
elseif(iq_grid==2 .or. iq_grid==3) then
  iaq1 = 1
  iaq2 = naq
end if
!
if(iq_prt >= 2) write(luq_prt,'(a,2i4)') &
&  'Q_INIT: Range of indices for loop over directions:',iaq1,iaq2
!
!  generate directions, given in degrees
!
q_sector = 0.5*(abs(q_dird1) + abs(q_dird2))
!
if(iq_grid==1 .or. iq_grid==2) then    ! define symmetric sector
  q_deltad = 2.*q_sector/real(naq-1.)  ! delta in degrees
  q_ang1 = -q_sector                   ! degrees
  q_ang2 =  q_sector                   ! degrees
  if(iq_prt>0) write(luq_prt,'(a)') 'Q_INIT: take care of q_dird1 and check if sector is OK'
!
elseif(iq_grid==3) then                ! full sector
  q_deltad = 360./real(naq)            ! degrees
  q_ang1 = 0                           ! degrees
  q_ang2 = 360.-q_delta                ! degrees
end if
!
q_delta = q_deltad*dera                ! directional step in radians
ncirc   = 2.00001*pi/q_delta           ! number of directions on circle
!
if(iq_prt >= 2) then
  write(luq_prt,'(a,3f10.3)')     'Q_INIT: d(1),d(n),dsector:',q_dird1,q_dird2,q_sector
  write(luq_prt,'(a,f6.2,a)')     'Q_INIT: Angular step     :',q_deltad,' degrees'
  write(luq_prt,'(a,2f8.2,i4,a)') 'Q_INIT: ang1 ang2 nang   :',q_ang1,q_ang2,naq,' degrees'
  write(luq_prt,'(a,i4)')         'Q_INIT: #Angles on circle:',ncirc
  write(luq_prt,*)
end if
!
!  generate directions arrays, given in degrees and radians
!
do iaq=1,naq
  q_ad(iaq) = q_ang1 + q_deltad*(iaq-1.)
  q_a(iaq)  = q_ad(iaq)*dera
  if(iq_prt >= 2) then
    write(luq_prt,'(a,i4,f10.4,f10.2)') 'Q_INIT: iaq q_a q_ad:',iaq,q_a(iaq),q_ad(iaq)
    if(iaq==naq) write(luq_prt,*)
  end if
end do
!
!  set loop indices for generation of grid
!  for sector grids and circle grids
!
if(iq_grid==1 .or. iq_grid==2) then
  iag1  = iaref
  iag2  = naq
!
!  circle grid
!
elseif(iq_grid==3) then
  iag1 = 1
  iag2 = naq/2+1
end if
!
iamax = iag2-iag1+1
!-------------------------------------------------------------------------
!
!
!!!!!!!!!!!!!!call q_stack('-q_init')
!
return
end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!--------------------------------!!!!!
!!-------------Q-CTRGRID----------------!!
!!!!!--------------------------------!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine q_ctrgrid(itask,igrid)
!------------------------------------------------------------------------------
!
!   +-------+    ALKYON Hydraulic Consultancy & Research
!   |       |    Gerbrant van Vledder
!   |   +---+
!   |   | +---+  Last update: 13 Sept. 2003
!   +---+ |   |  Release: 5.03
!         +---+
!
! do not use m_xnldata
use m_fileio
implicit none
!------------------------------------------------------------------------------
!  0. Update history
!
!     Version  Date    Modification
!
!     29/07/1999  Initial version
!     11/10/1999  Error messages via module
!     12/10/1999  File I/O of interaction grid files added, and consistency check
!     25/10/1999  Contents of q_header extended with iq_grid & iq_start & nlocus0
!     27/10/1999  Close statement after reading BQF file
!     01/11/1999  Close statments after call of Q_GRID
!     03/11/1999  Parameter IQ_MAKE included
!     22/11/1999  Use of z_fileio modified, use parameter IUERR if an attempt
!                 to open a non-existing file was made
!     30/11/1999  Extra messages to logging file when file are closed
!     03/01/2000  IQ_START replaced by IQ_LOCUS
!     12/01/2001  Output parameter IGRID added
!                 igrid=0: a proper grid exists or has been made, and will be read
!                      =1: grid file does not exist
!                      =2: grid file exists, but it is incorrect
!      8/08/2002  Upgrade to release 4.0
!     15/08/2002  Bug fixed ininitialisation of igrid
!     19/08/2002  header extended with parameter IQ_INTERP
!     20/08/2002  wave number array replaced by sigma array in grid file
!     22/08/2002  Extra i/o check when reading BQF file
!     23/08/2002  retrieve number of point on locus from BQF file
!     09/09/2002  aqfile and bqfile 5 units for depth
!     10/09/2002  new algorithm for coding depth
!                 Test added to avoid rereading of last read/generated BQF file
!     08/10/2002  Output to test file made conditional
!     05/09/2003  Water depth for creating and testing BQF file DSTEP dependend
!     09/09/2003  Bug fixed in assigning IGRID=0 when BQF still in memory
!     13/09/2003  When BFQ incorrupt, it is deleted and a new one is created
!                 Bug fixed in setting of s_depth  when iq_disp==1
!
!  1. Purpose:
!
!     Control of interaction grid administration
!
!  2. Method
!
!  3. Parameters used
!
integer, intent(in)  :: itask  !  task to perform by Q_CTRGRID
!                                 ==1: read and check header block
!                                 ==2: read and write grid file, according to
!                                      setting of IQ_MAKE
integer, intent(out) :: igrid  !  status of grid checking
!                                 ==0: a proper grid exists
!                                 ==1: grid file does not exist
!                                 ==2: grid file exists, but it is incorrect
!                                 ==3: read error in accessing grid information
!
!  4. Error messages
!
!  5. Called by:
!
!     XNL_INIT
!     Q_SEARCHGRID
!
!  6. Subroutines used
!
!     Q_STACK
!     Q_ERROR
!     Q_MAKEGRID
!
!  7. Remarks
!
!     The generation of the database file depend on the control varaible of IQ_MAKE
!     if IQ_MAKE==1, make a grid when needed
!                 2, always make grid
!                 3, make a grid and stop, useful for test purposes
!
!     The maximum number of points on the locus, as stored in the BQF file
!     is read from the header and stored in the variable NLOCUS
!
!  8. Structure
!
!     Make header of grid file
!     Construct name of grid file
!     Check existence of grid file
!     if grid file exists
!       read header
!       check header
!       read first data block
!       check first data block
!       - directions
!       - wave numbers
!       - depth
!       set status of grid file
!     else
!       set status of grid file
!     end if
!
!     set status of generating/reading grid file
!
!     if make new grid
!       compute grid parameters
!       write grid parameters to file
!     else
!       read grid parameters from file
!     end if
!
!
!  9. Switches
!
! 10. Source code
!-------------------------------------------------------------------------------
!     Local variables
!
integer iaz,ikz,jaz,jkz                ! counters for checking header of BQF file
integer iz_geom,iz_disp,iz_cple        ! values as read in BQF file
integer naz                            ! number of directions in BQF file
integer nkz                            ! number of wave numbers in BQF file
integer idep,jdep                      ! coding for depth in BQF file
!
logical lbqf                           ! flag for existence of BQF file
real s_depth                           ! (stepped) depth
real q_depth_saved                     ! save input water depth, needed after call to Q_MAKEGRID
real z_depth                           ! water depth in BQF file
real, allocatable :: z_ad(:),z_sig(:)  ! directions and radian frequencies of grid in BQF file
integer ierr,iuerr                     ! error variables
!------------------------------------------------------------------------------
!
!!!!!!!!!!!!!!!!call q_stack('+q_ctrgrid')
!
!  echo input arguments
!
!
q_depth_saved = q_depth
!
!  generate header of BQF file
!
q_header = '000000000000000000000'
!           123456789012345678901
!                    1         2
write(q_header,'(3i3.3,6i2.2)') naq,nkq,nlocus0,&
&  iq_grid,iq_geom,iq_disp,iq_cple,iq_locus,iq_interp
!
if(iq_prt>=2) then
  write(luq_prt,'(2a)')    'Q_CTRGRID: header info:',trim(q_header)
  write(luq_prt,'(a,3i5)') 'Q_CTRGRID: naq nkq nlocus0:',naq,nkq,nlocus0
  write(luq_prt,'(a,3i3)') 'Q_CTRGRID: iq_grid,iq_geom,iq_disp:',iq_grid,iq_geom,iq_disp
  write(luq_prt,'(a,3i3)') 'Q_CTRGRID: iq_cple,iq_locus,iq_interp:',iq_cple,iq_locus,iq_interp
end if
!
!------------------------------------------------------------------------------
! construct name of grid file
!
if(iq_disp==1) then
  bqname = 'quad99999.bqf'
  s_depth = q_maxdepth
!
elseif(iq_disp==2) then
!
!---------------------------------------------------------------------------------------------
! generate code for actual depth
!---------------------------------------------------------------------------------------------
  idep = int(q_depth/q_dstep+0.5)
  jdep = idep*int(10.*q_dstep)
  jdep = max(1,jdep)
  jdep = min(99999,jdep)
!
  s_depth = real(idep)*q_dstep
!
!
  bqname = 'quad00000.bqf'
  write(bqname(5:9),'(i5.5)') min(int(q_maxdepth*10),jdep)
!
else
  !!!!!!!!!!!!!!call q_error('e','DISPER','Incorrect value for IQ_DISP')
  print *, "ERRRROOOOR"
  write(luq_err,'(a,i4)') 'IQ_DISP=',iq_disp
  goto 9999
end if
!
!
!-------------------------------------------------------------------------------------------
!  Compare LASTQUADFILE with bqname
!  if equal skip reading of BQF file, including checks of header
!-------------------------------------------------------------------------------------------
!
if(lastquadfile==bqname) then
  if(iq_screen>0) write(iscreen,'(2a)')   'Q_CTRGRID: Rereading of bqfile skipped: ',lastquadfile
  igrid = 0
  goto 9999
end if
!-------------------------------------------------------------------------------------------
if(iq_prt >= 2) then
  write(luq_prt,'(2a)') 'Q_CTRGRID: Header line of grid file:',trim(q_header)
  write(luq_prt,'(2a)') 'Q_CTRGRID: Name of BINARY grid file:',trim(bqname)
end if
!------------------------------------------------------------------------------
!
!  check if binary data file exists
!
call z_fileio(bqname,'OU',iufind,luq_bqf,iuerr)   ! binary quadruplet file
!
if(iq_prt >= 2) write(luq_prt,'(2a,2x,2i4)') &
& 'Q_CTRGRID:  bqname:',trim(bqname),luq_bqf,iuerr
!
if(itask==2 .and. iq_make==2) luq_bqf=-1
!
!  if the file exists,
!     read header information
!     check header of file
!  end
!
!  If header is incorrect, set flag IQ_GRID to TRUE for generating new grid
!
if(luq_bqf > 0 .and. iuerr ==0) then
  if(iq_prt >= 2) then
    write(luq_prt,'(2a)')   'Q_CTRGRID: Binary grid file detected: ',trim(bqname)
    write(luq_prt,'(a,i4)') 'Q_CTRGRID: Connected to unit:',luq_bqf
  end if
!
!
! grid file exists, unless proven otherwise
!---------------------------------------------------------------------------------
!
  lq_grid = .false.
  igrid   = 0
  read(luq_bqf,iostat=ierr) r_header
  if(ierr/=0) then
    !!!!!!!!call q_error('w','READBQF','Read error for header in BQF file')
    print *, "Another error!"
    write(luq_err,'(a)') 'BQF file deleted'
    call z_fileio(bqname,'DU',iufind,luq_bqf,iuerr)   ! binary quadruplet file
    igrid = 3
    lq_grid = .true.
  else
    read(r_header,'(6x,i3)') nlocus
!
  end if
!-----------------------------------------------------------------------------
!   check header of grid file
!-----------------------------------------------------------------------------
!
  if(trim(r_header)/=trim(q_header).and. .not.lq_grid) then
    lq_grid = .true.
    igrid   = 2
    if(iq_prt >=2) then
      write(luq_prt,'(a,1x,a)') &
&     'Q_CTRGRID: Header in binary quad file         :',trim(r_header)
      write(luq_prt,'(a,1x,a)') &
&     'Q_CTRGRID: Expected header of binary quad file:',trim(q_header)
      write(luq_prt,'(a)') 'Q_CTRGRID: The file headers disagree'
      write(luq_prt,'(a)') 'Q_CTRGRID: A new grid will be generated'
    end if
  end if
!------------------------------------------------------------------------------
!  check other parts of binary grid file
!
  if(.not.lq_grid) then
    read(luq_bqf) naz,nkz
    allocate (z_sig(nkz),z_ad(naz))
    read(luq_bqf) z_sig
    read(luq_bqf) z_ad
    read(luq_bqf) iz_geom,iz_disp,iz_cple
    read(luq_bqf) z_depth
!
    if(iq_prt >=2) then
      write(luq_prt,'(a)')    'Q_CTRGRID: Contents of BQF file'
      write(luq_prt,'(2a)')   'Q_CTRGRID: Header:',trim(r_header)
      write(luq_prt,'(a,i4)') 'Q_CTRGRID: NK:',nkz
      write(luq_prt,'(a,i4)') 'Q_CTRGRID: NA:',naz
    end if
  end if
!---------------------------------------------------------------------------------------
! check spectral interaction grid and depth for consistency
!---------------------------------------------------------------------------------------
  if(.not. lq_grid) then
az: do iaz = 1,naz
      if(abs(q_ad(iaz)-z_ad(iaz)) > 0.01) then
        write(luq_prt,'(a)') 'Q_CTRGRID: Directions do not agree'
        do jaz=1,naz
          write(luq_prt,'(1x,a,i4,2f10.3)') 'iaz q_ad z_ad:',jaz,q_ad(jaz),z_ad(jaz)
        end do
        lq_grid = .true.
        igrid = 2
        exit az
      end if
    end do az
  end if
!
  if(.not. lq_grid) then
ak: do ikz = 1,nkz
      if(abs(q_sig(ikz)-z_sig(ikz)) > 0.01) then
        write(luq_prt,'(a)') 'Q_CTRGRID: Wave numbers do not agree'
        do jkz=1,nkz
          write(luq_prt,'(1x,a,i4,2f10.3)') 'ikz q_k z_sig:',jkz,q_sig(jkz),z_sig(jkz)
        end do
        lq_grid = .true.
        igrid = 2
        exit ak
      end if
    end do ak
  end if
!
!  compare water depths
!
  if(abs(z_depth-s_depth) > 0.09 .and. iq_disp > 1 .and. .not. lq_grid) then
    write(luq_prt,'(a)') 'Q_CTRGRID: Water depths do not agree'
    write(luq_prt,'(a,2f16.2)') 'Q_CTRGRID: q_depth z_depth:',q_depth,z_depth
    lq_grid = .true.
    igrid = 2
  end if
!
  if(lq_grid) then
    close(luq_bqf)
    if(iq_log >= 1) write(luq_log,'(a)') 'Q_CTRGRID: Existing BQF-file invalid, it will be closed'
  end if
!
else
  lq_grid = .true.
  igrid = 1
end if
!------------------------------------------------------------------------------
if(itask==1) then
  if(luq_bqf>0) call z_fclose(luq_bqf)
  goto 9999
end if
!-----------------------------------------------------------------------------
!  if lq_grid==true  a new grid has to be generated
!                    if not, read the grid information into memory
!  or iq_make==2     always make an interaction grid
!  or iq_make==3     as 2, plus stop after making grid
!
if(lq_grid .or. iq_make==2 .or. iq_make==3) then
!
  if(luq_bqf>0) call z_fclose(luq_bqf)
  call z_fileio(bqname,'UU',iufind,luq_bqf,iuerr)   ! binary quadruplet file
!
  if(iq_log >= 1) then
    write(luq_log,*)
    write(luq_log,'(a)') 'Q_CTRGRID: New grid will be generated'
    write(luq_log,'(a,a)') 'Q_CTRGRID: Name of BQF file:',trim(bqname)
    write(luq_log,'(a,i4)') 'Q_CTRGRID: '//trim(bqname)//' connected to :',luq_bqf
  end if
!
  if(iq_screen >= 1) write(iscreen,'(2a)') &
& 'Q_CTRGRID: Generating wave number grid for quadruplet interactions: ',trim(bqname)
!
  q_depth = s_depth
  call q_makegrid
  q_depth = q_depth_saved
!
  if(iq_err /=0) then
     lastquadfile = 'quad_err_.bqf'
     goto 9999
  end if
!
  igrid = 0
!
  close(luq_bqf)
!
  if(iq_log >=1) then
    write(luq_log,'(a,i4)') 'Q_CTRGRID: '//trim(bqname)//' disconnected from:',luq_bqf
  end if
!
  if(iq_screen >=1) write(iscreen,'(a)') 'Q_CTRGRID: Grid generation completed succesfully'
!----------------------------------------------------------------------------------------
!
!  check of header and spectral grid succesfull
!  such that data can be read from BQF file
!----------------------------------------------------------------------------------------
!
else
  if(iq_screen >= 1) write(iscreen,'(2a)') 'Q_CTRGRID: Reading existing grid: ',trim(bqname)
  if(iq_prt >= 1)    write(luq_prt,'(2a)')  'Q_CTRGRID: Existing grid will be read:',trim(bqname)
  if(iq_log >= 1)    write(luq_log,'(2a)')  'Q_CTRGRID: Existing grid will be read:',trim(bqname)
!
  read(luq_bqf) quad_nloc
  read(luq_bqf) quad_ik2
  read(luq_bqf) quad_ia2
  read(luq_bqf) quad_ik4
  read(luq_bqf) quad_ia4
  read(luq_bqf) quad_w1k2
  read(luq_bqf) quad_w2k2
  read(luq_bqf) quad_w3k2
  read(luq_bqf) quad_w4k2
  read(luq_bqf) quad_w1k4
  read(luq_bqf) quad_w2k4
  read(luq_bqf) quad_w3k4
  read(luq_bqf) quad_w4k4
  read(luq_bqf) quad_zz
!
  lastquadfile = bqname
!
  close(luq_bqf)
!

end if
!
9999 continue
!
if (allocated(z_ad)) deallocate(z_ad,z_sig)
!
!
!!!!!!!!!!!!call q_stack('-q_ctrgrid')
!
return
end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!--------------------------------!!!!!
!!-------------Q-MAKEGRID---------------!!
!!!!!--------------------------------!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine q_makegrid
!------------------------------------------------------------------------------
!
!   +-------+    ALKYON Hydraulic Consultancy & Research
!   |       |    Gerbrant van Vledder
!   |   +---+
!   |   | +---+  Last update: 10 June 2003
!   +---+ |   |  Release: 5.0
!         +---+
!
! do not use m_xnldata
use m_constants
use serv_xnl4v5
!
!  0. Update history
!
!     25/02/1999  Initial version
!     11/10/1999  Error handling improved; Bugs fixed when w1=w3
!     12/10/1999  Storage modified and non-geometric option included
!     16/10/1999  Equation for computing address of 2d array simplified
!     21/10/1999  Range of precomputed grid added to data file
!     22/10/1999  Renaming of some indices
!     25/10/1999  Header with grid info extended
!     12/11/1999  Output format modified of data to GRD file, adapted
!                 for use on UNIX systems at WES
!     08/12/1999  Interface with A_CMPLOC extended
!     28/12/1999  Routine A_CMPLOC renamed to Q_CMPLOC
!     03/01/2000  IQ_START replaced by IQ_LOCUS
!     05/01/2000  Interface with Q_CMPLOC modified
!     08/02/2000  Output to LUQLOC made conditional
!     09/08/2002  Name changed from Q_GRIDV1 to Q_MAKEGRID
!                 Upgrade to release 4.0
!     15/08/2002  Bug fixed in indexing bins below lowest wave number
!     20/08/2002  Sigma written to QUAD file, instead of wave numbers
!     22/08/2002  Data along locus compacted, elimate zero's
!     10/09/2002  Upgrade to release 5
!                 Value of LASTQUADFILE set
!     10/06/2003  Output to GRD file always without compacting
!
!  1. Purpose:
!
!     Set-up grid for computation of loci
!
!     Generate data file with basic loci for computation of
!     nonlinear quadruplet interactions
!
!  2. Method
!
!
!  3. Parameter list:
!
!     Name    I/O  Type  Description
!
!  4. Error messages
!
!  5. Called by:
!
!     Q_CTRGRID
!
!  6. Subroutines used
!
!     Q_STACK
!     Q_CPMLOCUS
!     Q_MODIFY
!     Q_WEIGHT
!     Q_CHKRES
!     Q_NEAREST
!
!  7. Remarks
!
!  8. Structure
!
!  9. Switches
!
! 10. Source code
!------------------------------------------------------------------------------
!     Local variables
!
integer iloc,jloc            ! counters
integer iaq,ikq              ! counters
integer iaq3,ikq1,ikq3,nkq1  ! counters
integer jaq1,jaq3            ! counters
integer amem,kmem            ! index of angle and wave number in grid
real aa1,aa3,kk1,kk3         ! temporary wave number variables
!
integer nzloc                ! counter for non-zero contributions along locus
integer nztot1,nztot2        ! total number of zero and non-zero points on locus
integer ik2,ia2              ! index of wave number k2
integer ik4,ia4              ! index of wave number k4
!
real wk,wa                   ! weights
real w1k2,w2k2,w3k2,w4k2     ! interpolation weights
real w1k4,w2k4,w3k4,w4k4     ! interpolation weights
!
real ka,kb       ! lower and higher wave number magnitude
real km          ! wave number at mid point
real kw          ! half width of locus
!
real tfac        ! combined tail factor
!
logical lwrite   ! indicator if binary interaction grid has been written successfully
real smax        ! maximum s-value
!
real, allocatable :: xloc(:),yloc(:)
real qq
!-------------------------------------------------------------------------------
!!!!!!!!!!!!!!!!call q_stack('+q_makegrid')
!
! initializations
!
lwrite  = .false.
nztot1  = 0
nztot2  = 0
!%
quad_nloc = -1    ! number of points on all loci
!%
if(allocated(xloc)) deallocate(xloc) ; allocate (xloc(mlocus))
if(allocated(yloc)) deallocate(yloc) ; allocate (yloc(mlocus))
!
!  write header to grid file
!
!
!  set range of do loops for computing interaction grid
!
if(iq_geom==0 .or. iq_disp/=1) then
  nkq1 = nkq           ! loop over all k1 wave numbers, since no geometric scaling can be used
else
  nkq1 = 1             ! use only first wave number for k1, since geometric scaling can be used
end if
!
jaq1 = 1               ! index of direction of k1 in grid matrix
!-------------------------------------------------------------------------------------
!  compute components of reference wave number,
!  for setting up interaction grid
!-------------------------------------------------------------------------------------
k1: do ikq1=1,nkq1
!
  if(iq_screen==2) write(iscreen,*) 'k1-ring:',ikq1
!
  aa1   = q_ad(iaref)
  kk1   = q_k(ikq1)
  krefx = kk1*cos(q_ad(iaref)*dera)
  krefy = kk1*sin(q_ad(iaref)*dera)
!
  k1x  = krefx
  k1y  = krefy
!

k3: do ikq3 = ikq1,nkq   !
   if(iq_screen==2) write(iscreen,*) 'k1-k3 indices:',ikq1,ikq3
!
    kk3 = q_k(ikq3)
!
!
a3: do iaq3 = iag1,iag2
!
      if(iaq3 == iag1 .and. ikq3 == ikq1) cycle
!
      aa3 = q_ad(iaq3)
      k3x = kk3*cos(aa3*dera)
      k3y = kk3*sin(aa3*dera)
!------------------------------------------------------------------------------
!   compute locus for a specified combination of k1 and k3
!
!-----------------------------------------------------------------------------
      ia_k1 = iaq1; ik_k1 = ikq1
      ia_k3 = iaq3; ik_k3 = ikq3
      call q_cmplocus(ka,kb,km,kw,crf1)
!
      if(iq_err/=0) goto 9999
!------------------------------------------------------------------------------
!     redistibute or filter data points along locus
!
      call q_modify
      if(iq_err > 0) goto 9999
!------------------------------------------------------------------------------
!     compute weights for interpolation in computational grid
!
      call q_weight
      if(iq_err > 0) goto 9999
!------------------------------------------------------------------------------
!    special storing mechanism for interactions per combination of k1 and k3
!
      kmem  = (ikq3-ikq1+1) - (ikq1-2*nkq-2)*(ikq1-1)/2;
      jaq3  = iaq3-iaref+1        ! ensure that data stored in matrix start at index (1,1)
      amem  = jaq3                ! index of direction
!
!
!-------------------------------------------------------------------------------
!     Convert real indices to integer indexing and real weights
!
!    3-----------4 ja2p         w1 = (1-wk)*(1-wa)
!    |    .      |              w2 = wk*(1-wa)
!    |. . + . . .| wa2   A      w3 = (1-wk)*wa
!    |    .      |       |      w4 = wk*wa
!    |    .      |       wa
!    |    .      |       |
!    1-----------2 ja2   V
!   jk2  wk2  jk2p
!
!    <-wk->
!
!-------------------------------------------------------------------------------
      nzloc = 0
!
loc:  do iloc = 1,nlocus
!
        ik2  = floor(wk_k2(iloc))
        ia2  = floor(wa_k2(iloc))
        wk   = wk_k2(iloc)-real(ik2)
        wa   = wa_k2(iloc)-real(ia2)
        w1k2 = (1.-wk)*(1.-wa)
        w2k2 = wk*(1.-wa)
        w3k2 = (1.-wk)*wa
        w4k2 = wk*wa
!
        ik4  = floor(wk_k4(iloc))
        ia4  = floor(wa_k4(iloc))
        wk   = wk_k4(iloc)-real(ik4)
        wa   = wa_k4(iloc)-real(ia4)
        w1k4 = (1.-wk)*(1.-wa)
        w2k4 = wk*(1.-wa)
        w3k4 = (1.-wk)*wa
        w4k4 = wk*wa
!
!  Take care of points that lie below lowest wave number
!  when no geometric scaling is applied, then modify weights
!  such that directional position is retained
!
        if(iq_geom==0) then
          if(ik2 ==0) then
            ik2  = 1
            w1k2 = w1k2 + w2k2
            w2k2 = 0.
            w3k2 = w3k2 + w4k2
            w4k2 = 0.
          end if
          if(ik4 ==0) then
            ik4  = 1
            w1k4 = w1k4 + w2k4
            w2k4 = 0.
            w3k4 = w3k4 + w4k4
            w4k4 = 0.
          end if
        end if
!
!  compute combined tail factor and product of coupling coefficient, step size,
!  symmetry factor, and tail factor divided by jacobian
!
        tfac = wt_k2(iloc)*wt_k4(iloc)
        quad_zz(kmem,amem,iloc)   = cple_mod(iloc)*ds_mod(iloc)*sym_mod(iloc)/jac_mod(iloc)*tfac
!
!----------------------------------------------------------------------------------------
!  compact data by elimating zero-contribution on locus
!----------------------------------------------------------------------------------------
!
        if(iq_compact==1 .and. abs(quad_zz(kmem,amem,iloc)) > 1.e-15) then
          nzloc = nzloc + 1
          jloc  = nzloc
          nztot1 = nztot1 + 1
        else
          jloc = iloc
        end if
        nztot2 = nztot2 + 1
!
!  shift data
!
        quad_zz(kmem,amem,jloc)  = quad_zz(kmem,amem,iloc)
!
        quad_ik2(kmem,amem,jloc) = ik2           ! lower wave number index of k2
        quad_ia2(kmem,amem,jloc) = ia2           ! lower direction index of k2
        quad_ik4(kmem,amem,jloc) = ik4           ! lower wave number index of k4
        quad_ia4(kmem,amem,jloc) = ia4           ! lower direction index of k4
!
        quad_w1k2(kmem,amem,jloc) = w1k2         ! weight 1 of k2
        quad_w2k2(kmem,amem,jloc) = w2k2         ! weight 2 of k2
        quad_w3k2(kmem,amem,jloc) = w3k2         ! weight 3 of k2
        quad_w4k2(kmem,amem,jloc) = w4k2         ! weight 4 of k2
!
        quad_w1k4(kmem,amem,jloc) = w1k4         ! weight 1 of k4
        quad_w2k4(kmem,amem,jloc) = w2k4         ! weight 2 of k4
        quad_w3k4(kmem,amem,jloc) = w3k4         ! weight 3 of k4
        quad_w4k4(kmem,amem,jloc) = w4k4         ! weight 4 of k4
!
!
      end do loc
!
      if(iq_compact==1) then
        quad_nloc(kmem,amem) = nzloc                ! store compacted number of points on locus
      else
        quad_nloc(kmem,amem) = nlocus               ! store number of points on locus
        nzloc = nlocus
      end if
!
!     write(luq_prt,'(a,4i5)') 'Q_MAKEGRID kmem amem nlocus:',kmem,amem,nlocus,nzloc
!
    end do a3
  end do k3
end do k1
!------------------------------------------------------------------------------
!  Write locus information to binary file
!------------------------------------------------------------------------------
!
write(luq_bqf) q_header
!
!------------------------------------------------------------------------------
! spectral interaction grid
!------------------------------------------------------------------------------
!
write(luq_bqf) naq,nkq
write(luq_bqf) q_sig
write(luq_bqf) q_ad
write(luq_bqf) iq_geom,iq_disp,iq_geom
write(luq_bqf) q_depth
!
!------------------------------------------------------------------------------
! interaction grid
!------------------------------------------------------------------------------
!
write(luq_bqf) quad_nloc
write(luq_bqf) quad_ik2
write(luq_bqf) quad_ia2
write(luq_bqf) quad_ik4
write(luq_bqf) quad_ia4
write(luq_bqf) quad_w1k2
write(luq_bqf) quad_w2k2
write(luq_bqf) quad_w3k2
write(luq_bqf) quad_w4k2
write(luq_bqf) quad_w1k4
write(luq_bqf) quad_w2k4
write(luq_bqf) quad_w3k4
write(luq_bqf) quad_w4k4
write(luq_bqf) quad_zz
!
!
lwrite = .true.
lastquadfile = bqname
!
if(iq_screen >= 1 .and. iq_test>=1) write(iscreen,'(2a)') 'q_makegrid: LASTQUADFILE: ',lastquadfile
!
9999 continue
!
if(allocated(xloc)) deallocate(xloc,yloc)
!
! check if BQF file has been written succesfully
! if not, deleted both the AQFILE and BQFILE
!
if(.not. lwrite) then
  close(luq_bqf,status='delete')
  if(iq_log > 0) then
    write(luq_log,*)
    write(luq_log,*) 'Q_MAKEGRID: Grid files ',trim(aqname),' and ',trim(bqname),' deleted'
    write(luq_log,*) 'Q_MAKEGRID: Since an error occurred during the generation'
    write(luq_log,*) 'Q_MAKEGRID: of the interaction grid'
  end if
end if
!-------------------------------------------------------------------------------
!  write statistics of compacting to print file
!
if(iq_prt >=1) then
  if(iq_compact==0) nztot1 = nztot2
  write(luq_prt,'(a,i10)') 'Total number of points on loci        :',nztot2
  write(luq_prt,'(a,i10)') 'Total number of stored points on locus:',nztot1
  write(luq_prt,'(a,i10)') 'Total number of zero points on locus  :',nztot2-nztot1
  write(luq_prt,'(a,f8.2)') 'Reduction factor (%):',real(nztot2-nztot1)/real(nztot2)*100.
end if
!
!!!!!!!!!!!!call q_stack('-q_makegrid')
!
return
end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!--------------------------------!!!!!
!!-------------Q-CMPLOCUS---------------!!
!!!!!--------------------------------!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


subroutine q_cmplocus(ka,kb,km,kw,loclen)
!------------------------------------------------------------------------------
!
!   +-------+    ALKYON Hydraulic Consultancy & Research
!   |       |    Gerbrant van Vledder
!   |   +---+
!   |   | +---+  Last update: 8 August 2003
!   +---+ |   |  Release: 5.0
!         +---+
!
! do not use m_xnldata
use m_constants
use serv_xnl4v5
implicit none
!-------------------------------------------------------------------------------
!  0. Update history
!
!     Date        Description
!
!     18/11/1999  Initial version
!     08/12/1999  Tracing of locus updated and Q_TRACE included
!     22/12/1999  Option Q_POLAR included
!     05/01/2000  LOCLEN added in interface
!     09/08/2002  Upgrade to release 4.0
!     10/09/2002  g added to interface with X_CPLE
!                 test output modified
!     12/06/2003  Call to Z_POYAREA added to check POLAR2
!     08/08/2003  Check on areas only for loci with k3m/k1m < 100
!                 Otherwise machine accuracy plays a role
!
!  1. Purpose:
!
!     Compute locus function used for the determination of the
!     resonnance condition
!
!  2. Method
!
!     See ALKYON, 1999
!
!  3. Parameter list:
!
!Type   I/O          Name          Description
!----------!----------------------------------------------------------------------------
real, intent(out) :: ka,kb       ! lowest and highest wave number magnitude of k2-locus
real, intent(out) :: km          ! wave number magnitude at mid point
real, intent(out) :: kw          ! half width of locus
real, intent(out) :: loclen      ! estimated length of locus
!
!  4. Error messages
!
!  5. Called by:
!
!    Q_MAKEGRID
!
!  6. Subroutines used
!
!     z_zero1
!
!
!
!  7. Remarks
!
!  8. Structure
!
!  9. Switches
!
! 10. Source code
!-------------------------------------------------------------------------------
!     Local variables
!-------------------------------------------------------------------------------
real k1m          ! magnitude of wave number k1
real k3m          ! magnitude of wave number k3
real pcos,psin    ! cosine and sine of normalize angle of P
real klen         ! total length of line locus for case w1=w3
!
real kx_beg       ! x-component at start point
real ky_beg       ! y-component at start point
real kx_end       ! x-component at end point
real ky_end       ! y-component at end point
!
real dsp,dsm      ! distances in plus and minus direction
real sum          ! total length of locus
!
real w1,w3        ! radian frequencies of wave numbers k1 and k3
real eps          ! local accuracy for determination of roots
real area1        ! area of locus as computed
real area2        ! area of locus as based on LOCPOS and ellipse
real ratio        ! maximum ratio between area1 and area2
!
integer ierr      ! local error level
integer iloc,jloc ! counters along locus
integer itest     ! local test level for test output to unit luqtst
integer ip1       ! index +1
integer im1       ! index -1
integer jj        ! counter
!------------------------------------------------------------------------------
!  function declarations
!
!!real x_disper                 ! dispersion relation
!!real x_cple                   ! coupling coefficient
!!real x_jacobian               ! Jacobian term
!------------------------------------------------------------------------------
!!!!!!!!!!!!!!!!call q_stack('+q_cmplocus')
!
!  set initial values
!
eps     = 10.*epsilon(1.)     ! set accuracy 10 times machine accuracy
itest   = iq_test             ! assign test level from overall setting
!! itest   = 1                ! (re)set local test level
!
! compute characteristics of configuration
!
px   = k1x - k3x
py   = k1y - k3y
pmag = sqrt(px**2 + py**2)
xang = atan2(-px,py)
pang = atan2(py,px)
k1m  = sqrt(k1x**2 + k1y**2)
k3m  = sqrt(k3x**2 + k3y**2)
w1   = x_disper(k1m,q_depth)
w3   = x_disper(k3m,q_depth)
q    = w1-w3
!
!  compute cosine and sine of direction of P-vector
!  reverse direction for the case q<0
!
if(q < 0) then
  sang = pang+pi
  pcos = cos(pang+pi)
  psin = sin(pang+pi)
else
  sang = pang
  pcos = cos(pang)
  psin = sin(pang)
end if
!
!
!  first solution along locus: k2 = k3
!
!  check for special case if q = 0
!
if (abs(q) < eps_q) then
!
  call q_loc_w1w3(k1x,k1y,k3x,k3y,nlocus0,x2_loc,y2_loc,x4_loc,y4_loc,s_loc)
  nlocus1 = nlocus0
  ds_loc  = s_loc(2)-s_loc(1)
  klen    = s_loc(nlocus0)
  ka      = 0.
  kb      = 0.
  km      = 0.
  kw      = 0.
  sang    = xang
!
else
!------------------------------------------------------------------------------
!  compute characteristics of locus, such as its position in
!  wave number space
!------------------------------------------------------------------------------
!
  call q_locpos(ka,kb,km,kw,loclen)
  if(iq_err/=0) goto 9999
!
!  compute position of start and end point for tracing
!  the locus
!
  kx_beg = ka*pcos
  ky_beg = ka*psin
  kx_end = kb*pcos
  ky_end = kb*psin
!
!  compute position of locus using polar method
!  see Van Vledder (2000)
!
!%  call q_polar(ka,kb,kx_beg,ky_beg,kx_end,ky_end,loclen,ierr)
  call q_polar2(ka,kb,kx_beg,ky_beg,kx_end,ky_end,loclen,ierr)
!
! check area of locus by a simple test  (added 12 June 2003)
!
  call z_polyarea(x2_loc,y2_loc,nlocus1,area1)
  area2 = pi*(kb-ka)*0.5*kw
  ratio = max(area1/area2,area2/area1)
!
!
  if(ratio>1.5 .and. k3m/k1m < 100.) then
    !!!!!!!!call q_error('e','LOCUS','Severe problem in POLAR2')
    print *, "Error 1234"
    write(luq_err,'(a)') 'Q_CMPLOCUS: ratio > 1.5'
!
    goto 9999
  end if
!
!  01/10/2001
!  compute position of k4 locus by a simple translation
!
  do iloc=1,nlocus1
    x4_loc(iloc) = x2_loc(iloc) + px
    y4_loc(iloc) = y2_loc(iloc) + py
  end do
!
end if
!
if (iq_test >=2) write(luq_tst,'(1x,a,4f12.5,i4)')&
&   'Q_CMPLOCUS: k1x/y k3x/y nlocus:',k1x,k1y,k3x,k3y,nlocus1
!----------------------------------------------------------------------------------
! compute characteristics around locus
!----------------------------------------------------------------------------------
!
s_loc(1) = 0.
sum      = 0
!
do iloc=1,nlocus1
!
!  compute step sizes
!
  if (abs(q) < eps_q) then
!
!  for this case the sum of ds_loc is unequal to s_loc(nlocus1)
!
    sum = s_loc(nlocus1)
  else
!
!   compute indices of previous and next index on locus
!
    ip1 = iloc+1
    if (ip1 > nlocus1) ip1 = 1
    im1 = iloc-1
    if (im1 < 1) im1 = nlocus1
!
    dsp = sqrt((x2_loc(iloc)-x2_loc(ip1))**2 + (y2_loc(iloc)-y2_loc(ip1))**2)
    dsm = sqrt((x2_loc(iloc)-x2_loc(im1))**2 + (y2_loc(iloc)-y2_loc(im1))**2)
    if(iloc < nlocus1) s_loc(iloc + 1) = s_loc(iloc) + dsp
    ds_loc(iloc) = 0.5*(dsp + dsm)
    sum = sum+ds_loc(iloc)
  end if
!
!  compute gradient/Jacobian terms along locus
!
   jac_loc(iloc) = x_jacobian(x2_loc(iloc),y2_loc(iloc),x4_loc(iloc),y4_loc(iloc))
!
!  compute coupling coefficients along locus
!
  k2x = x2_loc(iloc)
  k2y = y2_loc(iloc)
  k4x = x4_loc(iloc)
  k4y = y4_loc(iloc)
!
  cple_loc(iloc) = x_cple(k1x,k1y,k2x,k2y,k3x,k3y,k4x,k4y,iq_cple,q_depth,q_grav)
!
end do
!
!
9999 continue
!
!!!!!!!!!!!!!!!!!!!!call q_stack('-q_cmplocus')
!
return
end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!--------------------------------!!!!!
!!-------------Q-MODIFY-----------------!!
!!!!!--------------------------------!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


subroutine q_modify
!------------------------------------------------------------------------------
!
!   +-------+    ALKYON Hydraulic Consultancy & Research
!   |       |    Gerbrant van Vledder
!   |   +---+
!   |   | +---+  Last update: 11 June 2003
!   +---+ |   |  Release: 5.0
!         +---+
!
! do not use m_xnldata
use m_constants
use serv_xnl4v5
implicit none
!--------------------------------------------------------------------------------
!  0. Update history
!
!       9/04/1999  Initial version
!      13/04/1999  New intermediate variables *_mod introduced
!      11/10/1999  Check on error messages in interpolation added
!      18/10/1999  Bug fixed in assigning new ds values to array DS_MOD
!      27/10/1999  Checked added on allocated of SOLD
!       8/12/1999  Test output added
!      29/12/1999  Bug fixed in assigning DS_MOD for first and last point on locus
!       1/10/2001  Components of k4-locus added
!                  No interpolation and modification if q==0
!       9/08/2002  Upgrade to version 4.0
!      15/08/2002  Step sizing improved
!       4/06/2003  Bug fixed in computing slen (length of locus)
!                  Locus closed to enable interpolation to finer resolution
!       6/06/2003  Activate output to XDIA configuration file
!      10/06/2003  Conversion to new indexing and lumping debugged
!      11/06/2003  Call to subroutine Q_SYMMETRY added
!
!  1. Purpose:
!
!     Modify points along the locus, such that they are evenly distributed
!     Only when intented, i.e. when IQ_LOCUS==2
!
!  2. Method
!
!     Compute new spacing along locus
!     Redistribute points and coefficient at new spacing using linear interpolation
!     Output DIA configuration when also lumping active
!
!     If no redistribution is needed, then copy relevant data
!
!  3. Parameter list:
!
!     Name    I/O  Type  Description
!
!  4. Error messages
!
!  5. Called by:
!
!     Q_CMPLOCUS
!
!  6. Subroutines used
!
!     Q_STACK
!     Q_SYMMETRY
!     Z_INTP1
!
!  7. Remarks
!
!  8. structure
!
!  9. Switches
!
! 10. Source code
!------------------------------------------------------------------------------
!     Local parameters
!
integer ierr,jerr    ! error indicators
integer nold,nnew    ! old and new number of points on locus
integer iold,inew    ! counter for loop along points
integer iloc         ! counter for loop along locus
integer jloc         ! counter for loop over lumped locus
integer itest        ! local test level, by default equal to IQ_TEST
!
real k2a,k2m         ! angle (deg) and wave number magnitude of wave number k2
real k4a,k4m         ! angle (deg) and wave number magnitude of wave number k4
real w2,w4           ! radian frequencies of wave numbers
!
!
real dk13,dk14       ! difference wave number
real dsnew,slen      ! new step size and length of locus
real zero            ! 0
real q_eps           ! accuracy to distinguish special case, with q=0
real diold           ! 'real' old number of indices between succeeding lumped bins
real dinew           ! 'real' new number of indices between succeeding lumped bins
!
!!real x_disper        ! evaluate dispersion relation
real, allocatable :: sold(:)     ! old coordinate along locus
real, allocatable :: snew(:)     ! new coordinate along locus
!--------------------------------------------------------------------------
!!!!!!!!!!!!!!!!!call q_stack('+q_modify')
!
!  initialisations
!
zero  = 0.
q_eps = 1.e-5
itest = iq_test
!
! itest = 1   ! set local test level for test purposes
!
if(itest>=1) then
  write(luq_tst,'(a,i4)') 'Q_MODIFY: iq_mod   :',iq_mod
  write(luq_tst,'(a,i4)') 'Q_MODIFY: iq_xdia  :',iq_xdia
  write(luq_tst,'(a,i4)') 'Q_MODIFY: iq_lump  :',iq_lump
  write(luq_tst,'(a,i4)') 'Q_MODIFY: iq_gauleg:',iq_gauleg
end if
!------------------------------------------------------------------------------
!  do not modify data when IQ_MOD==0
!------------------------------------------------------------------------------
!
if(iq_mod==0) then
  nlocus   = nlocus1
  x2_mod   = x2_loc
  y2_mod   = y2_loc
  x4_mod   = x4_loc
  y4_mod   = y4_loc
  s_mod    = s_loc
  ds_mod   = ds_loc
  jac_mod  = jac_loc
  cple_mod = cple_loc
  call q_symmetry(k1x,k1y,k3x,k3y,x4_mod,y4_mod,sym_mod,nlocus)
else
!------------------------------------------------------------------------------
! Modify spacing along locus
!------------------------------------------------------------------------------
  nold = nlocus1
!
! close locus by adding one point, equal to first point
! only for normal locus
!
  if(abs(q)>q_eps) nold  = nold+1
!
!------------------------------------------------------------------------------
! Determine new number of points along locus
!------------------------------------------------------------------------------
!
  if(iq_gauleg > 0) then
    nnew = iq_gauleg
  elseif(iq_lump > 0) then
    nnew = iq_lump
  else
    nnew  = nlocus0
  end if
!
!
  allocate (sold(nold),snew(nnew))
!------------------------------------------------------------------------------
!  Compute circumference of locus, distinguish 2 case, open or closed
!------------------------------------------------------------------------------
!
  if(abs(q)<q_eps) then
    slen = s_loc(nold)
    sold = s_loc
  else
    slen = 0
    do iold=1,nold-1               ! loop length minus one, since locus is closed
      sold(iold) = s_loc(iold)
      slen = slen + ds_loc(iold)
    end do
!
!------------------------------------------------------------------------------
!  close locus by copying first value in last value
!------------------------------------------------------------------------------
!
    sold(nold)     = slen
    x2_loc(nold)   = x2_loc(1)
    y2_loc(nold)   = y2_loc(1)
    x4_loc(nold)   = x4_loc(1)
    y4_loc(nold)   = y4_loc(1)
    jac_loc(nold)  = jac_loc(1)
    cple_loc(nold) = cple_loc(1)
  end if
!
!------------------------------------------------------------------------------
! compute new spacing along loci and coordinates along locus
! Gauss-Legendre integration
!------------------------------------------------------------------------------
!
  if(iq_gauleg > 0) then
    if(iq_gauleg > nnew) stop 'Q_MODIFY: iq_gauleg > nlocus0'
    nnew = iq_gauleg
    call y_gauleg(zero,slen,snew,ds_mod,nnew)
!
  else
    if(abs(q)>q_eps) then
      dsnew  = slen/real(nnew)
      do inew=1,nnew
        snew(inew) = (inew-1.)*dsnew
      end do
    else
      dsnew  = slen/real(nnew-1.)
      do inew=1,nnew
        snew(inew) = (inew-1)*dsnew
      end do
    end if
    ds_mod = dsnew
  end if
!
!
  jerr = 0
!------------------------------------------------------------------------------
! Compute characteristics of locus for special case q=0
!------------------------------------------------------------------------------
!
  if(abs(q)<1.e-5) then
    call z_intp1(sold,x2_loc,snew,x2_mod,nold,nnew,ierr)
    if(ierr > 0) write(luq_err,*) 'Z_INTP1 x_loc, ierr=',ierr
    jerr = jerr + ierr
!
    call z_intp1(sold,y2_loc,snew,y2_mod,nold,nnew,ierr)
    if(ierr > 0) write(luq_err,*) 'Z_INTP1 x_loc, ierr=',ierr
    jerr = jerr + ierr
 !
    call z_intp1(sold,x4_loc,snew,x4_mod,nold,nnew,ierr)
    if(ierr > 0) write(luq_err,*) 'Z_INTP1 x_loc, ierr=',ierr
    jerr = jerr + ierr
!
    call z_intp1(sold,y4_loc,snew,y4_mod,nold,nnew,ierr)
    if(ierr > 0) write(luq_err,*) 'Z_INTP1 x_loc, ierr=',ierr
    jerr = jerr + ierr
!
    call z_intp1(sold,s_loc,snew,s_mod,nold,nnew,ierr)
    if(ierr > 0) write(luq_err,*) 'Z_INTP1 x_loc, ierr=',ierr
    jerr = jerr + ierr
!
    call q_symmetry(k1x,k1y,k3x,k3y,x4_mod,y4_mod,sym_mod,nnew)
!
! ---  lumping along locus --------------------------------------------
!
    if(iq_lump>0) then
      diold  = slen/real(nold)
      dinew  = slen/real(nnew)
      ds_mod = 0.
      call q_symmetry(k1x,k1y,k3x,k3y,x4_loc,y4_loc,sym_loc,nold)
!
      do iloc=1,nlocus1
        jloc = floor((iloc-1.)*diold/dinew)+1
        ds_mod(jloc)   = ds_mod(jloc) + cple_loc(iloc)*ds_loc(iloc)/jac_loc(iloc)*sym_loc(iloc)
        jac_mod(jloc)  = 1.
        cple_mod(jloc) = 1.
      end do
!
      sym_mod = 1          ! symmetry already taken account in lumping proces
!
! --- No lumping -------------------------------------------------------------
!
    else
      call z_intp1(sold,jac_loc,snew,jac_mod,nold,nnew,ierr)
      if(ierr > 0) write(luq_err,*) 'Z_INTP1 jac_loc, ierr=',ierr
      jerr = jerr + ierr
!
      call z_intp1(sold,cple_loc,snew,cple_mod,nold,nnew,ierr)
      if(ierr > 0) write(luq_err,*) 'Z_INTP1 cp_loc, ierr=',ierr
      jerr = jerr + ierr
    end if
!------------------------------------------------------------------------------------------------
!  compute characteristics for closed locus
!------------------------------------------------------------------------------------------------
  else
    call z_intp1(sold,x2_loc,snew,x2_mod,nold,nnew,ierr)
    if(ierr > 0) write(luq_err,*) 'Z_INTP1 x_loc, ierr=',ierr
    jerr = jerr + ierr
!
    call z_intp1(sold,y2_loc,snew,y2_mod,nold,nnew,ierr)
    if(ierr > 0) write(luq_err,*) 'Z_INTP1 y_loc, ierr=',ierr
    jerr = jerr + ierr
!
    call z_intp1(sold,x4_loc,snew,x4_mod,nold,nnew,ierr)
    if(ierr > 0) write(luq_err,*) 'Z_INTP1 x_loc, ierr=',ierr
    jerr = jerr + ierr
!
    call z_intp1(sold,y4_loc,snew,y4_mod,nold,nnew,ierr)
    if(ierr > 0) write(luq_err,*) 'Z_INTP1 y_loc, ierr=',ierr
    jerr = jerr + ierr
!
    call z_intp1(sold,s_loc,snew,s_mod,nold,nnew,ierr)
    if(ierr > 0) write(luq_err,*) 'Z_INTP1 s_loc, ierr=',ierr
    jerr = jerr + ierr
!
!
    call q_symmetry(k1x,k1y,k3x,k3y,x4_mod,y4_mod,sym_mod,nnew)
!
!  ----- Lumping along locus -----------------------------------
!
    if(iq_lump>0) then
      diold  = slen/real(nold-1)
      dinew  = slen/real(nnew)
      ds_mod = 0.
      call q_symmetry(k1x,k1y,k3x,k3y,x4_loc,y4_loc,sym_loc,nold)
!
      do iloc=1,nold-1
        jloc = floor((iloc-1.)*diold/dinew + 1.49999)
        jloc = mod(jloc-1+nnew,nnew)+1
        ds_mod(jloc)   = ds_mod(jloc) + cple_loc(iloc)*ds_loc(iloc)/jac_loc(iloc)*sym_loc(iloc)
        jac_mod(jloc)  = 1.
        cple_mod(jloc) = 1.
      end do
!
      sym_mod = 1          ! symmetry already taken account in lumping proces
!
!------------  No lumping along locus  --------------------------------
!
    else
      call z_intp1(sold,jac_loc,snew,jac_mod,nold,nnew,ierr)
      if(ierr > 0) write(luq_err,*) 'Z_INTP1 jac_loc, ierr=',ierr
      jerr = jerr + ierr
!
      call z_intp1(sold,cple_loc,snew,cple_mod,nold,nnew,ierr)
      if(ierr > 0) write(luq_err,*) 'Z_INTP1 cp_loc, ierr=',ierr
      jerr = jerr + ierr
    end if
!
    if(jerr > 0) then
      iq_err = iq_err + 1
      !!!!!!!!!!!!!call q_error('e','INTER','Problem in interpolation process')
      print *, "Error 21732"
      goto 9999
    end if
  end if
!
  nlocus = nnew
!
end if
!
!------------------------------------------------------------------------------
!
!
!------------------------------------------------------------------------------
!
!!  compute symmetry factor for reducing computational load
!!
!!call q_symmetry(k1x,k1y,k3x,k3y,x4_mod,y4_mod,sym,nnew)
!!
do iloc=1,nlocus
  k2x = x2_mod(iloc)
  k2y = y2_mod(iloc)
  k4x = x4_mod(iloc)
  k4y = y4_mod(iloc)
!
  k2m = sqrt(k2x**2 + k2y**2)
  k4m = sqrt(k4x**2 + k4y**2)
  k2a = atan2(k2y,k2x)*rade
  k4a = atan2(k4y,k4x)*rade
!
  k2m_mod(iloc) = k2m
  k4m_mod(iloc) = k4m
  k2a_mod(iloc) = k2a
  k4a_mod(iloc) = k4a
!
!
!
end do
!
!
9999 continue
!
if(allocated(sold)) deallocate(sold,snew)
!
!!!!!!!!!!!!!!!!!!!!!!!call q_stack('-q_modify')
!
return
end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!--------------------------------!!!!!
!!-------------Q-WEIGHT-----------------!!
!!!!!--------------------------------!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


subroutine q_weight
!------------------------------------------------------------------------------
!
!   +-------+    ALKYON Hydraulic Consultancy & Research
!   |       |    Gerbrant van Vledder
!   |   +---+
!   |   | +---+  Last update: 20 Aug. 2002
!   +---+ |   |  Release: 4.0
!         +---+
!
! do not use m_xnldata
implicit none
!
!  0. Update history
!
!     13/04/1999  Initial version
!     27/10/1999  Weight computed in the case that k2m < k(1)
!     01/11/1999  Use of Q_XK and Q_SK added to compute weights if k > kmax
!     26/11/1999  Bug fixed when checking conversion
!      8/12/1999  Use of SK_MAX introduced to handle very large loci
!     09/08/2002  Modification of weights
!     13/08/2002  storage of log-spacing replace by linear spacing
!     20/08/2002  Bug fixed when geometric scaling is assumed
!
!  1. Purpose:
!
!     Compute interpolation weights of locus
!
!  2. Method
!
!     Compute position of wave number in wave number grid
!     Usable for linear interpolation
!
!  3. Parameter list:
!
!     Name    I/O  Type  Description
!
!  4. Error messages
!
!  5. Called by:
!
!     Q_MAKEGRID
!
!  6. Subroutines used
!
!  7. Remarks
!
!     The tail factors wt_k2 and wt_k4 are valid for the decay of the action density spectrum
!     N(kx,ky). With p (qk_tail) the power p of the tail of the N(k) spectrum, and q the power
!     of the tail of the N(kx,ky) spectrum, we have q=p-1
!
!     Since N(k) = k N(kx,ky) with k the Jacobian
!     it follows that the tail functions are given by
!
!     k^p = k k^q   =>   k^p = k^(q+1) => p=q+1  => q=p-1
!
!  8. Structure
!
!     Initialisations
!     do for all points on locus
!       compute directional index for k2 and k4
!       if geometric scaling then
!         compute wave number index directly
!         convert log-scaling to linear scaling
!       else
!         search position of wave number in k-array
!         if k < kmin then
!           k-index = 1 and factor is 0
!         elsif k < kmax then
!           compute k-index for k2 and k4
!         else
!           compute tail factor
!         end if
!       end if
!     end do
!
!
!  9. Switches
!
!     /T  enable test output
!
! 10. Source code:
!------------------------------------------------------------------------------
!     Local variables
!
integer iloc      ! counter along locus loop
integer jpos      ! index for interpolation and tracking of position in wave numebr array
integer itest     ! local test level
real k2a,k2m      ! angle (radians) and magnitude of wave number k2
real k4a,k4m      ! angle (radians) and magnitude of wave number k2
real dk           ! difference between two wave numbers
real xtest        ! test value for checking computation of weights, by inversion test
real ff,gg        ! variables in transformation of log-spacing to linear spacing
!
! functions used
!
!!real x_kfunc      ! function to invert computation of wieghts
!------------------------------------------------------------------------------
!!!!!!!!!!!!!!!!!!call q_stack('+q_weight')
!
! initialisations
!
itest = iq_test        ! set local test level
itest = itest
!------------------------------------------------------------------------------
do iloc=1,nlocus
  k2m = k2m_mod(iloc)
  k2a = k2a_mod(iloc)
  k4m = k4m_mod(iloc)
  k4a = k4a_mod(iloc)
!
  wt_k2(iloc) = 0.
  wt_k4(iloc) = 0.
!
! compute directional weights
!
  wa_k2(iloc) = (k2a-q_ang1)/q_deltad+1
  wa_k4(iloc) = (k4a-q_ang1)/q_deltad+1
!------------------------------------------------------------------------------
!  compute position of k2 in wave number grid
!  and compute weight function
!-----------------------------------------------------------------------------
  if(iq_disp==1.and. iq_geom==1) then    ! deep water is assumed and loci have geometric scaling
!
    wk_k2(iloc) = 1.+alog(k2m/kqmin)/alog(q_kfac)
    wt_k2(iloc) = 1.
    wk_k4(iloc) = 1.+alog(k4m/kqmin)/alog(q_kfac)
    wt_k4(iloc) = 1.
!
!  Replace log-spacing by linear spacing
!
    ff = wk_k2(iloc)
    gg = floor(ff)
    wk_k2(iloc) = gg+(q_kfac**(ff-gg)-1.)/(q_kfac-1.)
!
!!/T    if(iq_test>=3) write(luq_tst,'(a,4f10.5)') 'Q_WEIGHT: wlog gg wlin2:', &
!!/T &  ff,gg,wk_k2(iloc),abs(wk_k2(iloc)-ff)/abs(ff)*100.
!
    ff = wk_k4(iloc)
    gg = floor(ff)
    wk_k4(iloc) = gg+(q_kfac**(ff-gg)-1.)/(q_kfac-1.)
!
!!/T    if(iq_test>=3) write(luq_tst,'(a,4f10.5)') 'Q_WEIGHT: wlog gg wlin4:', &
!!/T    ff,gg,wk_k4(iloc),abs(wk_k4(iloc)-ff)/abs(ff)*100.
!
!  for finite depth a search is carried out to compute
!  the position of the interacting wave number in the
!  non-geometric k-grid
!
  else
    jpos = 1
    do while (k2m > q_k(jpos))
      jpos = jpos + 1
      if(jpos > nkq) exit
    end do
!
    if(k2m <= q_k(1)) then
      wk_k2(iloc) = k2m/q_k(1)
      wt_k2(iloc) = 0.
    elseif(k2m < q_k(nkq) .and. k2m > q_k(1)) then
      dk          = q_k(jpos)-q_k(jpos-1)
      wk_k2(iloc) = real(jpos-1) + (k2m-q_k(jpos-1))/dk
      wt_k2(iloc) = 1.
    elseif(k2m >= q_k(nkq)) then
      wk_k2(iloc) = min(wk_max,real(nkq) + (k2m-q_k(nkq))/q_sk(nkq))
      wt_k2(iloc) = (k2m/q_k(nkq))**(qk_tail-1.)
!
! minus 1 to account for Jacobian from kx,ky to polar k-grid
!
    end if
!
!  compute position of k4 in wave number grid
!  and compute weight function
!
    jpos = 1
    do while (k4m > q_k(jpos))
      jpos = jpos + 1
      if(jpos > nkq) exit
    end do
!
    if(k4m <= q_k(1)) then
      wk_k4(iloc) = k4m/q_k(1)
      wt_k4(iloc) = 0.
    elseif(k4m < q_k(nkq) .and. k4m > q_k(1)) then
      dk   = q_k(jpos)-q_k(jpos-1)
      wk_k4(iloc) = real(jpos-1) + (k4m-q_k(jpos-1))/dk
      wt_k4(iloc) = 1.
    elseif(k4m >= q_k(nkq)) then
      wk_k4(iloc) = min(wk_max,real(nkq) + (k4m-q_k(nkq))/q_sk(nkq))
      wt_k4(iloc) = (k4m/q_k(nkq))**(qk_tail-1.)
    end if
!
  end if
!
!
end do
!
9999 continue
!
!!!!!!!!!!!!!!!!!!!!!!call q_stack('-q_weight')
!
return
end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!--------------------------------!!!!!
!!-------------Q_LOC_W1W3---------------!!
!!!!!--------------------------------!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


subroutine q_loc_w1w3(k1x,k1y,k3x,k3y,npts,k2x,k2y,k4x,k4y,s)
!-----------------------------------------------------------------
!
!   +-------+    ALKYON Hydraulic Consultancy & Research
!   |       |    Gerbrant van Vledder
!   |   +---+
!   |   | +---+  Last update: 11 June 2003
!   +---+ |   |  Release: 5.0
!         +---+
!
!
implicit none
!
!  0. Update history
!
!     15/04/2002  Initial version
!     20/08/2002  Direction of k1 may be non-zero
!     27/08/2002  Singular solution crosses origin
!     11/06/2003  Length of locus fixed to 3
!
!  1. Purpose:
!
!     Compute locus for the special case w1=w3
!
!  2. Method
!
!     For this case, the k2-locus consists of a straight line
!
!  3. Parameter used:
!
integer, intent(in) :: npts      ! Number of points
real, intent(in)    :: k1x       ! x-component of wave number k1
real, intent(in)    :: k1y       ! y-component of wave number k1
real, intent(in)    :: k3x       ! x-component of wave number k3
real, intent(in)    :: k3y       ! y-component of wave number k3
!
real, intent(out)   :: k2x(npts) ! x-component of wave number k2
real, intent(out)   :: k2y(npts) ! y-component of wave number k2
real, intent(out)   :: k4x(npts) ! x-component of wave number k4
real, intent(out)   :: k4y(npts) ! y-component of wave number k4
real, intent(out)   :: s(npts)   ! distance along locus
!
!  4. Error messages
!
!  5. Caled by:
!
!     Q_CMPLOCUS
!
!  6. Subroutines used
!
!  7. Remarks
!
!     Routine based on modified version of routine SHLOCX of Resio and Tracy
!     On 15/4/2002 a bug fixed in computation of THR when angle of k3 is larger than 90°
!
!     In addition, the assumption that k1y=0 and thus dir1=0 is removed
!     In bug fix of 20/8/2002 this restriction is removed.
!
!  8. Structure
!
!     Compute angle of symmetry axis
!     Compute distance between 2 lines of solution
!     compute wave numbers along locus
!     rotate angles
!
!  9. Switches
!
! 10. Source code
!------------------------------------------------------------------------------
!     Local variables
!
integer ipt   ! counter of points along locus
!
real dirs     ! angle of symmetry axis
real dir1     ! direction of wave number k1
real dir3     ! direction of wave number k3
real dk0      ! step size along locus
real xk0      ! x-component
real yk0      ! y-component
real w2       ! radian frequency
real xx2,yy2  ! values along k2-locus
real xx4,yy4  ! values along k4-locus
real k1m      ! magnitude of wave number k1
!------------------------------------------------------------------------------
!
!    dirs is the angle of rotation from the x-axis to the "bisecting" angle
!
dir1 = atan2(k1y,k1x)
dir3 = atan2(k3y,k3x)
dirs = 0.5*(180-abs(180-abs(dir3-dir1)))
k1m  = sqrt(k1x**2 + k1y**2)
!
!     k1x is the total length of the wavenumber vector
!     xk0 is the length of this vector in the rotated coordinate system
!
xk0 = k1m * cos(dirs)
yk0 = k1m * sin(dirs)
!
! Specify step size for solution of singular case
!
!! dk0 = 0.11    ! Removed on 11/6/2003, this value is used in original WRT code
!! dk0 = kqmax/real(npts-1.)             this value depends on actual grid
dk0 = 3./real(npts-1.)                 !  this is test value
!
!  modify rotation angle
!
dirs = dirs + dir1
!
!  generate sequence of parallel lines
!  rotate lines over modified angle DIRS
!
do ipt=1,npts
!  w2       = real(ipt-1.)*dk0         ! removed on Aug. 27 2002
!
  w2       = 2.*real(ipt-npts/2)*dk0   ! create line on both sides of origin
  xx2      = w2*xk0
  yy2      = yk0
  k2x(ipt) = xx2*cos(dirs) - yy2*sin(dirs)
  k2y(ipt) = yy2*cos(dirs) + xx2*sin(dirs)
  xx4      = xx2
  yy4      = -yy2
  k4x(ipt) = xx4*cos(dirs) - yy4*sin(dirs)
  k4y(ipt) = yy4*cos(dirs) + xx4*sin(dirs)
  s(ipt)   = real(ipt-1)*dk0*xk0
end do
!
return
end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!--------------------------------!!!!!
!!-------------Q-LOCPOS-----------------!!
!!!!!--------------------------------!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


subroutine q_locpos(ka,kb,km,kw,loclen)
!------------------------------------------------------------------------------
!
!   +-------+    ALKYON Hydraulic Consultancy & Research
!   |       |    Gerbrant van Vledder
!   |   +---+
!   |   | +---+  Last update: 14 Oct. 2002
!   +---+ |   |  Release: 5.0
!         +---+
!
! do not use m_xnldata
use m_constants
use serv_xnl4v5, only: z_root2
!
implicit none
!
!  0. Update history
!
!     Version Date        Description
!
!     03/12/1999  Initial version
!     09/08/2002  Upgrade to release 4.0
!     29/08/2002  Error handling z_root2 relaxed and some write statements modified
!     07/10/2002  Initialisation of QSQ replaced
!
!  1. Purpose:
!
!     Compute characteristics of locus used to optimize its acutal computation
!
!  2. Method
!
!  3. Parameter list:
!
!Type    I/O          Name     Description
!-----------------------------------------------------------------
real, intent (out) :: ka     ! minimum k along symmetry axis
real, intent (out) :: kb     ! maximum k along symmetry axis
real, intent (out) :: km     ! wave number at midpoint
real, intent (out) :: kw     ! half width of locus at midpoint
real, intent (out) :: loclen ! estimated length of locus
!
!  4. Error messages
!
!  5. Called by:
!
!     Q_CMPLOCUS
!
!  6. Subroutines used
!
!     z_zero2   Root finding method
!     x_locus1  Function of locus geometry, along symmetry axis
!     x_locus2  Function of locus geometry, perpendicular to symmetry axis
!     x_flocus  Locus function
!
!  7. Remarks
!
!  8. Structure
!
!  9. Switches
!
!     /S  enable subroutine tracing
!     /T  enable test output
!
! 10. Source code
!------------------------------------------------------------------------------
!     Local variables
!
real kp           ! wave number at peak
real kpx,kpy      ! wave number at peak maximum
real zp           ! value of locus function at maximum
real za,zb        ! (test) value of locus function at kmin & kmax
real zz1,zz2      ! intermediate function values in interation process
real kk1,kk2      ! start values for finding root of locus equation
real kk1x,kk1y    ! wave number components at one side of root
real kk2x,kk2y    ! wave number components at other side of root
real beta1,beta2  ! parameters specifying cross component
real betaw        ! parameter specifying iterated cross component
real kwx,kwy      ! wave number at side of locus
real zw           ! function value at (kwx,kwy)
real a1,a2,b1,b2  ! constants in polynomial approximation of elliptic function
real aa,bb,mm,mm1 ! semi-major exis of ellips and derived parameters
!
real eps         ! local machine accuracy for determination of roots
real bacc        ! accuracy for determination of beta
real kacc        ! accuracy for determination of wave number roots
real qs          ! (w1-w3)/sqrt(g)
real qsq         ! gs^2
!
! Function declaration
!
!real z_root2     ! root finding using Ridders method
!
integer ierr     ! local error indicator, used in function Z-ZERO1
integer itest    ! local test level for test output
integer lutest   ! unit for test output in service routines
integer iter     ! local iteration number
integer maxiter  ! maximum number of iteration for determining starting points
!
!  function declarations
!!real, external :: x_locus2    ! locus function perpendicular to symmetry axis
!!real x_flocus                 ! 2-d locus function
!---------------------------------------------------------------------------------
!  assign test options
!
itest  = iq_test              ! assign test level
lutest = 0                    ! assign default, no test output in service routines
!
itest  = 0                    ! reset local test level
if(itest > 0) lutest=luq_tst   ! assign unit for test output
!
!!!!!!!!!!!!!!!call q_stack('+q_locpos')
!
!  set initial values
!
eps     = epsilon(1.)         ! determine machine accurcy
maxiter = 20                  ! maximum number of iterations
!
! compute location of maximum, located at k_2 = P
!
kpx  = -px
kpy  = -py
kp   = sqrt(kpx**2 + kpy**2)
zp   = x_locus1(kp)
!
! find location of points A and B on locus
! for deep water, explicit relations are available
!
if(iq_disp==1) then
  qs = q/sqrtg
  qsq  = qs*qs
  if(qs < 0) then
    ka = 0.5*(-qs+sqrt(2.0*pmag-qsq))
    ka = ka**2
    kb = (pmag+qsq)/(2.*qs)
    kb = kb**2
    za = x_locus1(ka)
    zb = x_locus1(kb)
  else
    ka = 0.5*(-qs+sqrt(2.0*pmag-qsq))
    ka = -ka**2
    kb = (pmag-qsq)/(2.*qs)
    kb = kb**2
    za = x_locus1(ka)
    zb = x_locus1(kb)
  end if
!
!
!  find location of points A and B on locus
!  for water of finite depth, an iteration process is applied to
!  determine the zero-crossings of the locus function
!
else
!
  if(q<0) then
!
!   set two start points to locate position of wave number ka
!
    kk1 = 0.
    kk2 = kp
!
!   search root by Ridder's method
!
    kacc = 10.*max(kk1,kk2)*eps
    ka = z_root2(x_locus1,kk1,kk2,kacc,lutest,ierr)
!
!
!
!   determine start points to locate position of wave number kb
!
    kk1 = kp
    kk2 = kp
    zz1 = zp
    zz2 = zp
    iter = 0
!
!  ensure that two points are found on either side of zero-crossing
!
    do while (zz1*zz2 >= 0 .and. iter < maxiter)
      iter = iter + 1
      kk2 = kk2*2
      zz2 = x_locus1(kk2)
    end do
!
    if(iter>=maxiter) then
      !!!!!!!!!!!!!call q_error('e','Start kb','Too many iterations needed')
      print *, "Error 46129"
      goto 9999
    end if
!
!   search root by Ridders method
!
    kacc = 10.*max(kk1,kk2)*eps
    kb = z_root2(x_locus1,kk1,kk2,kacc,lutest,ierr)
!
!==================================================================
!   find positions for ka and kb for the case q > 0
!
  else
!
!   set two start points to locate position of wave number ka
!
    kk1  = 0.
    kk2  = -kp
    zz1  = x_locus1(kk1)
    zz2  = x_locus1(kk2)
    iter = 0
!
!  ensure that two points are found on either side of zero-crossing
!
    do while (zz1*zz2 >= 0 .and. iter < maxiter)
      iter = iter + 1
      kk2 = kk2*2
      zz2 = x_locus1(kk2)
    end do
!
    if(iter>=maxiter) then
      !!!!!!!!!!!!call q_error('e','Start ka','Too many iterations needed')
      print *, "Error 923571"
      goto 9999
    end if
!
!   search root by Ridder's method
!
    kacc = 10.*max(abs(kk1),abs(kk2))*eps
    ka = z_root2(x_locus1,kk1,kk2,kacc,lutest,ierr)
!
!   determine start points to locate position of wave number kb
!
    kk1  = 0
    kk2  = kp
    zz1  = x_locus1(kk1)
    zz2  = x_locus1(kk2)
    iter = 0
!
!  ensure that two points are found on either side of zero-crossing
!
    do while (zz1*zz2 >= 0 .and. iter < maxiter)
      iter = iter + 1
      kk2 = kk2*2
      zz2 = x_locus1(kk2)
    end do
!
    if(iter>=maxiter) then
      !!!!!!!!!!!!call q_error('e','Start kb','Too many iterations needed')
      print *, "Error 436781"
      goto 9999
    end if
!
!   search root by Ridders method
!
    kacc = 10.*max(kk1,kk2)*eps
    kb = z_root2(x_locus1,kk1,kk2,kacc,luq_tst,ierr)
!
!   find positions for ka and kb for the case q > 0
!
  end if
!
  za = x_locus1(ka)
  zb = x_locus1(kb)
!
end if
!
! compute position of mid point
!
kmid = 0.5*(ka+kb)
km   = kmid
!
if(q < 0) then
  kmidx = kmid*cos(pang+pi)
  kmidy = kmid*sin(pang+pi)
else
  kmidx = kmid*cos(pang)
  kmidy = kmid*sin(pang)
end if
!
!
! compute width of locus near mid point of locus
!
! set starting values for determination of crossing point
!
beta1 = 0.
kk1x  = kmidx
kk1y  = kmidy
beta2 = 0.5
kk2x  = kmidx - beta2*py
kk2y  = kmidy + beta2*px
zz1   = x_flocus(kk1x,kk1y)
zz2   = x_flocus(kk2x,kk2y)
!
!
iter = 0
do while (zz1*zz2 > 0 .and. iter < maxiter)
  iter = iter + 1
  kk2x = kmidx - beta2*py
  kk2y = kmidy + beta2*px
  zz1  = x_flocus(kk1x,kk1y)
  zz2  = x_flocus(kk2x,kk2y)
  beta2 = beta2*2
end do
!
! call Ridders method to locate position of zero-crossing
!
!
bacc = 10.*max(beta1,beta2)*eps
betaw = z_root2(x_locus2,beta1,beta2,bacc,lutest,ierr)
!
!
kwx = kmidx - betaw*py
kwy = kmidy + betaw*px
zw  = x_flocus(kwx,kwy)
kw  = betaw*pmag
!
!
! estimate circumference of locus, assuming it to be an ellips
! estimate axis, this seems to be a rather good estimate
!
aa = 0.5*abs(ka-kb)
bb = kw
!
if (aa > bb) then
  mm = 1-(bb/aa)**2
else
  mm = 1-(aa/bb)**2
end if
!
mm1 = 1.-mm
a1 = 0.4630151;  a2 = 0.1077812;
b1 = 0.2452727;  b2 = 0.0412496;
!
if (mm1==0) then
  loclen = 4.*max(aa,bb)
else
  loclen = 4.*max(aa,bb)*((1. + a1*mm1 + a2*mm1**2) + (b1*mm1 + b2*mm1**2)*log(1/mm1))
end if
!
!
9999 continue
!
!!!!!!!!!!!!!!!!!call q_stack('-q_locpos')
!
return
end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!--------------------------------!!!!!
!!-------------Q-POLAR2-----------------!!
!!!!!--------------------------------!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine q_polar2(kmin,kmax,kx_beg,ky_beg,kx_end,ky_end,loclen,ierr)
!------------------------------------------------------------------------------
!
!   +-------+    ALKYON Hydraulic Consultancy & Research
!   |       |    Gerbrant van Vledder
!   |   +---+
!   |   | +---+  Last update: 8 Aug. 2003
!   +---+ |   |  Release: 5.0
!         +---+
!
! do not use m_xnldata
use m_constants
use serv_xnl4v5, only: z_wnumb
!
implicit none
!
!  0. Update history
!
!     Date        Description
!
!     03/12/1999  Initial version
!     09/08/2002  Geometric spacing of k added
!                 Upgrade to release 4.0
!     13/08/2002  reorganisation of loops generating points on locus
!     08/08/2003  Check included for maximum number of IPOL by using MPOL
!                    MPOL=MLOCUS/2+1-1  (-1 added regarding IPOL=IPOL+1 in Q_MODIFY)
!                 Check included on ARG=0 for IQ_LOCUS=2 and parameter dke added
!
!  1. Purpose:
!
!     Compute position of locus for given k1-k3 vector
!
!  2. Method
!
!     Explicit polar method, see Van Vledder 2000, Monterey paper
!     Optionally using a fixed k-step, geometric k-step or adaptive stepping
!
!  3. Parameters used:
!
!Type    I/O        Name             Description
!------------------------------------------------------------------------------
real, intent(in) :: kmin           ! minimum wave number on locus
real, intent(in) :: kmax           ! maximum wave number on locus
real, intent(in) :: kx_beg         ! x-coordinate of begin point
real, intent(in) :: ky_beg         ! y-coordinate of begin point
real, intent(in) :: kx_end         ! x-coordinate of end point
real, intent(in) :: ky_end         ! y-coordinate of end point
real, intent(in) :: loclen         ! estimated length of locus
integer, intent (out)  :: ierr     ! error condition
!
!     Parameters with module
!
!     nlocus0   Preferred number of points on locus
!     q         w1-w3, difference of radian frequencies
!     pmag      |k1-k3| (vector form)
!     pdir      direction of difference vector k1-k3
!
!  4. Error messages
!
!  5. Called by:
!
!     Q_CPMLOCUS
!
!  6. Subroutines used:
!
!     X_COSK
!
!  7. Remarks
!
!     The type of locus computation is controlled by the parameter IQ_LOCUS
!     Set in Q_SETCFG
!
!  8. Structure
!
!  9. Switches
!
!     /S  enable subroutine tracing
!     /T  enable test output
!
! 10. Source code
!------------------------------------------------------------------------------
!     Local variabels
!
integer ipol      ! counter
integer jpol      ! counter
integer iend      ! indicates end of locus computation
integer ipass     ! counter for passes
integer npol      ! number of points on locus
integer npass     ! number of passes in iteration process
integer mpol      ! maximum number of points on locus, related to MLOCUS
!
real kold         ! temporary wave number
real knew         ! temporary wave number
real cosold       ! 'old' cosine of angle
real cosnew       ! 'new' cosine of angle
real dkpol        ! step in wave number
real dkold        ! 'old' step in wave number
real ang1         ! 'old' angle
real ang2         ! 'new' angle
real kk1          ! 'old' wave number
real kk2          ! 'new' wave number
real kratio       ! ratio between succesive k-values when IQ_LOCUS=3
real arg          ! argument
real dk           ! step in wave number
real dke          ! estimate of new dk
real dsnew        ! new step size along locus
real dsz          ! estimated step size along locus
!
integer itest    ! local test level
integer lutest   ! unit number for test output in service routine
!
!  function declarations
!!!real    z_wnumb  ! compute wave number, via module SERV_XNL4V4
!!real    x_disper ! dispersion relation
!
!------------------------------------------------------------------------------
! initialisations
!------------------------------------------------------------------------------
!!!!!!!!!!!!!!!!!call q_stack('+q_polar2')
!
ierr = 0                      ! set error code to zero
npol = (nlocus0+1)/2+1        ! first estimate of number k-values along symmetry axis
mpol = mlocus/2               ! set maximum number of points along locus axis
!
!-------------------------------------------------------------------------------
!
select case(iq_locus)
!------------------------------------------------------------------------------
! CASE = 1: Linear spacing of wave numbers along symmetry axis
!------------------------------------------------------------------------------
  case(1)
!
  dk = (kmax-kmin)/real(npol-1)
  do ipol=1,npol
    k_pol(ipol) = kmin + (ipol-1)*dk
    c_pol(ipol) = x_cosk(k_pol(ipol))
  end do
!------------------------------------------------------------------------------
!  Case = 2: Variable k-stepping along symmetry axis,
!            such that step along locus is more or less constant
!------------------------------------------------------------------------------
  case(2)
!
! set first point on locus
!
  ipol        = 1
  k_pol(ipol) = kmin
  c_pol(ipol) = -1.
  kold        = kmin
  cosold      = -1.
!
! compute initial step size of polar wave number
!
  dk0   = (kmax - kmin)/real(npol)      ! estimate of step size of equidistant radii
  dsz   = loclen/real(nlocus0)          ! estimate of step size along locus
  npass = 3                             ! set number of passes in iteration
  dk0   = dk0/2                         ! reduce initial step
  dk    = dk0
  iend  = 0
!
!
  do while (k_pol(ipol) < kmax .and. iend==0 .and. ipol < mpol)
    do ipass=1,npass
      knew  = min(kmax,k_pol(ipol)+dk)
      dkold = knew - k_pol(ipol)
      cosnew = x_cosk(knew)
      ang1  = pang + acos(cosold)
      ang2  = pang + acos(cosnew)
      kk1   = kold
      kk2   = knew
      arg   = kk1**2 + kk2**2 -2.*kk1*kk2*cos(ang1-ang2)
      dsnew = sqrt(abs(arg))
      if(dsnew>0) dke   = dk*dsz/dsnew
      dk    = dke
    end do
!----------------------------------------------------------------------------------------------
!  assign new estimate and check value of IPOL
!----------------------------------------------------------------------------------------------
    ipol        = ipol + 1
    k_pol(ipol) = k_pol(ipol-1) + dkold
    c_pol(ipol) = cosnew
    kold        = knew
    cosold      = cosnew
    if (abs(dkold) < 0.0005*(kmax-kmin)) iend=1
  end do
!
! fill last bin with coordinates of end point
!
  if(k_pol(ipol) < kmax .and. ipol <  mpol) then
    ipol = ipol + 1
    c_pol(ipol) = -1.
    k_pol(ipol) = kmax
  end if
!
!  update the number of k-points on symmetry axis
!
  npol = ipol
!
!-------------------------------------------------------------------------------
!  Case 3: Geometric spacing of wave numbers along symmetry axis
!-------------------------------------------------------------------------------
  case(3)
  kratio = (kmax/kmin)**(1./(npol-1.))
  do ipol=1,npol
    k_pol(ipol) = kmin*kratio**(ipol-1.)
    c_pol(ipol) = x_cosk(k_pol(ipol))
  end do
!
end select
!
!------------------------------------------------------------------------------
!
!  compute actual number of points on locus
!  this will always be an even number
!  mirror image the second half of the locus
!
nlocus1 = 2*npol-2
!
a_pol(1) = pang + acos(c_pol(1))
c_pol(1) = cos(a_pol(1))
!
do ipol=2,npol
  jpol = 2*npol-ipol
  a_pol(ipol) = pang + acos(c_pol(ipol))
  a_pol(jpol) = pang - acos(c_pol(ipol))
  c_pol(jpol) = cos(a_pol(jpol))
  k_pol(jpol) = k_pol(ipol)
end do
!
! compute x- and y-position along locus
!
do ipol=1,nlocus1
  x2_loc(ipol) = k_pol(ipol)*cos(a_pol(ipol))
  y2_loc(ipol) = k_pol(ipol)*sin(a_pol(ipol))
end do
!
!
9999 continue
!
!!!!!!!!!!!!!!!!!call q_stack('-q_polar2')
!
return
end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!--------------------------------!!!!!
!!-------------Q-SYMMETRY---------------!!
!!!!!--------------------------------!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


subroutine q_symmetry(k1x,k1y,k3x,k3y,k4x,k4y,symfac,nloc)
!------------------------------------------------------------------------------
!
!   +-------+    ALKYON Hydraulic Consultancy & Research
!   |       |    Gerbrant van Vledder
!   |   +---+
!   |   | +---+  Last update: 16 June 2003
!   +---+ |   |  Release: 5.0
!         +---+
!
implicit none
!--------------------------------------------------------------------------------
!  0. Update history
!
!     10/06/2003  Initial version
!     16/06/2003  Switch iq_sym added
!
!  1. Purpose:
!
!     Compute symmetry factor to reduce integration
!
!  2. Method
!
!     Compute distance between k1 and k3, and between k4 and k1
!
!  3. Parameter list:
!
! Type   i/o             Name           Description
!----------------------------------------------------------------------------------
integer, intent(in)   :: nloc         ! number of points in array with wave number
real, intent(in)      :: k1x          ! x-component  of wave number k1
real, intent(in)      :: k1y          ! y-component  of wave number k1
real, intent(in)      :: k3x          ! x-component  of wave number k3
real, intent(in)      :: k3y          ! y-component  of wave number k3
real, intent(in)      :: k4x(nloc)    ! x-components of wave number k4
real, intent(in)      :: k4y(nloc)    ! y-components of wave number k4
real, intent(out)     :: symfac(nloc) ! symmetry factor
!----------------------------------------------------------------------------------
!  4. Error messages
!
!  5. Called by:
!
!     Q_MODIFY
!
!  6. Subroutines used
!
!     Q_STACK
!
!  7. Remarks
!
!  8. structure
!
!  9. Switches
!
! 10. Source code
!------------------------------------------------------------------------------
integer iloc      ! counter
real dk13         ! distance between k1 and k3
real dk14         ! distance between k1 and k4
!------------------------------------------------------------------------------
!
!!!!!!!!!!!!!!!!!!!call q_stack('+q_symmetry')
!
!
! evaluate criterion |k3-k1| < |k4-k1|
! if true then symfac=1
!
symfac = 1.
if(iq_sym==1) then
  dk13 = (k1x-k3x)**2 + (k1y-k3y)**2
  do iloc=1,nloc
    dk14 = (k1x-k4x(iloc))**2 + (k1y-k4y(iloc))**2
    if (dk13 >= dk14) symfac(iloc) = 0.
  end do
end if
!
!!!!!!!!!!!!!!!!!!call q_stack('-q_symmetry')
!
return
end subroutine


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!--------------------------------!!!!!
!!-------------Q-SEARCHGRID-------------!!
!!!!!--------------------------------!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine q_searchgrid(depth,igrid)
!------------------------------------------------------------------------------
!
!   +-------+    ALKYON Hydraulic Consultancy & Research
!   |       |    Gerbrant van Vledder
!   |   +---+
!   |   | +---+  Last update: 9 September 2003
!   +---+ |   |  Release: 5.03
!         +---+
!
! do not use m_xnldata
implicit none
!------------------------------------------------------------------------------
!  0. Update history
!
!     Version  Date    Modification
!
!     20/08/2002  Initial version
!     29/08/2002  Write statements made conditionsl
!      5/09/2003  Search algorithm improved
!     09/09/2003  factor ID_FACMAX introduced and extra test output created
!                 Input water depth saved for output
!
!  1. Purpose:
!
!     Search nearest valid grid, read grid file and scale factor
!
!  2. Method
!
!     Using the actual water depth
!     all possible interaction grids are checked
!     in upward and downward direction
!
!  3. Parameters used
!
real, intent(in)     :: depth  !  depth for which grid file must be found
integer, intent(out) :: igrid  !  status of grid checking
!                                 ==0: a proper grid exists
!                                 ==1: grid file does not exist
!                                 ==2: grid file exists, but it is incorrect
!                                 ==3: read error in accessing grid information
!
!  4. Error messages
!
!  5. Called by:
!
!     Q_XNL4V4
!
!  6. Subroutines used
!
!     Q_CTRGRID
!     Q_STACK
!
!  7. Remarks
!
!
!  8. Structure
!
!
!  9. Switches
!
! 10. Source code
!---------------------------------------------------------------------------
!     Local variables
!
integer id        ! counter
integer idepth    ! integer depth
integer id_upper  ! upper limit of search
integer id_lower  ! lower limit of depth search
!
real d_lower      ! lower valid depth
real d_upper      ! upper valid depth
real r_lower      ! ratio with lower valid depth
real r_upper      ! ratio with upper valid depth
real s_depth      ! target depth in m, saved in this variable
real dfac1,dfac2  ! depth scale factors
real eps          ! accuracy
!------------------------------------------------------------------------------
!
!!!!!!!!!!!!!!!!call q_stack('+q_searchgrid')
!
eps = 0.0001
!
!------------------------------------------------------------------------------
!  check if a depth exists for current grid
!------------------------------------------------------------------------------
!
!
q_depth = depth + eps

call q_ctrgrid(1,igrid)
!
!
if(igrid==0) then
  if(iq_screen>=1) write(iscreen,'(a)') 'Q_SEARCHGRID: grid accepted, read whole database'
!
  call q_ctrgrid(2,igrid)
  goto 9999
end if
!
! save depth for which nearest grid file is to be found
!
s_depth  = depth
idepth   = int(s_depth*10+eps)
id_lower = int(q_mindepth*10+eps)
id_upper = int(q_maxdepth*10+eps)
!
id_upper = min(id_facmax*idepth,id_upper)
!
!  set 'not found' condition
!
d_lower = -1.
d_upper = -1.
!
!------------------------------------------------------------------------------
! search downwards until a valid grid is found
!------------------------------------------------------------------------------
!
do id = idepth-1,id_lower,-1
  q_depth = real(id)/10.+eps

  call q_ctrgrid(1,igrid)


  if(igrid==0) then
    d_lower = q_depth
    exit
  end if
end do
!
!------------------------------------------------------------------------------
!  seach upwards until a valid grid is found
!------------------------------------------------------------------------------
!
do id = idepth+1,id_upper
  q_depth = real(id)/10.+eps


  call q_ctrgrid(1,igrid)


  if(igrid==0) then
    d_upper = q_depth
    exit
  end if
end do
if(iq_prt>=1) write(luq_prt,*)
!------------------------------------------------------------------------------
!
!  determine nearest grid
!------------------------------------------------------------------------------
!
if(d_lower > 0) then
  r_lower = s_depth/d_lower
else
  r_lower = -1.
end if
!
if(d_upper > 0) then
  r_upper = d_upper/s_depth
else
  r_upper = -1.
end if
!
if(iq_prt>=1) then
  write(luq_prt,'(a,3f8.2)') 'Q_SEARCHGRID: d_lower d_target d_upper      :',d_lower,s_depth,d_upper
  write(luq_prt,'(a,2f8.2)') 'Q_SEARCHGRID: r_lower r_upper               :',r_lower,r_upper
end if
!------------------------------------------------------------------------------
!  select nearest valid grid
!------------------------------------------------------------------------------
if(r_lower>0 .and. r_upper>0) then
  if(r_lower < r_upper) then
    q_depth = d_lower
  else
    q_depth = d_upper
  end if
!
elseif(r_lower > 0 .and. r_upper <0 ) then
  q_depth = d_lower
elseif(r_lower < 0 .and. r_upper > 0) then
  q_depth = d_upper
else
  !!!!!!!!!!!!!!!!call q_error('e','SEARCHGRID','No valid nearest grid could be found')
  print *, "Error 127643"
  goto 9999
end if
!
!-----------------------------------------------------------------------------------------------
! compute depth scaling factors
!------------------------------------------------------------------------------
!
call q_dscale(a,q_sig,q_a,nkq,naq,s_depth,q_grav,dfac1)
call q_dscale(a,q_sig,q_a,nkq,naq,q_depth,q_grav,dfac2)
!
q_scale = dfac1/dfac2
!
if(iq_prt>=1) then
  write(luq_prt,'(a,2f8.4)') 'Q_SEARCHGRID: target and nearest scale factors:',dfac1,dfac2
  write(luq_prt,'(a,f8.4)')  'Q_SEARCHGRID: compound scale factor           :',q_scale
end if
!
!  Read BQF for nearest valid water depth
!
call q_ctrgrid(2,igrid)
if(iq_prt>=2) then
  write(luq_prt,'(a,f12.2)') 'Q_SEARCHGRID: Q_CTRGRID called with depth:',q_depth
  write(luq_prt,'(a,i4)') 'Q_SEARCHGRID: igrid of nearest grid operation:',igrid
end if
!
9999 continue
!
!  restore water depth
!
q_depth = s_depth
!
!!!!!!!!!!!!!call q_stack('-q_searchgrid')
!
return
end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!--------------------------------!!!!!
!!-------------Q-DSCALE-----------------!!
!!!!!--------------------------------!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine q_dscale(n,sigma,angle,nsig,nang,depth,grav,q_dfac)
!------------------------------------------------------------------------------
!
!   +-------+    ALKYON Hydraulic Consultancy & Research
!   |       |    Gerbrant van Vledder
!   |   +---+
!   |   | +---+  Last update: 23 Aug. 2002
!   +---+ |   |  Release: 5.0
!         +---+
!
use serv_xnl4v5
implicit none
!
!  0. Update history
!
!      Date        Modification
!
!      25/02/1999  Initial version
!       2/12/1999  Result modified if total energy <= 0
!                  Cosmetic changes
!      13/08/2002  Upgrade to release 4.0
!      23/09/2002  Mean wave number multiplied by 0.75
!
!  1. Purpose:
!
!     Compute scaling factor for nonlinear transfer in finite depth
!
!  2. Method
!
!     Compute mean wave number km
!
!     Compute scale factor based on parameterized function of (km*d)
!     according to Herterich and Hasselmann
!     and parameterisation from WAM model
!
!
!  3. Interface parameter list:
!
! Type          I/O     Name           Description
!-------------------------------------------------------------------------
integer, intent (in) :: nsig          ! Number of sigma-values
integer, intent (in) :: nang          ! Number of directions
real, intent(in)     :: n(nsig,nang)  ! N(nsig,nang) Action density
real, intent(in)     :: sigma(nsig)   ! sigma values
real, intent(in)     :: angle(nang)   ! directions in (radians)
real, intent(in)     :: depth         ! Depth (m)
real, intent(in)     :: grav          ! Gravitational acceleration
real, intent(out)    :: q_dfac        ! scale factor
!
!  4. Error messages
!
!  5. Called by:
!
!     XNL_MAIN
!
!  6. Subroutines used
!
!     x_wnumb
!     z_steps
!     q_stack
!
!  7. Remarks
!
!  8. Structure
!
!  9. Switches
!
! 10. Source code
!------------------------------------------------------------------------------
!     local variables
!
real w          ! radian frequency
real kk         ! local wave number
real sqkk       ! square root of local wave number
real dnn        ! summation quantity
real kms        ! mean wave number
real kd         ! depth*mean wave number product
real sum0       ! summation variable for total energy
real sumk       ! summation variable for wave number
real delta      ! directional step, in radians
!
integer isig    ! counter over sigma loop
integer iang    ! counter over direction loop
!
! functions
!!!real z_wnumb    ! function to compute wave number
!
! temporary data
!
real dsigma(nsig) ! step size of sigma array, used for integration
!------------------------------------------------------------------------------
!
!!!!!!!!!!!!!!!call q_stack('+q_dscale')
!
call z_steps(sigma,dsigma,nsig)   !  compute step size of sigma's
delta = angle(2)-angle(1)         !  compute directional step (radians)
!
sum0 = 0.
sumk = 0.
!
!  compute sums for total energy andwave number
!
do isig = 1,nsig
  w    = sigma(isig)
  kk   = z_wnumb(w,depth,grav)   ! compute wave number for given sigma,depth
  sqkk = sqrt(kk)
  do iang=1,nang
    dnn  = n(isig,iang)*dsigma(isig)*delta
    sum0 = sum0 + dnn
    sumk = sumk + 1./sqkk*dnn
  end do
end do
!
!  compute mean wave number and scale factor based
!  on the WAM approximation
!
if(sum0 > 0) then
  kms = (sum0/sumk)**2
  kd = max(0.5,0.75*kms*depth)
  q_dfac = 1+5.5/kd*(1.-5./6.*kd)*exp(-5./4.*kd)
!  pause
else
  kms = 0.
  kd  = 0.
  q_dfac = 1.
end if
!
!!!!!!!!!!!!!!!call q_stack('-q_dscale')
!
return
end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!----Z--SUBROUTINES---------------------!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!--------------------------------!!!!!
!!-------------Z-CMPCG------------------!!
!!!!!--------------------------------!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine z_cmpcg(sigma,depth,grav,cg)
!-----------------------------------------------------------------------------!
!
!   +-------+    ALKYON Hydraulic Consultancy & Research
!   |       |    Gerbrant van Vledder
!   |   +---+
!   |   | +---+
!   +---+ |   |
!         +---+
!
implicit none
!
!  0. Update history
!
!     12/01/2001  Initial version
!     11/04/2001  Check included for the cases tat sigma < 0 or depth <0
!                 Result is cg = -10
!
!  1. Purpose:
!
!     Compute group velocity for a given radian frequency and depth
!
!  2. Method
!
!     Linear wave theory
!
!  3. Parameter list:
!
!Type   I/O           Name    Description
!------------------------------------------------------------------------------
real, intent(in)  :: sigma  !  radian frequency (rad)
real, intent(in)  :: depth  !  water depth (m)
real, intent(in)  :: grav   !  gravitational acceleration (m/s^2)
real, intent(out) :: cg     !  group velocity (m/s)
!
real k                      !  wave number
!/A
!!real z_wnumb                !  compute wave number
!/Z
!-----------------------------------------------------------------------------
k = z_wnumb(sigma,depth,grav)
!
if(depth <= 0. .or. sigma <= 0.) then
  cg = -10.
else
  if(depth*k > 30.) then
    cg = grav/(2.*sigma)
  else
    cg = sigma/k*(0.5+depth*k/sinh(2.*depth*k))
  end if
end if
!
return
end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!--------------------------------!!!!!
!!-------------Z-STEPS------------------!!
!!!!!--------------------------------!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine z_steps(x,dx,nx)
!-----------------------------------------------------------------------------!
!
!
!   +-------+    ALKYON Hydraulic Consultancy & Research
!   |       |    Gerbrant van Vledder
!   |   +---+
!   |   | +---+  Creation date:  September 28, 1998
!   +---+ |   |  Last Update:    march 19, 2003
!         +---+
!
!  0. Update history
!
!     19/03/2003   Input argument defined using intent option
!                  check included nx > 0
!
!  1. Purpose
!
!     Compute bandwidth of spectral discretization
!
implicit none
!
integer, intent(in) :: nx     ! Number of elements in array
real, intent(in)    :: x(nx)  ! Input data array with elements
real, intent(out)   :: dx(nx) ! Output array with step sizes
!
integer ix                    ! counter
!------------------------------------------------------------------------------
if (nx<1) then
  return
!
elseif (nx==1) then
  dx = 0
else
  do ix=2,nx-1
    dx(ix) = 0.5 * (x(ix+1) - x(ix-1))
  end do
!
  if (nx >= 4) then
    dx(1)  = dx(2)*dx(2)/dx(3)
    dx(nx) = dx(nx-1)*dx(nx-1)/dx(nx-2)
  else
    dx(1)  = dx(2)
    dx(nx) = dx(nx-1)
  end if
end if
!
return
end subroutine


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!--------------------------------!!!!!
!!-------------Z_POLYAREA---------------!!
!!!!!--------------------------------!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


subroutine z_polyarea(xpol,ypol,npol,area)
!-----------------------------------------------------------------------------!
!
!   +-------+    ALKYON Hydraulic Consultancy & Research
!   |       |    P.O. Box 248
!   |   +---+    8300 AE Emmeloord
!   |   | +---+  Tel: +31 527 620909
!   +---+ |   |  Fax: +31 527 610020
!         +---+  http://www.alkyon.nl
!
!         Gerbrant van Vledder
!
!  0. Update history
!
!     0.01  12/06/2003  Initial version
!
!  1. Purpose
!
!     Computes area of a closed polygon
!
!  2. Method
!
!     The area of the polygon
!
!  3. Parameter list
!
!     Name    I/O  Type  Description
!
integer, intent(in)  ::  npol       ! Number of points of polygon
real, intent(in)     ::  xpol(npol) ! x-coodinates of polygon
real, intent(in)     ::  ypol(npol) ! y-coordinates of polygon
real, intent(out)    ::  area       ! area of polygon
!
!  4. Subroutines used
!
!  5. Error messages
!
!  6. Remarks
!
integer ipol,ipol1         ! counters
real xmin,xmax,ymin,ymax   ! minima and maxima of polygon
real xmean,ymean           ! mean values
real xa,ya,xb,yb           ! temporary variables
real sumx,sumy             ! sums
real darea                 ! piece of area
!-------------------------------------------------------------------------------
if(npol<=1) then
  crf  = 0.
  xz   = 0.
  yz   = 0.
  area = 0.
  return
end if
!
! compute minimum and maximum coordinates
!
xmin = minval(xpol)
xmax = maxval(xpol)
ymin = minval(ypol)
ymax = maxval(ypol)
!
!  compute mean of range of x- and y-coordinates
!
xmean = 0.5*(xmin + xmax)
ymean = 0.5*(ymin + ymax)
!
! compute area and center of gravity
! do loop over all line pieces of polygon
!
area = 0.
sumx = 0.
sumy = 0.
!
do ipol=1,npol
  ipol1 = ipol + 1
  if(ipol==npol) ipol1 = 1
  xa = xpol(ipol)
  ya = ypol(ipol)
  xb = xpol(ipol1)
  yb = ypol(ipol1)
!
  darea = 0.5*((xa-xmean)*(yb-ymean) - (xb-xmean)*(ya-ymean))
  area  = area + darea
  sumx  = sumx + darea*(xa+xb+xmean)/3.
  sumy  = sumy + darea*(ya+yb+ymean)/3.
end do
!
return
end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!--------------------------------!!!!!
!!-------------Z-INTP1------------------!!
!!!!!--------------------------------!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


subroutine z_intp1(x1,y1,x2,y2,n1,n2,ierr)                                    !
!-----------------------------------------------------------------------------!
!
!   +-------+    ALKYON Hydraulic Consultancy & Research
!   |       |    Gerbrant van Vledder
!   |   +---+
!   |   | +---+
!   +---+ |   |
!         +---+
!
implicit none
!
!  0. Update history
!
!     30/03/1999 Initical version
!      9/04/1999 Check included for monotonicity of x-data
!     11/10/1999 Error messages added and updated
!     18/01/2001 Check include if n1==1
!     24/01/2001 Check for equality of y2 data loosened if n2==1
!     13/09/2001 Documentation updated
!
!  1. Purpose
!
!     Interpolate function values
!
!  2. Method
!
!     Linear interpolation

!     If a requested point falls outside the input domain, then
!     the nearest point is used (viz. begin or end point of x1/y1 array
!
!     If the input array has only one point. A constant value is assumed
!
!  3. Parameter list
!
!     Name    I/O  Type  Description
!
integer, intent(in) ::  n1   !   number of data points in x1-y1 arrays
integer, intent(in) ::  n2   !   number of data points in x2-y2 arrays
real, intent(in) ::  x1(n1)  !   x-values of input data
real, intent(in) ::  y1(n1)  !   y-values of input data
real, intent(in) ::  x2(n2)  !   x-values of output data
real, intent(out) :: y2(n2)  !   y-values of output data
integer, intent(out) :: ierr !   Error indicator
!
!  4. Subroutines used
!
!  5. Error messages
!
!     ierr = 0    No errors detected
!          = 1    x1-data not monotonic increasing
!          = 10   x2-data not monotonic increasing
!          = 11   x1- and x2 data not monotonic increasing
!          = 2    x1-data not monotonic decreasing
!          = 20   x1-data not monotonic decreasing
!          = 22   x1- and x2 data not monotonic decreasing
!
!          = 2    No variation in x1-data
!          = 3    No variation in x2-data is allowed if n2=1
!
!  6. Remarks
!
!     It is assumed that the x1- and x2-data are either
!     monotonic increasing or decreasing
!
!     If a requested x2-value falls outside the range of x1-values
!     it is assumed that the corresponding y2-value is equal to
!     the nearest boundary value of the y1-values
!
!     Example: x1 = [0 1 2 3]
!              y1 = [1 2 1 0]
!
!              x2 = -1,  y2 = 1
!              x2 =  5,  y2 = 0
!
!------------------------------------------------------------------------------
integer i1,i2        ! counters
!
real ds            ! step size
real fac           ! factor in linear interpolation
real s1,s2         ! search values
real xmin1,xmax1   ! minimum and maximum of x1-data
real xmin2,xmax2   ! minimum and maximum of x2-data
!
real, parameter :: eps=1.e-20
!------------------------------------------------------------------------------
!   initialisation
!
ierr = 0
!
!  check number of points of input array
!
if(n1==1) then
  y2 = y1(1)
  goto 9999
end if
!
!  check minimum and maximum data values
!
xmin1 = minval(x1)
xmax1 = maxval(x1)
xmin2 = minval(x2)
xmax2 = maxval(x2)
!
if (abs(xmin1-xmax1) < eps .or. abs(x1(1)-x1(n1)) < eps) then
  ierr = 2
  goto 9999
end if
!
if ((abs(xmin2-xmax2) < eps .or. abs(x2(1)-x2(n2)) < eps) .and. n2 > 1) then
  ierr = 3
  goto 9999
end if
!
! check input data for monotonicity
!
if(x1(1) < x1(n1)) then             ! data increasing
  do i1=1,n1-1
    if(x1(i1) > x1(i1+1)) then
      ierr=1
      write(*,*) 'z_intp1: i1 x1(i1) x1(i1+1):',i1,x1(i1),x1(i1+1)
      goto 9999
    end if
  end do
!
  do i2=1,n2-1
    if(x2(i2) > x2(i2+1)) then
      ierr=ierr+10
      write(*,*) 'z_intp1: i2 x2(i2) x2(i2+1):',i2,x2(i2),x2(i2+1)
      goto 9999
    end if
  end do
!
else                                 ! data decreasing
  do i1=1,n1-1
    if(x1(i1) < x1(i1+1)) then
      ierr=2
      write(*,*) 'z_intp1: i1 x1(i1) x1(i1+1):',i1,x1(i1),x1(i1+1)
      goto 9999
    end if
  end do
!
  do i2=1,n2-1
    if(x2(i2) < x2(i2+1)) then
      ierr=ierr + 20
      write(*,*) 'z_intp1: i2 x2(i2) x2(i2+1):',i2,x2(i2),x2(i2+1)
      goto 9999
    end if
  end do
end if
!
!------------------------------------------------------------------------------
! initialize
!------------------------------------------------------------------------------
if(ierr==0) then
  i1 = 1
  s1 = x1(i1)
!
  do i2 = 1,n2
    s2 = x2(i2)
    do while (s1 <= s2 .and. i1 < n1)
      i1 = i1 + 1
      s1 = x1(i1)
    end do
!
!  special point
!  choose lowest s1-value if x2(:) < x1(1)
!
    if(i1 ==1) then
      y2(i2) = y1(i1)
    else
      ds = s2 - x1(i1-1)
      fac = ds/(x1(i1)-x1(i1-1))
      y2(i2) = y1(i1-1) + fac*(y1(i1)-y1(i1-1))
    end if
!
! special case at end: choose s2(n2) > s1(n1), choose last value of y1(1)
!
    if(i2==n2 .and. s2>s1) y2(n2) = y1(n1)
  end do
end if
!
9999 continue
!
return
end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!----Y--SUBROUTINES---------------------!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!--------------------------------!!!!!
!!-------------Y-GAULEG-----------------!!
!!!!!--------------------------------!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


SUBROUTINE y_gauleg(x1,x2,x,w,n)
!-------------------------------------------------------------------
INTEGER, intent(in) ::  n    ! Number of intervals
real, intent(in)    ::  x1   ! lower limit of integration interval
real, intent(in)    ::  x2   ! upper limit of integration interval
real, intent(out)   ::  x(n) ! Positions for function evaluations
real, intent(out)   ::  w(n) ! Weights
!
!-----------------------------------------------------------------------
DOUBLE PRECISION EPS
PARAMETER (EPS=3.d-14)
INTEGER i,j,m
DOUBLE PRECISION p1,p2,p3,pp,xl,xm,z,z1
!-----------------------------------------------------------------------
      m=(n+1)/2
      xm=0.5d0*(x2+x1)
      xl=0.5d0*(x2-x1)
      do 12 i=1,m
        z=cos(3.141592654d0*(i-.25d0)/(n+.5d0))
1       continue
          p1=1.d0
          p2=0.d0
          do 11 j=1,n
            p3=p2
            p2=p1
            p1=((2.d0*j-1.d0)*z*p2-(j-1.d0)*p3)/j
11        continue
          pp=n*(z*p1-p2)/(z*z-1.d0)
          z1=z

          z=z1-p1/pp
        if(abs(z-z1).gt.EPS)goto 1
        x(i)=xm-xl*z
        x(n+1-i)=xm+xl*z
        w(i)=2.d0*xl/((1.d0-z*z)*pp*pp)
        w(n+1-i)=w(i)
12    continue
!
return
END subroutine
