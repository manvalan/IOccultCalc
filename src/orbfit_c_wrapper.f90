!===============================================================================
! orbfit_c_wrapper.f90
!
! Wrapper Fortran per esporre funzioni OrbFit con interfaccia C-compatible
! Questo permette di chiamare le routine OrbFit da C++
!
! Author: IOccultCalc project
! Date: November 2025
! License: GPL (same as OrbFit)
!===============================================================================

module orbfit_c_interface
  use iso_c_binding
  use propag_state, only: inipro
  ! use prelim  ! Commentato per ora - gaussn e vaisala hanno problemi con common blocks
  implicit none

  ! Dichiarazioni esplicite per funzioni prelim (con common blocks)
  interface
    subroutine gaussn(tobs, alpha, delta, obscod, elem, eletyp, &
                      t0, nroots, nsol, rtop, fail, msg, debug, multi)
      double precision, intent(in) :: tobs(3), alpha(3), delta(3)
      integer, intent(in) :: obscod(3)
      double precision, intent(out) :: elem(6,3), t0(3), rtop(3)
      character(len=*), intent(out) :: eletyp(3), msg
      logical, intent(in) :: debug, multi
      logical, intent(out) :: fail
      integer, intent(out) :: nroots, nsol
    end subroutine gaussn
    
    subroutine vaisala(tobs, alpha, delta, obscod, solar_r, &
                       elem, eletyp, t0, rtop, fail, msg, debug)
      double precision, intent(in) :: tobs(2), alpha(2), delta(2), solar_r
      integer, intent(in) :: obscod(2)
      double precision, intent(out) :: elem(6), t0, rtop
      character(len=*), intent(out) :: eletyp, msg
      logical, intent(in) :: debug
      logical, intent(out) :: fail
    end subroutine vaisala
    
    subroutine iodini()
    end subroutine iodini
  end interface

contains

  !-----------------------------------------------------------------------------
  ! Gauss Method Wrapper
  !-----------------------------------------------------------------------------
  subroutine gauss_c_wrapper(tobs, alpha, delta, obscod, &
                              elem, eletyp, t0, nroots, nsol, rtop, &
                              fail, errmsg, debug) &
       bind(C, name="orbfit_gauss_c")
    
    ! Interfaccia C
    real(c_double), dimension(3), intent(in) :: tobs    ! MJD
    real(c_double), dimension(3), intent(in) :: alpha   ! rad
    real(c_double), dimension(3), intent(in) :: delta   ! rad
    integer(c_int), dimension(3), intent(in) :: obscod  ! observatory codes
    real(c_double), dimension(6,3), intent(out) :: elem ! orbital elements [6,3]
    character(kind=c_char), dimension(4,3), intent(out) :: eletyp ! element types
    real(c_double), dimension(3), intent(out) :: t0     ! epochs
    integer(c_int), intent(out) :: nroots               ! number of roots
    integer(c_int), intent(out) :: nsol                 ! number of solutions
    real(c_double), dimension(3), intent(out) :: rtop   ! topocentric distances
    integer(c_int), intent(out) :: fail                 ! failure flag
    character(kind=c_char), dimension(200), intent(out) :: errmsg ! error message
    integer(c_int), intent(in) :: debug                 ! debug flag
    
    ! Variabili Fortran interne
    double precision :: tobs_f90(3), alpha_f90(3), delta_f90(3)
    integer :: obscod_f90(3)
    double precision :: elem_f90(6,3), t0_f90(3), rtop_f90(3)
    character(len=4) :: eletyp_f90(3)
    character(len=200) :: msg_f90
    logical :: fail_f90, debug_f90, multi
    integer :: nroots_f90, nsol_f90
    integer :: i, j
    
    ! Converti input da C a Fortran
    tobs_f90 = tobs
    alpha_f90 = alpha
    delta_f90 = delta
    obscod_f90 = obscod
    debug_f90 = (debug /= 0)
    multi = .true.
    
    ! Inizializza output
    fail_f90 = .false.
    msg_f90 = ' '
    elem_f90 = 0.d0
    eletyp_f90 = '    '
    t0_f90 = 0.d0
    rtop_f90 = 0.d0
    nroots_f90 = 0
    nsol_f90 = 0
    
    ! Chiama la routine Gauss di OrbFit
    ! NOTA: Richiede che OrbFit sia stato inizializzato (iodini chiamato)
    call gaussn(tobs_f90, alpha_f90, delta_f90, obscod_f90, &
                elem_f90, eletyp_f90, t0_f90, nroots_f90, nsol_f90, &
                rtop_f90, fail_f90, msg_f90, debug_f90, multi)
    
    ! Converti output da Fortran a C
    elem = elem_f90
    t0 = t0_f90
    rtop = rtop_f90
    nroots = nroots_f90
    nsol = nsol_f90
    fail = merge(1, 0, fail_f90)
    
    ! Converti stringhe element type (Fortran → C)
    do i = 1, 3
      do j = 1, 4
        eletyp(j, i) = eletyp_f90(i)(j:j)
      end do
    end do
    
    ! Converti messaggio errore (Fortran → C)
    do i = 1, min(200, len_trim(msg_f90))
      errmsg(i) = msg_f90(i:i)
    end do
    ! Null-terminate la stringa C
    if (len_trim(msg_f90) < 200) then
      errmsg(len_trim(msg_f90) + 1) = c_null_char
    end if
    
  end subroutine gauss_c_wrapper

  !-----------------------------------------------------------------------------
  ! Vaisala Method Wrapper
  !-----------------------------------------------------------------------------
  subroutine vaisala_c_wrapper(tobs, alpha, delta, obscod, solar_r, &
                                elem, eletyp, t0, rtop, fail, errmsg, debug) &
       bind(C, name="orbfit_vaisala_c")
    
    ! Interfaccia C
    real(c_double), dimension(2), intent(in) :: tobs    ! MJD
    real(c_double), dimension(2), intent(in) :: alpha   ! rad
    real(c_double), dimension(2), intent(in) :: delta   ! rad
    integer(c_int), dimension(2), intent(in) :: obscod
    real(c_double), intent(in) :: solar_r               ! heliocentric distance guess
    real(c_double), dimension(6), intent(out) :: elem
    character(kind=c_char), dimension(4), intent(out) :: eletyp
    real(c_double), intent(out) :: t0
    real(c_double), intent(out) :: rtop
    integer(c_int), intent(out) :: fail
    character(kind=c_char), dimension(200), intent(out) :: errmsg
    integer(c_int), intent(in) :: debug
    
    ! Variabili Fortran
    double precision :: tobs_f90(2), alpha_f90(2), delta_f90(2)
    integer :: obscod_f90(2)
    double precision :: solar_r_f90
    double precision :: elem_f90(6), t0_f90, rtop_f90
    character(len=4) :: eletyp_f90
    character(len=200) :: msg_f90
    logical :: fail_f90, debug_f90
    integer :: i
    
    ! Converti input
    tobs_f90 = tobs
    alpha_f90 = alpha
    delta_f90 = delta
    obscod_f90 = obscod
    solar_r_f90 = solar_r
    debug_f90 = (debug /= 0)
    
    ! Inizializza
    fail_f90 = .false.
    msg_f90 = ' '
    
    ! Chiama Vaisala di OrbFit
    call vaisala(tobs_f90, alpha_f90, delta_f90, obscod_f90, solar_r_f90, &
                 elem_f90, eletyp_f90, t0_f90, rtop_f90, fail_f90, msg_f90, debug_f90)
    
    ! Converti output
    elem = elem_f90
    t0 = t0_f90
    rtop = rtop_f90
    fail = merge(1, 0, fail_f90)
    
    ! Converti element type string
    do i = 1, 4
      eletyp(i) = eletyp_f90(i:i)
    end do
    
    ! Converti error message
    do i = 1, min(200, len_trim(msg_f90))
      errmsg(i) = msg_f90(i:i)
    end do
    if (len_trim(msg_f90) < 200) then
      errmsg(len_trim(msg_f90) + 1) = c_null_char
    end if
    
  end subroutine vaisala_c_wrapper

  !-----------------------------------------------------------------------------
  ! RA15 Propagator Initialization
  !-----------------------------------------------------------------------------
  subroutine ra15_init_c_wrapper() bind(C, name="orbfit_ra15_init_c")
    
    ! Chiama inizializzazione propagatore OrbFit
    call inipro()
    
  end subroutine ra15_init_c_wrapper

  !-----------------------------------------------------------------------------
  ! RA15 Propagator Wrapper - DISABLED
  ! NOTE: propag() requires orbit_elem type which is complex to wrap
  ! For now we only expose Gauss and Vaisala methods
  !-----------------------------------------------------------------------------
  ! subroutine ra15_propag_c_wrapper(...) bind(C, name="orbfit_ra15_propag_c")
  ! DISABLED - requires orbit_elem type
  ! end subroutine

  !-----------------------------------------------------------------------------
  ! Initialization Wrapper (loads data files)
  !-----------------------------------------------------------------------------
  subroutine orbfit_init_c_wrapper(lib_path, lib_path_len) &
       bind(C, name="orbfit_init_c")
    
    ! Interfaccia C
    integer(c_int), intent(in), value :: lib_path_len
    character(kind=c_char), dimension(lib_path_len), intent(in) :: lib_path
    
    ! Variabili Fortran
    character(len=lib_path_len) :: lib_path_f90
    integer :: i
    
    ! Converti path da C a Fortran
    do i = 1, lib_path_len
      if (lib_path(i) == c_null_char) exit
      lib_path_f90(i:i) = lib_path(i)
    end do
    
    ! TODO: Set environment variables per OrbFit data path
    ! call setenv('ORBFIT_DATA', lib_path_f90, .true.)
    
    ! Inizializza IOD (Initial Orbit Determination)
    call iodini()
    
    ! Inizializza propagatore
    call inipro()
    
    print *, 'OrbFit initialized from Fortran wrapper'
    
  end subroutine orbfit_init_c_wrapper

end module orbfit_c_interface


!===============================================================================
! STUBS per funzioni OrbFit mancanti
!
! Se OrbFit non è linkato completamente, questi stub evitano errori di linking.
! Vanno rimossi quando si linka con le vere librerie OrbFit.
!===============================================================================

! NOTA: Commentare/rimuovere questa sezione quando si linkano le vere librerie

! subroutine gaussn(tobs, alpha, delta, obscod, elem, eletyp, &
!                   t0, nroots, nsol, rtop, fail, msg, debug, multi)
!   double precision, intent(in) :: tobs(3), alpha(3), delta(3)
!   integer, intent(in) :: obscod(3)
!   double precision, intent(out) :: elem(6,3), t0(3), rtop(3)
!   character(len=*), intent(out) :: eletyp(3), msg
!   logical, intent(in) :: debug, multi
!   logical, intent(out) :: fail
!   integer, intent(out) :: nroots, nsol
!   
!   fail = .true.
!   msg = 'gaussn: STUB - link with real OrbFit libgauss.a'
!   nroots = 0
!   nsol = 0
! end subroutine
! 
! subroutine vaisala(tobs, alpha, delta, obscod, solar_r, &
!                    elem, eletyp, t0, rtop, fail, msg, debug)
!   double precision, intent(in) :: tobs(2), alpha(2), delta(2), solar_r
!   integer, intent(in) :: obscod(2)
!   double precision, intent(out) :: elem(6), t0, rtop
!   character(len=*), intent(out) :: eletyp, msg
!   logical, intent(in) :: debug
!   logical, intent(out) :: fail
!   
!   fail = .true.
!   msg = 'vaisala: STUB - link with real OrbFit libgauss.a'
! end subroutine
! 
! subroutine inipro()
!   print *, 'inipro: STUB - link with real OrbFit libprop.a'
! end subroutine
! 
! subroutine propag(el0, t0, t1, el1, dum_mat, dum_mat2, nd)
!   double precision, intent(in) :: el0(6), t0, t1
!   double precision, intent(out) :: el1(6)
!   double precision, intent(in) :: dum_mat, dum_mat2
!   integer, intent(in) :: nd
!   
!   el1 = el0  ! Dummy: no propagation
!   print *, 'propag: STUB - link with real OrbFit libprop.a'
! end subroutine
! 
! subroutine iodini()
!   print *, 'iodini: STUB - link with real OrbFit libgauss.a'
! end subroutine
