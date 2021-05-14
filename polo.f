!     --------------------------------------------------------------------------
!     Pole search routine
!     -------------------------------------------------------------------------
      subroutine csearch ( shoot, gwa, ka, kb, ierr, irs )
!     Local search for complex zero of a complex function.
!     input: ka = starting value !will be overwritten!!!
!     kb = value in the vicinity of ka !will be overwritten
!     output: ka= zero of the function shoot
!     gwa=value of shoot at ka !(to control performance)
      implicit none
!     parameters
      integer,  parameter :: niterations = 1000 ! maximum number of iterations
      real (8), parameter :: e_tol = 0.0005d0 ! Use with determinant search
!     real (8), parameter :: e_tol = 500.d0 ! Use with 1/T search
      ! printing
      logical, parameter :: lprint = .false.
      ! input
      integer :: irs
      complex (8) :: gwb, kb
      ! input/output
      complex (8) :: gwa, ka
      ! output
      integer     :: ierr
      ! internal variables
      integer     :: j
      real (8)    :: gwamag, gwbmag, gw0mag
      complex (8) :: gwc, kc, gw0, k0
      ! external functions
      complex (8), external :: shoot
      ! code
      gwa    = shoot (ka,irs)
      gwamag =  dble(gwa*dconjg(gwa))
      gwamag = dsqrt(gwamag)
      !
      gwb    = shoot (kb,irs)
      gwbmag =  dble(gwb*dconjg(gwb))
      gwbmag = dsqrt(gwbmag)
      j      = -1
      !
      ! question:  gwa > gwb  ?
      if ( gwamag.lt.gwbmag ) then
        gwc = gwa
        kc  = ka
        gwa = gwb
        ka  = kb
        gwb = gwc
        kb  = kc
      endif
      ! gwa > gwb  : regula falsi
      !
      do j = 1, niterations ! the local search starts.
        ierr = j
        if (gwbmag.le.e_tol) then
          ka  = kb
          gwa = gwb
          exit
        end if
        ! question: bad convergence??? ( if too slow, stop! )
        gw0    = gwa - gwb
        gw0mag =  dble(gw0*DCONJG(gw0))
        gw0mag = dsqrt(gw0mag)
        if (gw0mag.le.e_tol) then
          ka  = kb
          gwa = gwb
          if (lprint) write (3,'(e12.4,a20)') gw0mag,'=gwb-gwa, STOP'
          exit
        end if
        ! convergence is fast enough, we proceed with a linear approximation
        ! to the complex function and guess the next argument.
        k0     = kb - gwb*(kb-ka)/(gwb-gwa)
        gw0    = shoot(k0,irs)
        gw0mag = dble(gw0*DCONJG(gw0))
        gw0mag = dsqrt(gw0mag)
        if (lprint) write (1124,'(a14,i3,3f10.0)')
     c        'Iteration No. ', j, k0, gw0mag
        ! the search continues. Reorder arguments.
        if (gw0mag.lt.gwbmag) then
          gwamag = gwbmag
          gwa    = gwb
          ka     = kb
          gwbmag = gw0mag
          gwb    = gw0
          kb     = k0
        else
          gwamag = gw0mag
          gwa    = gw0
          ka     = k0
        end if
        if (lprint) write(6,'(i3,5e16.8)') j, kb, gwb, gwbmag
      end do                    ! j
      return
      end subroutine csearch

!     ---------------------------------------------------------------
!     Aquí comienza la busquesa de polos en las hojas de Riemann
!     ---------------------------------------------------------------
      subroutine poles_search(irs,iteraciones, ierr,
     +     valorabsoluto, s1, s2, polo)
        implicit none
        ! input
        integer :: irs
        complex (8) :: s1, s2
        ! internal variables
      	integer :: irs1, ierr
      	complex (8) :: valorfuncion, uppervalue, lowervalue
        ! output
        integer  :: iteraciones
        real (8) :: valorabsoluto
        complex (8) :: polo
        ! external functions
        complex (8), external :: denominadorpolo
        ! code
        ierr = 0.d0
        irs1 = irs
        lowervalue = s1
        uppervalue = s2
        valorfuncion =  dcmplx(0.d0,0.d0)
          call csearch ( denominadorpolo, valorfuncion,
     +      lowervalue, uppervalue, ierr, irs1 )
!csearch ( shoot, gwa, ka, kb, ierr, irs )
        iteraciones = ierr
        valorabsoluto = dble(valorfuncion*dconjg(valorfuncion))
        polo = lowervalue
        return
      end subroutine poles_search
!     ---------------------------------------------------------------
      complex (8) function denominadorpolo (spoint,irs) result (out)
        use constantes
        implicit none
        ! input
        integer :: irs
        complex (8) :: spoint
        ! internal variables
        complex (8) :: rho1, rho2, tau1, tau2, RT1, RT2
        complex (8) :: Pa, Pab
        complex (8) :: Ma, Mb,g1, g2, h1, h2
        complex (8) :: A, B, C
        ! output
        complex (8) :: salida
        ! external functions
        real (8),    external :: respacio
        complex (8), external :: phasespace, cespacio
        !variables
        ! 1 polo
        Ma = dcmplx(mLambda1**2,0.d0)
        g1 = dcmplx(g1Lambda1**2,0.d0)
        g2 = dcmplx(g2Lambda1**2,0.d0)
        ! 2 polos
        Mb = dcmplx(mLambda1**2,0.d0)
        h1 = dcmplx(g1Lambda2**2,0.d0)
        h2 = dcmplx(g1Lambda2**2,0.d0)
        !Coeficientes del polinomio para 1 polos en K
        Pa = Ma - spoint
        !Coeficientes del polinomio para 2 polos en K
        Pab = (Ma - spoint)*(Mb - spoint)
        A =  (Ma - spoint)*g1 + (Mb - spoint)*h1
        B =  (Ma - spoint)*h1 + (Mb - spoint)*h2
        C = g1*h2 + g2*h1 - 2*g1Lambda1*g1Lambda2*g2Lambda1*g2Lambda2
        ! code
        rho1 = phasespace(mpion, msigma, spoint)
        rho2 = phasespace (mkaon, mnucleon, spoint)
        tau1 = cespacio (mpion,msigma, spoint)
        tau2 = cespacio (mkaon,mnucleon, spoint)

        !cruze en las hojas de Riemann
        ! rho1 I = rho1 + 2*t1... etc...
        RT1  = (rho1 + 2*tau1)
        RT2  = (rho2 + 2*tau2)

        select case (irs)
        ! 1 polo en la matriz K
        ! Pol(s) + A' rho1 + B' rho2

        case (11)
        salida = Pa - rho1*g1 - rho2*g2
!        salida = 0.1d0
!PRUEBA DE POLO
        !salida = Pa
        case (12)
        salida = Pa - (rho1*g1) - (RT2*g2)
        case (13)
        salida = Pa -RT1*g1 - RT2*g2
        case (14)
        salida = Pa - RT1*g1 - rho2*g2
        ! 2 polos en la matriz K
        ! Pol(s) + A rho1 + B rho2 + C rho1rho2
        case (21)
        salida = Pab + A*rho1 + B*rho2 + C*rho1*rho2
        case (22)
        salida = Pab + A*rho1 + B*RT2 + C*rho1*RT2
        case (23)
        salida = Pab + A*RT1 + B*RT2 + C*RT1*RT2
        case (24)
        salida = Pab + A*RT1 + B*rho2 + C*RT1*rho2
        case default
          stop 'wrong RS'
        end select
        out = salida
        return
      end function denominadorpolo
