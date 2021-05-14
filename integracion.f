!     -------------------------------------------------------------------------
!
!     -------------------------------------------------------------------------
!
!     -------------------------------------------------------------------------
!     -------------------------------------------------------------------------
!     integration routine
!     -------------------------------------------------------------------------
!     -------------------------------------------------------------------------
!
!     ------------------------------------------------------------------
!     ------------------------------------------------------------------
      module parameters_quad8
        integer,  parameter :: isize = 2, rsize = 5, csize = 2
        real (8), parameter :: tol = 1.d-8
      end module parameters_quad8
!     ------------------------------------------------------------------
c     gauss quadrature weights and kronrod quadrature abscissae and weights
c     as evaluated with 80 decimal digit arithmetic by l. w. fullerton,
c     bell labs, nov. 1981.
c
      module gauss_kronrod15
      integer,  parameter :: nkronrod = 8
      real (8), dimension (nkronrod/2) :: wg
      real (8), dimension (nkronrod) :: xgk, wgk

      data wg  (1) / 0.1294849661 6886969327 0611432679 082d0 /
      data wg  (2) / 0.2797053914 8927666790 1467771423 780d0 /
      data wg  (3) / 0.3818300505 0511894495 0369775488 975d0 /
      data wg  (4) / 0.4179591836 7346938775 5102040816 327d0 /
c
      data xgk (1) / 0.9914553711 2081263920 6854697526 329d0 /
      data xgk (2) / 0.9491079123 4275852452 6189684047 851d0 /
      data xgk (3) / 0.8648644233 5976907278 9712788640 926d0 /
      data xgk (4) / 0.7415311855 9939443986 3864773280 788d0 /
      data xgk (5) / 0.5860872354 6769113029 4144838258 730d0 /
      data xgk (6) / 0.4058451513 7739716690 6606412076 961d0 /
      data xgk (7) / 0.2077849550 0789846760 0689403773 245d0 /
      data xgk (8) / 0.0000000000 0000000000 0000000000 000d0 /
c
      data wgk (1) / 0.0229353220 1052922496 3732008058 970d0 /
      data wgk (2) / 0.0630920926 2997855329 0700663189 204d0 /
      data wgk (3) / 0.1047900103 2225018383 9876322541 518d0 /
      data wgk (4) / 0.1406532597 1552591874 5189590510 238d0 /
      data wgk (5) / 0.1690047266 3926790282 6583426598 550d0 /
      data wgk (6) / 0.1903505780 6478540991 3256402421 014d0 /
      data wgk (7) / 0.2044329400 7529889241 4161999234 649d0 /
      data wgk (8) / 0.2094821410 8472782801 2999174891 714d0 /
      end module gauss_kronrod15
!     ----------------------------------------------------
!     ----------------------------------------------------
!     Recursive Gauss-kronrod 15 points rule
!
!     n is an integer parameter control
!     sxi is a real parameter
!     (more parameters can be added)
!
!     (a,b) is the range
!
!     integrando is the external function to integrate
!
      real (8) function quadgk8v (ivoid,rvoid,cvoid,a,b,integrando )
     +     result (out)
        use gauss_kronrod15 ! points and weights of the G-K rule
        use parameters_quad8 ! tolerance parameters
        implicit none
        !     input
        integer,  dimension (isize) :: ivoid
        real (8), dimension (rsize) :: rvoid
        complex (8), dimension (csize) :: cvoid
        real (8)    :: a, b
        complex (8) :: c1
        !     internal variables
        integer   :: fcnt, i
        real (8) :: bma, bpa, q
        !     external functions
        real (8), external :: integrando, gkstp8v
        !     code
        fcnt = 3
        bma  = ( b - a )/2.d0
        bpa  = ( b + a )/2.d0
        q    = wgk(nkronrod)
     +         *integrando(ivoid,rvoid,cvoid,bma*xgk(nkronrod)+bpa)
        do i = 1, nkronrod-1
           q = q + wgk(i)*integrando(ivoid,rvoid,cvoid, bma*xgk(i)+bpa)
     +           + wgk(i)*integrando(ivoid,rvoid,cvoid,-bma*xgk(i)+bpa)
        enddo
        q   = bma*q
        out = gkstp8v(ivoid,rvoid,cvoid,a,b,q,fcnt,integrando)
        return
      end function quadgk8v
!     ------------------------------------------------------------
      recursive real (8) function gkstp8v( ivoid, rvoid, cvoid,
     +     a, b, q0, fcnt, integrando )
     +     result (out)
        use parameters_quad8
        use gauss_kronrod15
        implicit none
        !     input
        integer  :: fcnt
        integer,  dimension (isize) :: ivoid
        real (8), dimension (rsize) :: rvoid
        complex (8), dimension (csize) :: cvoid
        real (8) :: a, b, q0
        !     internal variables
        integer  :: i
        real (8) :: c, cma, cpa, bmc, bpc, q, q1, q2
        !     external functions
        real (8), external :: integrando
        !     code
        if (fcnt.gt.1000000) then
           out = q0
        else
          c   = ( a + b )/2.d0
          cma = ( c - a )/2.d0
          cpa = ( c + a )/2.d0
          bmc = ( b - c )/2.d0
          bpc = ( b + c )/2.d0
          q1  = wgk(nkronrod)
     +        *integrando(ivoid,rvoid,cvoid,cma*xgk(nkronrod)+cpa)
          q2  = wgk(nkronrod)
     +        *integrando(ivoid,rvoid,cvoid,bmc*xgk(nkronrod)+bpc)
          do i = 1, nkronrod-1
             q1 = q1
     +           + wgk(i)*integrando(ivoid,rvoid,cvoid, cma*xgk(i)+cpa)
     +           + wgk(i)*integrando(ivoid,rvoid,cvoid,-cma*xgk(i)+cpa)
              q2 = q2
     +           + wgk(i)*integrando(ivoid,rvoid,cvoid, bmc*xgk(i)+bpc)
     +           + wgk(i)*integrando(ivoid,rvoid,cvoid,-bmc*xgk(i)+bpc)
          enddo
          q1   = cma*q1
          q2   = bmc*q2
          q    = q1 + q2
          fcnt = fcnt + 2
          if ((fcnt.eq.5).or.(abs(q-q0).gt.tol)) then
            q1 = gkstp8v (ivoid,rvoid,cvoid,a,c,q1,fcnt,integrando)
            q2 = gkstp8v (ivoid,rvoid,cvoid,c,b,q2,fcnt,integrando)
             q = q1 + q2
          endif
          out = q
        endif
        return
      end function gkstp8v
