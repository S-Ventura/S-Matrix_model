!     ---------------------------------------------------------------
!     Kallen funtion, phasespace functions, irho(s) integrando funtions
!     ---------------------------------------------------------------
      real (8) function rlambda (a,b,c) result (out)
        implicit none
        ! input
        real (8) :: a, b, c
        ! code
        out = a*a + b*b + c*c - 2.d0*( a*b + a*c + b*c )
        return
      end function rlambda
!     ---------------------------------------------------------------
      complex (8) function clambda (a,b,c) result (out)
        implicit none
        ! input
        complex (8) :: a, b, c
        ! code
        out = a*a + b*b + c*c - 2.d0*( a*b + a*c + b*c )
        return
      end function clambda
!     ---------------------------------------------------------------
!     phasespace functions tau(s) = {lambda^1/2}/2s
!     ---------------------------------------------------------------
      real (8) function respacio (m1,m2,s) result (out)
        use dimensiones
        implicit none
        ! input
        real (8) :: s, m1, m2
        ! internal variables
        integer  :: i
        real (8) :: l1, m1sq, m2sq
        ! external functions
        real (8), external :: rlambda
        ! code
        m1sq = m1*m1
        m2sq = m2*m2
        l1   = rlambda(m1sq,m2sq,s)
        if (l1.lt.0.d0) l1 = 0.d0
        out  = dsqrt(l1)/(4.d0*s)
        return
      end function respacio
!     ---------------------------------------------------------------
        complex (8) function cespacio (m1,m2,s) result (out)
        use dimensiones
        implicit none
        ! input
        real (8) :: m1, m2
        complex (8) :: s
        ! internal variables
        integer :: i
        complex (8) :: l1, m1sq, m2sq
        ! external functions
        complex (8), external :: clambda
        ! code
        m1sq = dcmplx(m1**2,0.d0)
        m2sq = dcmplx(m2**2,0.d0)
        l1   = clambda(m1sq,m2sq,s)
        out  = sqrt(l1)/(dcmplx(4.d0,0)*s)
        return
      end function cespacio
!     ---------------------------------------------------------------
!     Continuación analitíca del espacio de fases - Chew Mandelstam
!     ---------------------------------------------------------------
      complex (8) function phasespace (m1,m2,s) result (out)
        use dimensiones
        use constantes
        use parameters_quad8
        implicit none
        ! input
        real (8) :: m1, m2
        complex (8) :: s
        ! internal variables
        integer  :: i
        real (8) :: st, a, b
        real (8) :: partereal, parteimag
        integer,  dimension (isize) :: ivoid
        real (8), dimension (rsize) :: rvoid
        complex (8), dimension (csize) :: cvoid
        ! external functions
        real (8), external :: quadgk8v
        real (8), external :: integrando_re, integrando_im
        ! code
        ivoid = 0.d0
        rvoid = 0.d0
        cvoid = xzero
        rvoid (1) = m1
        rvoid (2) = m2
        cvoid (1) = s
        st = (m1+m2)**2 + epsilonhigh
        a = datan(st)
        b = pi/2.d0 - epsilonhigh
        partereal = quadgk8v(ivoid,rvoid,cvoid,a,b,integrando_re)
        parteimag = quadgk8v(ivoid,rvoid,cvoid,a,b,integrando_im)
        out = (s-st)*dcmplx(partereal,parteimag)/pi
        return
      end function phasespace
!     ---------------------------------------------------------------
      real (8) function integrando_re(ivoid,rvoid,cvoid,x) result (out)
        use dimensiones
        use parameters_quad8
        implicit none
        ! input
        real (8) :: x
        integer,  dimension (isize) :: ivoid
        real (8), dimension (rsize) :: rvoid
        complex (8), dimension (csize) :: cvoid
        ! internal variables
        integer  :: i
        real (8) :: sp, wp, st, m1, m2
        complex (8) :: s
        ! output
!        real (8) :: sal
        complex (8) :: salida
        ! external functions
        real (8), external :: respacio
        ! code
        m1 = rvoid (1)
        m2 = rvoid (2)
        s  = cvoid (1)
        st = ( m1+m2 )**2
        wp = 1.d0/cos(x)**2
        sp = dtan(x)
        salida = wp*respacio(m1,m2,sp)/((sp-s)*(sp-st))
!        sal = dble(salida)
!        if (isnan(sal)) print*,sp,respacio(m1,m2,sp),sp-s,sp-st
        out    = realpart(salida)
        return
      end function integrando_re
!     ---------------------------------------------------------------
      real (8) function integrando_im(ivoid,rvoid,cvoid,x) result (out)
        use dimensiones
        use parameters_quad8
        implicit none
        ! input
        real (8) :: x
        integer,  dimension (isize) :: ivoid
        real (8), dimension (rsize) :: rvoid
        complex (8), dimension (csize) :: cvoid
        ! internal variables
        integer  :: i
        real (8) :: sp, wp, st, m1, m2
        complex (8) :: s
        ! output
        complex (8) :: salida
        ! external functions
        real (8), external :: respacio
        ! code
        s  = cvoid (1)
        m1 = rvoid (1)
        m2 = rvoid (2)
        st = ( m1+m2 )**2
        wp = 1.d0/cos(x)**2
        sp = dtan(x)
        salida = wp*respacio(m1,m2,sp)/((sp-s)*(sp-st))
        out    = imagpart(salida)
        return
      end function integrando_im
!     ---------------------------------------------------------------
