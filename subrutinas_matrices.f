!     ---------------------------------------------------------------
      module dimensiones
        integer,  parameter :: matsize  = 2
        integer,  parameter :: masasize = 3
        integer,  parameter :: iosize   = 4
      end module dimensiones
!     ---------------------------------------------------------------
      module constantes
        integer,  parameter :: izero = 0
        real (8), parameter :: pi = dacos(-1.d0)
        real (8), parameter :: epsilon = 0.0001d0
        real (8), parameter :: epsilonhigh = 0.00001d0
        real (8), parameter :: mpion = 0.139d0
        real (8), parameter :: mkaon = 0.498d0
        real (8), parameter :: msigma = 1.192d0
        real (8), parameter :: mnucleon = 0.939d0
        real (8), parameter :: mnorm = 1.405d0
        complex (8), parameter :: xzero = dcmplx(0.d0,0.d0)
        complex (8), parameter :: xr = dcmplx(1.d0,0.d0)
        complex (8), parameter :: xi = dcmplx(0.d0,1.d0)
!     ---------------------------------------------------------------

!Parámetros de la matriz K
        real (8), parameter :: mLambda1 = 1.4d0
        real (8), parameter :: mLambda2 = 1.25d0
        real (8), parameter :: mLambda3 = 1.432d0
        real (8), parameter :: g1Lambda1 = 0.25d0, g2Lambda1 = 0.3d0
        real (8), parameter :: g1Lambda2 = 0.55d0, g2Lambda2 = 0.35d0
        real (8), parameter :: g1Lambda3 = 0.75d0, g2Lambda3 = 0.45d0
      end module constantes
!     ---------------------------------------------------------------
!     ---------------------------------------------------------------
!     Subrutinas que definen matrices
!     ---------------------------------------------------------------
!     ---------------------------------------------------------------
      subroutine matrizT (indice,s,masas,acoplos,Tmatriz)
        use dimensiones
        use constantes
        implicit none
        ! input
        integer  :: indice
        complex (8) :: s
        real (8), dimension (masasize) :: masas
        real (8), dimension (masasize,matsize,matsize) :: acoplos
        ! internal variables
        complex (8), dimension (matsize,matsize) :: kmatriz
        complex (8), dimension (matsize,matsize) :: rho
        complex (8), dimension (matsize,matsize) :: Tdenom, Tdenominv
        ! output
        complex (8), dimension (matsize,matsize) :: Tmatriz
        ! code
        call kmatrix (indice,s,masas,acoplos,kmatriz)
        call rhomatriz (s,rho)
        call Tdenominador (rho,kmatriz,Tdenom)
        call invertirmatriz(Tdenom,Tdenominv)
        Tmatriz = matmul(Tdenominv,kmatriz)
        return
      end subroutine matrizT
!     ---------------------------------------------------------------
      subroutine rhomatriz (s,rho)
        use dimensiones
        use constantes
        implicit none
        ! input
        complex (8) :: s
        ! output
        complex (8), dimension (matsize,matsize) :: rho
        ! external functions
        complex (8) :: phasespace
        ! code
        rho (1,1) = phasespace (mpion,msigma,s)
        rho (1,2) = xzero
        rho (2,1) = xzero
        rho (2,2) = phasespace (mkaon,mnucleon,s)
        return
      end subroutine rhomatriz
!     ---------------------------------------------------------------
      subroutine Tdenominador (rhomatriz,kmatriz,Tdenom)
        use dimensiones
        use constantes
        implicit none
        ! input
        complex (8), dimension (matsize,matsize) :: kmatriz
        ! internal variables
        complex (8), dimension (matsize,matsize) :: identidad, rhomatriz
        ! output
        complex (8), dimension (matsize,matsize) :: Tdenom
        ! code
        identidad = xzero
        identidad (1,1) = xr
        identidad (2,2) = xr
        Tdenom = identidad - matmul(rhomatriz,kmatriz)
        return
      end subroutine Tdenominador
!     ---------------------------------------------------------------
      subroutine invertirmatriz (A,B)
        use dimensiones
        implicit none
        ! input
        complex (8), dimension (matsize,matsize) :: A
        ! internal variables
        complex (8) :: det
        ! output
        complex (8), dimension (matsize,matsize) :: B
        ! code
        det    =   A(1,1)*A(2,2) - A(1,2)*A(2,1)
        B(1,1) =   A(2,2)/det
        B(1,2) = - A(1,2)/det
        B(2,1) = - A(2,1)/det
        B(2,2) =   A(1,1)/det
        return
      end subroutine invertirmatriz
!     ---------------------------------------------------------------
      subroutine kmatrix(indice,s,masas,acoplos,kmatriz)
        use dimensiones
        use constantes
        implicit none
        ! input
        integer  :: indice
        complex (8) :: s
        real (8), dimension (masasize) :: masas
        real (8), dimension (masasize,matsize,matsize) :: acoplos
        ! internal variables
        integer  :: i, j, k, ikmat
        real (8) :: norma
        complex (8) :: kterm
        ! output
        complex (8), dimension (matsize,matsize) :: kmatriz
        ! external functions
        real (8), external :: respacio
        ! code
        ! input k matrix
        select case (ikmat)
        case (1)
        ikmat = 1
        masas (ikmat) = mLambda1
        acoplos (ikmat,1,1) = g1Lambda1
        acoplos (ikmat,2,2) = g2Lambda1
        acoplos (ikmat,1,2) = g1Lambda1
        acoplos (ikmat,2,1) = g2Lambda1
        case (2)
        ikmat = 2
        masas (ikmat) = mLambda1
        acoplos (ikmat,1,1) = g1Lambda1
        acoplos (ikmat,2,2) = g2Lambda1
        acoplos (ikmat,1,2) = g1Lambda1
        acoplos (ikmat,2,1) = g2Lambda1
        case(3)
        ikmat = 3
        masas (ikmat) = mLambda1
        acoplos (ikmat,2,2) = g2Lambda1
        acoplos (ikmat,1,2) = g1Lambda1
        acoplos (ikmat,2,1) = g2Lambda1
        case default
          stop 'wrong RS'
        end select
        norma = respacio(mpion,msigma,mnorm*mnorm)/4.d0
        do i = 1, matsize
          do j = 1, matsize
            kmatriz (i,j) = 0.d0
            do k = 1, masasize
              kterm = acoplos(k,i,j)*acoplos(k,j,i)/(masas(k)**2-s)
              kmatriz (i,j) = kmatriz (i,j) + kterm/norma
            enddo
          enddo
        enddo
        return
      end subroutine Kmatrix
!     ---------------------------------------------------------------
!     Subrutinas del programa pricipal
!     ---------------------------------------------------------------
      subroutine plotter_amplitud
        use dimensiones
        use constantes
        implicit none
        ! parameters
        real (8), parameter :: st1 = ( mpion + msigma)**2
        real (8), parameter :: st2 = ( mkaon + mnucleon)**2
        ! input kinematics
        integer,  parameter :: npointsx = 500, npointsy = 50
        real (8), parameter :: sini= 1.d0, sfin = 4.d0
        real (8), parameter :: simini= -1.d0 , simfin = 1.d0
        ! input kmatrix
        integer :: indice
        real (8), dimension (masasize) :: masas
        ! internal variables
        integer  :: i, j
        real (8) :: stepre, stepim
        complex (8) :: s
        real (8), dimension (masasize,matsize,matsize) :: acoplos
        complex (8), dimension (matsize,matsize) :: rho
        complex (8), dimension (matsize,matsize) :: Tmatriz
        ! output
        real (8) :: sreal, simag
        complex (8) :: theory11p, theory22p, theory12p
        complex (8) :: theory11m, theory22m, theory12m
        ! input k matrix
        ! code
        indice = 1
        stepre = ( sfin - sini )/dble(npointsx-1)
        stepim = ( simfin - simini )/dble(npointsy-1)
        open (1,file='espaciodefases3D1.txt')
        sreal = sini
        do i = 1, npointsx
          simag = simini
          do j = 1, npointsy
            s = sreal + xi*simag
            call matrizT (indice,s,masas,acoplos,Tmatriz)
            theory11p = Tmatriz (1,1)
            theory12p = Tmatriz (1,2)
            theory22p = Tmatriz (2,2)
            s = sreal - xi*simag
            call matrizT (indice,s,masas,acoplos,Tmatriz)
            theory11m = Tmatriz (1,1)
            theory12m = Tmatriz (1,2)
            theory22m = Tmatriz (2,2)
            write (1,'(2F10.5,10F15.8)') sreal, simag,
     +         realpart(theory11p), imagpart(theory11p),
     +         realpart(theory11m), imagpart(theory11m),
     +         realpart(theory22p), imagpart(theory22p),
     +         realpart(theory22m), imagpart(theory22m),
     +         realpart(theory12p), imagpart(theory12m)
            simag = simag + stepim
          enddo
          sreal = sreal + stepre
        enddo
        close (1)
        print*, mLambda1, g1Lambda1,g2Lambda1
        return
        return
      end subroutine plotter_amplitud
!     ---------------------------------------------------------------
      subroutine plotter3d_espaciodefases
        use dimensiones
        use constantes
        implicit none
        ! parameters
        real (8), parameter :: st1 = ( mpion + msigma)**2
        real (8), parameter :: st2 = ( mkaon + mnucleon)**2
        ! input kinematics
        integer,  parameter :: npointsx = 300, npointsy = 50
        real (8), parameter :: sini=(1.3d0)**2 , sfin =(1.6d0)**2
        real (8), parameter :: simini=-1.d0 , simfin = 1.d0
        ! internal variables
        integer  :: i, j
        real (8) :: stepre, stepim
        complex (8) :: s
        complex (8), dimension (matsize,matsize) :: rho
        ! output
        real (8) :: sreal, simag
        complex (8) :: theory11p, theory22p, theory11m, theory22m
        ! code
        stepre = ( sfin - sini )/dble(npointsx-1)
        stepim = ( simfin - simini )/dble(npointsy-1)
        open (1,file='amplitud3D.txt')
        sreal = sini
        do i = 1, npointsx
          simag = simini
          do j = 1, npointsy
            s = sreal + xi*simag
            call rhomatriz (s,rho)
            theory11p = rho (1,1)
            theory22p = rho (2,2)
            s = sreal - xi*simag
            call rhomatriz (s,rho)
            theory11m = rho (1,1)
            theory22m = rho (2,2)
            write (1,'(2F10.5,8F15.8)') sreal, simag,
     +         realpart(theory11p), imagpart(theory11p),
     +         realpart(theory11m), imagpart(theory11m),
     +         realpart(theory22p), imagpart(theory22p),
     +         realpart(theory22m), imagpart(theory22m)
            simag = simag + stepim
          enddo
          sreal = sreal + stepre
        enddo
        close (1)
        return
      end subroutine plotter3d_espaciodefases
!     ---------------------------------------------------------------
      subroutine plotter_espaciodefases
        use dimensiones
        use constantes
        implicit none
        ! parameters
        real (8), parameter :: st1 = ( mpion + msigma)**2
        real (8), parameter :: st2 = ( mkaon + mnucleon)**2
        ! input
        integer,  parameter :: npoints = 1000
        real (8), parameter :: sini=(1.3d0)**2 , sfin =(1.6d0)**2
        ! internal variables
        real (8) :: step
        complex (8) :: s
        complex (8), dimension (matsize,matsize) :: rho
        ! output
        integer  :: i
        real (8) :: sreal, theory11, theory22
        complex (8) :: theory11p, theory22p, theory11m, theory22m
        ! external functions
        real (8), external :: respacio
        ! code
        step = ( sfin - sini )/dble(npoints-1)
        open (1,file='espaciodefases.txt')
        sreal = sini
        do i = 1, npoints
          if (sreal.gt.st1) then
            theory11 = respacio(mpion,msigma,sreal)
          else
            theory11 = 0.d0
          endif
          if (sreal.gt.st2) then
            theory22 = respacio(mkaon,mnucleon,sreal)
          else
            theory22 = 0.d0
          endif
          s = sreal + xi*epsilon
          call rhomatriz (s,rho)
          theory11p = rho (1,1)
          theory22p = rho (2,2)
          s = sreal - xi*epsilon
          call rhomatriz (s,rho)
          theory11m = rho (1,1)
          theory22m = rho (2,2)
          write (1,'(I7,F10.5,10F15.8)') i, sreal,
     +         theory11, realpart(theory11p), imagpart(theory11p),
     +                   realpart(theory11m), imagpart(theory11m),
     +         theory22, realpart(theory22p), imagpart(theory22p),
     +                   realpart(theory22m), imagpart(theory22m)
          sreal = sreal + step
        enddo
        close (1)
        return
      end subroutine plotter_espaciodefases
!     ---------------------------------------------------------------
