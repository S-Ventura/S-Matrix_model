!     --------------------------------------------------------------

      include 'integracion.f'
      include 'subrutinas_matrices.f'
      include 'phasespacefunctions.f'
      include 'polo.f'

!     ---------------------------------------------------------------
!     Main program
!     ---------------------------------------------------------------
      program amplitud
        implicit none
!        call plotter_espaciodefases
!        call plotter3d_espaciodefases
        call plotter_amplitud
!        call poles_search
      end program amplitud
