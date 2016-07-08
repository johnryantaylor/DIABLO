      SUBROUTINE GHOST_CHAN_MPI
! This subroutine is part of the MPI package for the channel flow
! Diablo package.

      RETURN
      END

      SUBROUTINE END_RUN_MPI
! This subroutine is part of the MPI package for the channel flow
! Diablo package.

      RETURN
      END


      SUBROUTINE GHOST_GRID_MPI
! This subroutine is part of the MPI package for the channel flow
! Diablo package.

      RETURN
      END

C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE THOMAS_FORWARD_REAL_MPI(A,B,C,G,NY,NX)
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
C This subroutine performs the forward sweep of the Thomas algorithm

      RETURN
      END

C----*|--.---------.---------.---------.---------.---------.---------.-|------
      SUBROUTINE THOMAS_FORWARD_COMPLEX_MPI(A,B,C,G,NY,NX)
C----*|--.---------.---------.---------.---------.---------.---------.-|-----|
C This subroutine performs the forward sweep of the Thomas algorithm

      RETURN
      END

C----*|--.---------.---------.---------.---------.---------.---------.-|------
      SUBROUTINE THOMAS_BACKWARD_REAL_MPI(A,B,C,G,NY,NX)
C----*|--.---------.---------.---------.---------.---------.---------.-|------
C This subroutine performs the backward sweep of the Thomas algorithm

      RETURN
      END

C----*|--.---------.---------.---------.---------.---------.---------.-|------
      SUBROUTINE THOMAS_BACKWARD_COMPLEX_MPI(A,B,C,G,NY,NX)
C----*|--.---------.---------.---------.---------.---------.---------.-|-----|
C This subroutine performs the backward sweep of the Thomas algorithm

      RETURN
      END

!----*|--.---------.---------.---------.---------.---------.---------.-|------
      SUBROUTINE INIT_MPI
C----*|--.---------.---------.---------.---------.---------.---------.-|-----|

      RETURN
      END


!----*|--.---------.---------.---------.---------.---------.---------.-|------
      SUBROUTINE INIT_CHAN_MPI
C----*|--.---------.---------.---------.---------.---------.---------.-|-----|
C Initialize any constants here

      RETURN
      END

      SUBROUTINE APPLY_BC_TH_MPI(MATL,MATD,MATU,VEC,N)
! This subroutine applies the boundary conditions to the
! scalar fields prior to the implicit solve
      RETURN
      END

      SUBROUTINE APPLY_BC_TH_MPI_C(MATL_C,MATD_C,MATU_C,VEC_C,N)
! This subroutine applies the boundary conditions to the
! scalar fields prior to the implicit solve
      RETURN
      END

      SUBROUTINE APPLY_BC_U2_MPI(MATL,MATD,MATU,VEC)
! This subroutine applies the boundary conditions to the
        RETURN
        END

      SUBROUTINE APPLY_BC_U2_MPI_C(MATL_C,MATD_C,MATU_C,VEC_C)
! This subroutine applies the boundary conditions to the
        RETURN
        END

      SUBROUTINE APPLY_BC_U1_MPI(MATL,MATD,MATU,VEC)
! This subroutine applies the boundary conditions to the
! velocity field prior to the implicit solve
        RETURN
        END

      SUBROUTINE APPLY_BC_U1_MPI_C(MATL_C,MATD_C,MATU_C,VEC_C)
! This subroutine applies the boundary conditions to the
        RETURN
        END

      SUBROUTINE APPLY_BC_U3_MPI(MATL,MATD,MATU,VEC)
! This subroutine applies the boundary conditions to the
        RETURN
        END

      SUBROUTINE APPLY_BC_U3_MPI_C(MATL_C,MATD_C,MATU_C,VEC_C)
! This subroutine applies the boundary conditions to the
        RETURN
        END

      SUBROUTINE APPLY_BC_REM_DIV_MPI(MATL_C,MATD_C,MATU_C,VEC_C,K)
! This subroutine applies the boundary conditions for the Poisson Eq.
! Note, MATL, MATD, etc. are dimensioned in header

        RETURN
        END

      SUBROUTINE APPLY_BC_POISSON_MPI(MATL_C,MATD_C,MATU_C,VEC_C,K)
! This subroutine applies the boundary conditions for the Poisson Eq.
! Note, MATL, MATD, etc. are dimensioned in header

      RETURN
      END

      SUBROUTINE APPLY_BC_VEL_MPI
! This subroutine applies the boundary conditions for the Poisson Eq.
! Note, MATL, MATD, etc. are dimensioned in header
 
      RETURN
      END
 
      SUBROUTINE get_minimum_mpi(val)                                          
! This subroutine is part of the MPI package for the channel flow
! Diablo package. 

      RETURN   
      END    

