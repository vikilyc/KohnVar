MODULE Integrators
    use Parameters
    CONTAINS
        SUBROUTINE SimpsonWeights(Points, Weights, StartValue, EndValue, N)
            ! This subroutine is adopted from ePolyScat with some modifications.
            ! https://epolyscat.droppages.com
            ! ** 
            ! ** Compute the points and weights for a repeated Simpson's rule integration
            ! **
            IMPLICIT NONE
            REAL (KIND = XR), ALLOCATABLE, INTENT(OUT), DIMENSION(:) :: Points, Weights
            REAL (KIND = XR) , INTENT(IN) :: StartValue ! First point in the grid
            REAL (KIND = XR), INTENT(IN) :: EndValue ! Last point in the grid //spacing between points on an evenly spaced 
            INTEGER :: N ! number of grid points, must be an odd number

        ! ** Local Variables
            INTEGER :: iPoint
            REAL (KIND = XR) :: DeltaValue
            IF (ALLOCATED(Points)) DEALLOCATE(Points)
            IF (ALLOCATED(Weights)) DEALLOCATE(Weights)
            ALLOCATE (Weights(N))
            ALLOCATE (Points(N))

            IF (N > 1) THEN
                DeltaValue = (EndValue - StartValue) / REAL(N - 1, KIND = XR)
            ELSE
                DeltaValue = 0.0_XR
            END IF
            Points = [(iPoint * DeltaValue + StartValue, iPoint = 0, N-1)]
            IF (N > 1) THEN
                Weights(1) = DeltaValue / 3.0_XR
                Weights(N) = DeltaValue / 3.0_XR
                Weights(2:N-1:2) = 4.0_XR * DeltaValue / 3.0_XR
                Weights(3:N-2:2) = 2.0_XR * DeltaValue / 3.0_XR
            ELSE
                Weights(1) = 1.0_XR
            END IF
        END SUBROUTINE SimpsonWeights
END MODULE INTEGRATORS

