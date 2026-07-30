module IGLOO_VectorModule
    use, intrinsic :: iso_fortran_env, only : I4 => int32, R8 => real64
    implicit none

contains

    !> CROSS PRODUCT BETWEEN VECTORS
    pure function normal2D(a) result(c)
        real(R8), intent(in) :: a(3)
        real(R8) :: c(3)
        c(1) = - a(2)
        c(2) =   a(1)
        c(3) = 0._R8
    end function normal2D


    !> CROSS PRODUCT BETWEEN VECTORS
    pure function cross(a, b) result(c)
        real(R8), intent(in) :: a(3), b(3)
        real(R8) :: c(3)
        c(1) = a(2)*b(3) - a(3)*b(2)
        c(2) = a(3)*b(1) - a(1)*b(3)
        c(3) = a(1)*b(2) - a(2)*b(1)
    end function cross


    !> FUNCTION TO ROTATE A VECTOR GIVEN AXIS & ANGLE
    pure function rotateVector(vect,axis,theta) result(rotatedVect)
        implicit none
        real(R8), intent(in) :: vect(3)
        real(R8), intent(in) :: axis(3)
        real(R8), intent(in) :: theta
        real(R8) :: cost, sint
        real(R8), dimension(3, 3) :: rotation
        real(R8) :: rotatedVect(3)

        cost = cos(theta)
        sint = sin(theta)
        rotation = reshape([ &
            axis(1)**2   *(1-cost)+        cost, axis(1)*axis(2)*(1-cost)-axis(3)*sint, axis(1)*axis(3)*(1-cost)+axis(2)*sint, &
            axis(2)*axis(1)*(1-cost)+axis(3)*sint, axis(2)**2   *(1-cost)+        cost, axis(2)*axis(3)*(1-cost)-axis(1)*sint, &
            axis(3)*axis(1)*(1-cost)-axis(2)*sint, axis(3)*axis(2)*(1-cost)+axis(1)*sint, axis(3)**2   *(1-cost)+        cost  &
            ], shape=[3, 3])
        rotatedVect = matmul(rotation, vect)

    end function rotateVector


    !> FUNCTION TO CHOOSE THE POINT WITH MINIMUM x, THEN MINIMUM y, THEN MINIMUM z
    pure integer function chooseVector(vectors)
        implicit none
        real(R8), intent(in) :: vectors(3,4)
        integer :: i, iMax

        iMax = 1
        do i = 2, 4
            if (compare_vectors(vectors(:,i), vectors(:,iMax))) then
                iMax = i
            end if
        end do
        chooseVector = iMax

    contains

        pure logical function compare_vectors(v1, v2)
        implicit none
        real(R8), intent(in) :: v1(3), v2(3)

        compare_vectors = (v1(1) < v2(1)) .or. &
                            ((v1(1) == v2(1)) .and. (v1(2) < v2(2))) .or. &
                            ((v1(1) == v2(1)) .and. (v1(2) == v2(2)) .and. (v1(3) < v2(3)))
        end function compare_vectors

    end function chooseVector

end module IGLOO_VectorModule
