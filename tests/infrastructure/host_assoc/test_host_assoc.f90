module host_assoc_ext
    use, intrinsic :: iso_fortran_env, only: R8 => real64
    implicit none
    type :: rec_t
        real(R8) :: a = -1.0_R8
    end type rec_t
contains
    !> Takes an array through an assumed-shape dummy (a descriptor) and leaves it unchanged.
    subroutine pass_through(x)
        real(R8), intent(inout) :: x(:)
        x(1) = x(1) + 0.0_R8
    end subroutine pass_through

    subroutine pass_through_rec(r)
        type(rec_t), intent(inout) :: r(:)
        r(1)%a = r(1)%a + 0.0_R8
    end subroutine pass_through_rec
end module host_assoc_ext

program test_host_assoc
    !
    ! An internal procedure must address the main program's arrays, not a copy of them.
    !
    ! With the Intel compilers' OpenMP option (local variables automatic), a main-program array of explicit shape that
    ! main hands to an assumed-shape dummy lives on main's stack, while an internal procedure that reads or writes it by
    ! host association addresses a static copy: the helper reads zeros (or the type's default initialization) and its
    ! writes never reach main. The test programs are built with static arrays for that reason (-auto-scalar,
    ! tests/CMakeLists.txt); this program fails when a test executable is built without it.
    !
    !   H1  the internal procedure reads the values main wrote
    !   H2  main reads the values the internal procedure wrote
    !   H3  as H1 for an array of a derived type with a default-initialized component
    !
    use host_assoc_ext, only: R8, rec_t, pass_through, pass_through_rec
    implicit none

    integer, parameter :: N = 10

    real(R8)    :: tab(N)
    type(rec_t) :: recs(N)
    integer     :: k, nfail, nbad

    nfail = 0
    do k = 1, N
        tab(k) = real(k, R8)
        recs(k)%a = real(k, R8)
    enddo
    call pass_through(tab)
    call pass_through_rec(recs)

    call report('H1 the internal procedure reads the values main wrote', misread_tab())

    call overwrite_tab()
    nbad = 0
    do k = 1, N
        if (tab(k) /= -real(k, R8)) nbad = nbad + 1
    enddo
    write(*,'(a,es10.3,a)') '        main reads tab(5) = ', tab(5), ' (the internal procedure wrote -5)'
    call report('H2 main reads the values the internal procedure wrote', nbad)

    call report('H3 the internal procedure reads a derived-type array main wrote', misread_recs())

    if (nfail == 0) then
        write(*,'(a)') '[PASS] test_host_assoc: internal procedures address the main program''s arrays.'
        stop 0
    else
        write(*,'(a,i0,a)') '[FAIL] test_host_assoc: ', nfail, ' check(s) failed.'
        stop 1
    endif

contains

    integer function misread_tab() result(nb)
        integer :: i
        nb = 0
        do i = 1, N
            if (tab(i) /= real(i, R8)) nb = nb + 1
        enddo
        write(*,'(a,es10.3,a)') '        the internal procedure reads tab(5) = ', tab(5), ' (main wrote 5)'
    end function misread_tab

    subroutine overwrite_tab()
        integer :: i
        do i = 1, N
            tab(i) = -real(i, R8)
        enddo
    end subroutine overwrite_tab

    integer function misread_recs() result(nb)
        integer :: i
        nb = 0
        do i = 1, N
            if (recs(i)%a /= real(i, R8)) nb = nb + 1
        enddo
        write(*,'(a,es10.3,a)') '        the internal procedure reads recs(5)%a = ', recs(5)%a, ' (main wrote 5)'
    end function misread_recs

    subroutine report(what, nb)
        character(len=*), intent(in) :: what
        integer,          intent(in) :: nb
        if (nb == 0) then
            write(*,'(a,a)') '[PASS] ', what
        else
            write(*,'(a,a,a,i0,a,i0,a)') '[FAIL] ', what, ': ', nb, ' of ', N, ' values differ'
            nfail = nfail + 1
        endif
    end subroutine report

end program test_host_assoc
