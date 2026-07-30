module verif_dump
    ! Curve dumps for the (non-gating) unit-family plot pipeline.
    ! Files land in curves/<family>/<id>.dat under the ctest working dir
    ! (build_verif) and are rendered by tools/plot_curves.py into one SVG per
    ! family. Column 1 is x; data column 1 is the PRODUCTION curve (drawn as
    ! markers), the rest are references (lines). Labels are |-separated
    ! (mathtext allowed). Best-effort: any I/O problem is reported and skipped,
    ! never failing the test.
    use, intrinsic :: iso_fortran_env, only: real64
    implicit none
    private

    public :: dump_curve

contains

    subroutine dump_curve(family, id, title, xlabel, ylabel, labels, x, ys, xlog, ylog)
        character(len=*), intent(in) :: family   ! e.g. 'standard-drag'
        character(len=*), intent(in) :: id       ! file stem, e.g. 'SD1'
        character(len=*), intent(in) :: title, xlabel, ylabel
        character(len=*), intent(in) :: labels   ! 'production|ref 1|ref 2'
        real(real64), intent(in) :: x(:)
        real(real64), intent(in) :: ys(:,:)      ! (size(x), ncurves)
        logical, intent(in), optional :: xlog, ylog
        integer :: u, i, ios
        logical :: lx, ly

        lx = .false.; ly = .false.
        if (present(xlog)) lx = xlog
        if (present(ylog)) ly = ylog
        if (size(ys, 1) /= size(x)) then
            write(*, '(a)') 'verif_dump: size mismatch for '//trim(id)//' (skipped)'
            return
        end if

        call execute_command_line('mkdir -p curves/'//trim(family), wait=.true.)
        open(newunit=u, file='curves/'//trim(family)//'/'//trim(id)//'.dat', &
             status='replace', action='write', iostat=ios)
        if (ios /= 0) then
            write(*, '(a)') 'verif_dump: cannot open dump for '//trim(id)//' (skipped)'
            return
        end if
        write(u, '(a)') '# id: '//trim(id)
        write(u, '(a)') '# title: '//trim(title)
        write(u, '(a)') '# xlabel: '//trim(xlabel)
        write(u, '(a)') '# ylabel: '//trim(ylabel)
        write(u, '(a)') '# labels: '//trim(labels)
        if (lx) write(u, '(a)') '# xlog: 1'
        if (ly) write(u, '(a)') '# ylog: 1'
        do i = 1, size(x)
            write(u, '(*(es16.8e3,1x))') x(i), ys(i, :)
        end do
        close(u)
    end subroutine dump_curve

end module verif_dump
