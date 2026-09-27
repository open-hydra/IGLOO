!> Integer-temperature property tables (cp, h, rho, sigma, mu, psat), their lookup, and the header,
!> row, node and column checks of INPUT/<prefix>properties.dat (the same checks as ICE's reader).
module IGLOO_Lib_Properties
  use, intrinsic :: iso_fortran_env, only: R8 => real64
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  private :: table_h_offset
  integer :: Tmin, Tmax    ! Extreme temperatures in tables
  real(kind=8), dimension(:,:), allocatable :: cp_tab, h_tab, rho_tab, sig_tab, mup_tab, psat_tab

  integer, parameter :: ntok_max = 32, tok_len = 64, nzone_max = 64

  !> Codes of the header, node and column checks; table_reason gives the text.
  integer, parameter :: TAB_OK = 0, TAB_NO_TEMPERATURE = 1, TAB_NO_CP = 2, TAB_NO_DENSITY = 3, &
                        TAB_NO_ENTHALPY = 4, TAB_TWO_ENTHALPY = 5, TAB_DUPLICATE = 6, &
                        TAB_FEW_ROWS = 7, TAB_OFF_NODE = 8, TAB_NONFINITE = 9, &
                        TAB_RHO_NONPOSITIVE = 10, TAB_CP_NONPOSITIVE = 11, TAB_H_NONMONOTONE = 12, &
                        TAB_H_CP_MISMATCH = 13, TAB_DATUM_MISMATCH = 14, TAB_NEGATIVE_T = 15

  !> Codes of the Psat column checks; psat_reason gives the text.
  integer, parameter :: PSAT_OK = 0, PSAT_ABSENT = 1, PSAT_NONFINITE = 2, PSAT_NEGATIVE = 3, &
                        PSAT_NONMONOTONE = 4, PSAT_CONSTANT = 5, PSAT_TBOIL_RANGE = 6, &
                        PSAT_TBOIL_MISMATCH = 7

  character(len=*), parameter :: grammar = &
    'expected: VARIABLES = "Temperature", "Cp", "Density", "Enthalpy" (or "Enthalpy_abs")[, "Psat"], '// &
    'Temperature first and the others in any order, one zone per material, rows on consecutive integer kelvins'

contains

  !> Linear interpolation in a 1D table indexed by integer temperature.
  pure function lookupTab(tab, T) result(val)
    implicit none
    real(8), intent(in) :: tab(Tmin:Tmax), T
    real(8) :: val
    real(8) :: Vij, Viij, Tdiff
    integer :: T_i

    !> Index clamped to the tabulated range.
    T_i   = min(max(idint(T), Tmin), Tmax-1)
    Tdiff = T - T_i
    Vij   = tab(T_i)       ! int(T)
    Viij  = tab(T_i + 1)   ! int(T)+1
    val   = Vij + (Viij-Vij)*Tdiff

  end function lookupTab

  !> Inverse of lookupTab on an increasing table: the end segments extended, the bracketing node by bisection.
  pure function comp_TfromTab(tab,prop) result(T)
    implicit none
    real(8), intent(in) :: prop, tab(Tmin:Tmax)
    real(8) :: T
    integer :: i, j, k

    if (prop /= prop) then
      T = prop
    elseif (prop < tab(Tmin)) then
      T = real(Tmin,8) + (prop-tab(Tmin))/(tab(Tmin+1)-tab(Tmin))
    elseif (prop >= tab(Tmax)) then
      T = real(Tmax,8) + (prop-tab(Tmax))/(tab(Tmax)-tab(Tmax-1))
    else
      i = Tmin
      j = Tmax
      do while (j-i > 1)
        k = (i+j)/2
        if (tab(k) <= prop) then
          i = k
        else
          j = k
        endif
      enddo
      T = real(i,8) + (prop-tab(i))/(tab(i+1)-tab(i))
    endif

  end function comp_TfromTab

  !> Forward difference of a property table at T.
  pure function comp_derivativeTab(tab,T) result(dprop)
    implicit none
    real(8), intent(in) :: T, tab(Tmin:Tmax)
    real(8) :: dprop
    integer :: Tprev

    Tprev = idint(T)
    dprop = tab(Tprev+1)-tab(Tprev)

  end function comp_derivativeTab

  !> Linear between the nodes, the end value outside; a NaN temperature returns NaN.
  pure function tableValue(tab, lo, hi, T) result(v)
    integer,  intent(in) :: lo, hi
    real(R8), intent(in) :: tab(lo:hi), T
    real(R8) :: v
    real(R8) :: Tc
    integer  :: i
    if (T /= T) then
      v = T
      return
    endif
    Tc = min(max(T, real(lo, R8)), real(hi, R8))
    if (Tc >= real(hi, R8)) then
      v = tab(hi)
      return
    endif
    i  = int(Tc)
    v  = tab(i) + (tab(i+1) - tab(i))*(Tc - real(i, R8))
  end function tableValue


  !> Maps the header names to column indices of vars (column 1 is Temperature, read as the mesh).
  pure subroutine classify_table_tokens(tokens, icp, irho, ih, ips, relative, code)
    character(len=*), intent(in)  :: tokens(:)
    integer,          intent(out) :: icp, irho, ih, ips, code
    logical,          intent(out) :: relative
    integer :: t

    icp = 0; irho = 0; ih = 0; ips = 0; relative = .true.
    code = TAB_NO_TEMPERATURE
    if (size(tokens) < 1) return
    if (tokens(1) /= 'Temperature' .and. tokens(1) /= 'temperature') return
    code = TAB_OK
    do t = 2, size(tokens)
      select case (trim(tokens(t)))
      case ('Cp', 'cp')
        if (icp > 0) code = TAB_DUPLICATE
        icp = t - 1
      case ('Density', 'density')
        if (irho > 0) code = TAB_DUPLICATE
        irho = t - 1
      case ('Enthalpy', 'enthalpy', 'Enthalpy_abs', 'enthalpy_abs')
        if (ih > 0) code = merge(TAB_DUPLICATE, TAB_TWO_ENTHALPY, &
                                 relative .eqv. (tokens(t) == 'Enthalpy' .or. tokens(t) == 'enthalpy'))
        ih = t - 1
        relative = (tokens(t) == 'Enthalpy' .or. tokens(t) == 'enthalpy')
      case ('Psat', 'psat', 'PSAT')
        if (ips > 0) code = TAB_DUPLICATE
        ips = t - 1
      end select
      if (code /= TAB_OK) return
    enddo
    if (icp == 0) code = TAB_NO_CP
    if (irho == 0) code = TAB_NO_DENSITY
    if (ih == 0) code = TAB_NO_ENTHALPY
  end subroutine classify_table_tokens


  !> At least two rows, none below 0 K, each within 1e-6 K of the integer node Tmin+i-1.
  pure function check_table_nodes(T) result(code)
    real(R8), intent(in) :: T(:)
    integer :: code, i, Tmin

    code = TAB_FEW_ROWS
    if (size(T) < 2) return
    code = TAB_NONFINITE
    if (.not. all(ieee_is_finite(T))) return
    code = TAB_NEGATIVE_T
    Tmin = nint(T(1))
    if (Tmin < 0) return
    code = TAB_OFF_NODE
    do i = 1, size(T)
      if (abs(T(i) - real(Tmin+i-1, R8)) > 1.e-6_R8) return
    enddo
    code = TAB_OK
  end function check_table_nodes


  !> Finite, positive density and cp, strictly increasing enthalpy consistent with cp (a
  !  varying cp: the trapezoid within 1e-3 of the step or 1e-6 of |h|, which absorbs the seams
  !  of a polynomial fit), and a relative enthalpy without offset when cp is constant.
  pure function check_table_columns(T, cp, rho, h, relative) result(code)
    real(R8), intent(in) :: T(:), cp(:), rho(:), h(:)
    logical,  intent(in) :: relative
    integer  :: code, i, n
    real(R8) :: hOff, dh

    n = size(T)
    code = TAB_NONFINITE
    if (.not. (all(ieee_is_finite(cp)) .and. all(ieee_is_finite(rho)) .and. all(ieee_is_finite(h)))) return
    code = TAB_RHO_NONPOSITIVE
    if (any(rho <= 0._R8)) return
    code = TAB_CP_NONPOSITIVE
    if (any(cp <= 0._R8)) return
    code = TAB_H_NONMONOTONE
    if (any(h(2:n) <= h(1:n-1))) return
    code = TAB_H_CP_MISMATCH
    if (all(cp == cp(1))) then
      hOff = table_h_offset(T, cp, h)
      do i = 1, n
        if (abs(h(i) - cp(1)*T(i) - hOff) > 1.e-6_R8*max(abs(h(i)), cp(1)*T(i))) return
      enddo
      code = TAB_DATUM_MISMATCH
      if (relative .and. abs(hOff) > 1.e-6_R8*cp(1)*T(1)) return
    else
      do i = 1, n-1
        dh = h(i+1) - h(i)
        if (abs(dh - 0.5_R8*(cp(i) + cp(i+1))*(T(i+1) - T(i))) > &
            max(1.e-3_R8*dh, 1.e-6_R8*max(abs(h(i)), abs(h(i+1))))) return
      enddo
    endif
    code = TAB_OK
  end function check_table_columns


  !> Enthalpy at 0 K along the first segment: h(Tmin) - cp*Tmin for a constant cp.
  pure function table_h_offset(T, cp, h) result(hOff)
    real(R8), intent(in) :: T(:), cp(:), h(:)
    real(R8) :: hOff
    if (all(cp == cp(1))) then
      hOff = h(1) - cp(1)*T(1)
    else
      hOff = h(1) - T(1)*(h(2) - h(1))
    endif
  end function table_h_offset


  !> Finite; all zero means absent; non-negative; non-decreasing; not constant; the boiling
  !  temperature inside [Tmin, Tmax-1] and psat there within a factor 2 of one atmosphere.
  pure function validate_psat_column(tab, lo, hi, Tboil, Patm) result(code)
    integer,  intent(in) :: lo, hi
    real(R8), intent(in) :: tab(lo:hi), Tboil, Patm
    integer  :: code, i
    real(R8) :: ratio
    code = PSAT_NONFINITE
    if (.not. all(ieee_is_finite(tab))) return
    code = PSAT_ABSENT
    if (all(tab == 0._R8)) return
    code = PSAT_NEGATIVE
    if (any(tab < 0._R8)) return
    code = PSAT_NONMONOTONE
    if (any(tab(lo+1:hi) < tab(lo:hi-1))) return
    code = PSAT_CONSTANT
    if (all(tab == tab(lo))) return
    code = PSAT_TBOIL_RANGE
    if (.not. (Tboil >= real(lo, R8) .and. Tboil <= real(hi-1, R8))) return
    i = int(Tboil)
    ratio = (tab(i) + (tab(i+1) - tab(i))*(Tboil - real(i, R8)))/Patm
    code = PSAT_TBOIL_MISMATCH
    if (ratio < 0.5_R8 .or. ratio > 2._R8) return
    code = PSAT_OK
  end function validate_psat_column


  pure function psat_reason(code) result(txt)
    integer, intent(in) :: code
    character(len=:), allocatable :: txt
    select case (code)
    case (PSAT_NONFINITE);      txt = 'a value that is not finite'
    case (PSAT_NEGATIVE);       txt = 'a negative pressure'
    case (PSAT_NONMONOTONE);    txt = 'a pressure that decreases with T'
    case (PSAT_CONSTANT);       txt = 'a constant pressure'
    case (PSAT_TBOIL_RANGE);    txt = 'boiling-temperature outside [Tmin, Tmax-1] of the table'
    case (PSAT_TBOIL_MISMATCH); txt = 'psat(boiling-temperature) is not within a factor 2 of one atmosphere'
    case default;               txt = 'unknown Psat error'
    end select
  end function psat_reason


  pure function table_reason(code) result(txt)
    integer, intent(in) :: code
    character(len=:), allocatable :: txt
    select case (code)
    case (TAB_NO_TEMPERATURE); txt = 'the first column is not "Temperature"'
    case (TAB_NO_CP);          txt = 'no "Cp" column'
    case (TAB_NO_DENSITY);     txt = 'no "Density" column'
    case (TAB_NO_ENTHALPY);    txt = 'no "Enthalpy" or "Enthalpy_abs" column'
    case (TAB_TWO_ENTHALPY);   txt = 'both "Enthalpy" and "Enthalpy_abs" (the enthalpy datum is ambiguous)'
    case (TAB_DUPLICATE);      txt = 'a column named twice'
    case (TAB_FEW_ROWS);       txt = 'fewer than two rows'
    case (TAB_OFF_NODE);       txt = 'rows not on consecutive integer kelvins'
    case (TAB_NEGATIVE_T);     txt = 'a temperature below 0 K'
    case (TAB_NONFINITE);      txt = 'a value that is not finite'
    case (TAB_RHO_NONPOSITIVE); txt = 'a density that is not positive'
    case (TAB_CP_NONPOSITIVE); txt = 'a cp that is not positive'
    case (TAB_H_NONMONOTONE);  txt = 'an enthalpy that does not increase with T'
    case (TAB_H_CP_MISMATCH);  txt = 'an enthalpy that disagrees with cp'
    case (TAB_DATUM_MISMATCH); txt = 'a relative "Enthalpy" column with an offset (h(Tmin) /= cp*Tmin: tag it "Enthalpy_abs")'
    case default;              txt = 'unknown table error'
    end select
  end function table_reason


  !> The quoted names of the first line that contains VARIABLES.
  subroutine read_variables_line(file, tokens, ntok, found)
    character(len=*),       intent(in)  :: file
    character(len=tok_len), intent(out) :: tokens(ntok_max)
    integer,                intent(out) :: ntok
    logical,                intent(out) :: found
    character(len=4096) :: line
    integer :: unit, ios, i, j, k

    ntok = 0; found = .false.; tokens = ''
    open(newunit=unit, file=trim(file), status='old', action='read', iostat=ios)
    if (ios /= 0) return
    do
      read(unit, '(A)', iostat=ios) line
      if (ios /= 0) exit
      if (index(line, 'VARIABLES') == 0) cycle
      found = .true.
      k = 1
      do while (ntok < ntok_max .and. k < len(line))
        i = index(line(k:), '"')
        if (i == 0) exit
        i = k + i - 1
        j = index(line(i+1:), '"')
        if (j == 0) exit
        j = i + j
        ntok = ntok + 1
        tokens(ntok) = line(i+1:j-1)
        k = j + 1
      enddo
      exit
    enddo
    close(unit)
  end subroutine read_variables_line


  !> What ORION's point reader does not check, line by line: the ZONE lines (its zone count), the
  !  lines holding I= (it gives each one a block), the rows each of them announces and each block
  !  of rows holds, every data row holding ntok numbers (ORION keeps only the last row's status),
  !  and text after the last row.
  subroutine scan_rows(file, ntok, nzone, nsize, ndata, nannounced, nrows, ntrail, badline)
    character(len=*), intent(in)  :: file
    integer,          intent(in)  :: ntok
    integer,          intent(out) :: nzone, nsize, ndata, nannounced(nzone_max), nrows(nzone_max), ntrail, badline
    character(len=4096) :: line
    real(R8) :: x, v(ntok_max)
    integer  :: unit, ios, rd, iline
    logical  :: numeric, prev_numeric

    nzone = 0; nsize = 0; ndata = 0; nannounced = -1; nrows = 0; ntrail = 0; badline = 0
    prev_numeric = .false.; iline = 0
    open(newunit=unit, file=trim(file), status='old', action='read', iostat=ios)
    if (ios /= 0) return
    do
      read(unit, '(A)', iostat=ios) line
      if (ios /= 0) exit
      iline = iline + 1
      if ((index(line, 'ZONE') > 0 .and. index(line, 'ZONETYPE') == 0) .or. index(line, 'Zone') > 0 .or. &
          index(line, 'ZONE T') > 0) nzone = nzone + 1     ! ORION's own rule
      if (index(line, 'I=') > 0) then
        nsize = nsize + 1
        if (nsize <= nzone_max) nannounced(nsize) = announced_rows(line)
      endif
      read(line, *, iostat=rd) x
      numeric = (rd == 0 .and. index(line, 'DATA') == 0)
      if (numeric) then
        if (.not. prev_numeric) ndata = ndata + 1
        if (ndata <= nzone_max) nrows(ndata) = nrows(ndata) + 1
        ntrail = 0
        if (badline == 0) then
          read(line, *, iostat=rd) v(1:ntok)
          if (rd /= 0 .or. index(line, '/') > 0) badline = iline
        endif
      else if (ndata > 0) then
        ntrail = ntrail + 1
      endif
      prev_numeric = numeric
    enddo
    close(unit)
  end subroutine scan_rows


  !> The row count of a size line as ORION takes it: I= in one of its first two comma-separated fields.
  pure function announced_rows(line) result(n)
    character(len=*), intent(in) :: line
    integer :: n, f, c, k, ios
    character(len=len(line)) :: rest, field
    n = -1
    rest = line
    do f = 1, 2
      c = index(rest, ',')
      if (c > 0) then
        field = rest(:c-1)
        rest  = rest(c+1:)
      else
        field = rest
        rest  = ''
      endif
      k = index(field, 'I=')
      if (k > 0) then
        read(field(k+2:), *, iostat=ios) n
        if (ios /= 0) n = -1
        return
      endif
    enddo
  end function announced_rows

endmodule IGLOO_Lib_Properties