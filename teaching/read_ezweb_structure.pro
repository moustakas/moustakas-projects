; $Source: /usr/users/townsend/cvsroot/www/resource/download/ez-web/read_ezweb_structure.pro,v $
;
; read_ezweb_structue()
;
; Read an EZ-Web structure file and
; return the data as an IDL structure
;
; Example:
;
;  s = read_ezweb_structure('structure_000001.txt')
;
;  PLOT, s.r, s.T, XTITLE='r (m)', YTITLE='T (K)'
;
; $Log: read_ezweb_structure.pro,v $
; Revision 1.4  2010/12/04 18:58:39  townsend
; Fixed bug in calculation of l_m and nabla_E
;
; Revision 1.3  2010/03/20 17:43:13  townsend
; Updated
;
; Revision 1.1  2009/12/04 19:40:34  townsend
; Added sources
;

; Read tabular data

PRO read_tab, filename, n_columns, data, magic=magic, double=double

; Open the file and (if required) work out the number 
; of columns

  OPENR, lun, /GET_LUN, filename

  line = ' '

  if(n_columns eq 0) then begin
     READF, lun, line
     n_columns = N_ELEMENTS(STRSPLIT(line, /EXTRACT))
     POINT_LUN, lun, 0
  endif

; Initialize the data array

  n_rows = ULONG(0)

  if(KEYWORD_SET(double)) then begin
     data = DBLARR(n_columns,64)
  endif else begin
     data = FLTARR(n_columns,64)
  endelse

; Now read the file until eof

  while not EOF(lun) do begin

; Check that there is enough space in the data
; array for the next row - if not, double
; the size of the array

     if(n_rows eq (SIZE(data))[2]) then begin
        temp = data[*,0:n_rows-1]
        if(KEYWORD_SET(double)) then begin
           data = DBLARR(n_columns,n_rows*2)
        endif else begin
           data = FLTARR(n_columns,n_rows*2)
        endelse
        data[*,0:n_rows-1] = temp
     endif

; Read in a line of data

     READF, lun, line

; Store the line in the data array

     data[*,n_rows] = STRSPLIT(line, /EXTRACT)

; Increment the number of rows

     n_rows = n_rows+1

  endwhile

  FREE_LUN, lun

; Trim the data array down to those entries which
; actually contrain data

  data = data[*,0:n_rows-1]

; Replace magic values with NaN

  if(KEYWORD_SET(magic)) then begin
     indices = WHERE(data eq magic, count)
     if(count ne 0) then data[indices] = !VALUES.F_NAN
  endif

; Finish

END

; Read EZ-Web structure file

FUNCTION read_ezweb_structure, filename

; Constants (SI)

  G_GRAVITY = 6.6742D-11
  C_LIGHT = 2.99792458D8

  M_SUN = 1.9891D30
  R_SUN = 6.96D8
  L_SUN = 3.826D26

  AMU = 1.66053886D-27

; Read the data

  read_tab, filename, 0, d, /DOUBLE

  n = N_ELEMENTS(d[0,*])

; Set up the structure

  s = CREATE_STRUCT('n_grid', n,               $
                    'm', reform(d[0,*])*M_sun,         $
                    'r', reform(d[1,*])*R_sun,         $
                    'l', reform(d[2,*])*L_sun,         $
                    'P', reform(d[3,*]),               $
                    'rho', reform(d[4,*]),             $
                    'T', reform(d[5,*]),               $
                    'u_e', reform(d[6,*]),             $
                    's', reform(d[7,*]),               $
                    'c_P', reform(d[8,*]),             $
                    'Gamma_1', reform(d[9,*]),         $
                    'nabla_ad', reform(d[10,*]),       $
                    'mu', reform(d[11,*]),             $
                    'n_e', reform(d[12,*]),            $
                    'P_e', reform(d[13,*]),            $
                    'P_rad', reform(d[14,*]),          $
                    'nabla_rad', reform(d[15,*]),      $
                    'nabla', reform(d[16,*]),          $
                    'v_c', reform(d[17,*]),            $
                    'kappa', reform(d[18,*]),          $
                    'epsilon_nuc', reform(d[19,*]),    $
                    'epsilon_pp', reform(d[20,*]),     $
                    'epsilon_CNO', reform(d[21,*]),    $
                    'epsilon_3alpha', reform(d[22,*]), $
                    'epsilon_nu_nuc', reform(d[23,*]), $
                    'epsilon_nu', reform(d[24,*]),     $
                    'epsilon_grav', reform(d[25,*]),   $
                    'X_H', reform(d[26,*]),            $
                    'X_H2', reform(d[27,*]),           $
                    'X_Hp', reform(d[28,*]),           $
                    'X_He', reform(d[29,*]),           $
                    'X_Hep', reform(d[30,*]),          $
                    'X_Hepp', reform(d[31,*]),         $
                    'X_C', reform(d[32,*]),            $
                    'X_N', reform(d[33,*]),            $
                    'X_O', reform(d[34,*]),            $
                    'psi', reform(d[35,*]),            $
                    'g', FLTARR(n),            $
                    'H_P', FLTARR(n),          $
                    'delta', FLTARR(n),        $
                    'F', FLTARR(n),            $
                    'F_R', FLTARR(n),          $
                    'F_C', FLTARR(n),          $
                    'l_m', FLTARR(n),          $
                    'nabla_E', FLTARR(n),      $
                    'n', FLTARR(n),            $
                    'P_gas', FLTARR(n),        $
                    'U', FLTARR(n),            $
                    'V', FLTARR(n))

; Set up quantities not defined in the file

  s.g = d[0,*]/d[1,*]^2*G_GRAVITY*M_SUN/R_SUN^2

  s.H_P = s.P/(s.g*s.rho)

  s.delta = s.rho*s.T*s.c_P*s.nabla_ad/s.P

  s.F = s.l/(4*!PI*s.r^2)
  s.F_C = s.F*(1. - s.nabla/s.nabla_rad)
  s.F_R = s.F - s.F_C

  s.l_m = 4.*s.rho*s.v_c^3*s.T*s.c_P/(s.g*s.delta*s.F_C)

  s.nabla_E = s.nabla - 2*s.F_C*s.H_P/(s.rho*s.v_c*s.T*s.c_P*s.l_m)

  s.n = s.rho/(s.mu*AMU)

  s.P_gas = s.P - s.P_rad

  s.U = 4.*!PI*s.rho*s.r^3/s.m

  s.V = s.rho*s.g*s.r/s.P

; Return the structure

  return, s

; Finish

END
