; $Source: /usr/users/townsend/cvsroot/www/resource/download/ez-web/read_ezweb_summary.pro,v $
;
; read_ezweb_summary()
;
; Read an EZ-Web summary file and
; return the data as an IDL structure
;
; Example: plot a Hertzsprung-Russell diagram
;          from a summary file
;
;  s = read_ezweb_structure('structure_000001.txt')
;
;  PLOT, ALOG10(s.T), ALOG10(s.L/3.84d26), $
;     XTITLE='log T (K)', YTITLE='log L/L_sun', $
;     XRANGE=[MAX(ALOG10(s.T)),MIN(ALOG10(s.T))]
;
; $Log: read_ezweb_summary.pro,v $
; Revision 1.2  2010/03/20 17:43:13  townsend
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

; Read EZ-Web summary

FUNCTION read_ezweb_summary, filename

; Constants (SI)

  M_SUN = 1.9891D30
  R_SUN = 6.96D8
  L_SUN = 3.826D26

  YEAR = 365.25*24.*3600.

; Read the data

  read_tab, filename, 0, d, /DOUBLE

  n = N_ELEMENTS(d[0,*])

; Set up the structure

  s = CREATE_STRUCT('n_time', n,            $
                    'i', reform(ROUND(d[0,*])),     $
                    't', reform(d[1,*]),       $
;                   't', reform(d[1,*]*YEAR),       $
                    'M', reform(d[2,*]*M_SUN),      $
                    'L', reform(10^d[3,*]*L_SUN),   $
                    'R', reform(10^d[4,*]*R_SUN),   $
                    'T_s', reform(10^d[5,*]),       $
                    'T_c', reform(10^d[6,*]),       $
                    'rho_c', reform(10^d[7,*]),     $
                    'P_c', reform(d[8,*]),          $
                    'Psi_c', reform(d[9,*]),        $
                    'X_c', reform(d[10,*]),         $
                    'Y_c', reform(d[11,*]),         $
                    'X_Cc', reform(d[12,*]),        $
                    'X_Nc', reform(d[13,*]),        $
                    'X_Oc', reform(d[14,*]),        $
                    't_dyn', reform(d[15,*]*YEAR),  $
                    't_KH', reform(d[16,*]*YEAR),   $
                    't_nuc', reform(d[17,*]*YEAR),  $
                    'L_PP', reform(d[18,*]*L_SUN),  $
                    'L_CNO', reform(d[19,*]*L_SUN), $
                    'L_3a', reform(d[20,*]*L_SUN),  $
                    'L_Z', reform(d[21,*]*L_SUN),   $
                    'L_nu', reform(d[22,*]*L_SUN),  $
                    'M_He', reform(d[23,*]*M_SUN),  $
                    'M_C', reform(d[24,*]*M_SUN),   $
                    'M_O', reform(d[25,*]*M_SUN),   $
                    'R_He', reform(d[26,*]*R_SUN),  $
                    'R_C', reform(d[27,*]*R_SUN),   $
                    'R_O', reform(d[28,*]*R_SUN))

; Return the structure

  return, s

; Finish

END

