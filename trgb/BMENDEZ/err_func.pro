FUNCTION ERR_FUNC, objname, mag, plot = plot, $ 
                   hst=hst, halo=halo, core=core


; Convert to Data Numbers -----------------------------------------------------
flux = 10.^((30-mag)/2.5)

; Read in polynomial fit values -----------------------------------------------

	paths = strarr(3)
        paths[0] = '/deepscr1/ioannis/trgb/'	; HST data
        paths[1] = '/deepscr1/ioannis/archive/'	; HST archival data
        paths[2] = '/deep3/marc/trgb/data/'	; Keck data

        if keyword_set(hst) then begin ; HST data

            pushd, paths[0]+objname
            datapath = paths[0]+objname

            if keyword_set(halo) then filename = datapath+'/'+objname+'_halo_error_fit.dat' else $
            if keyword_set(core) then filename = datapath+'/'+objname+'_core_error_fit.dat' else $
              filename = datapath+'/'+objname+'_error_fit.dat'

            print & print, 'Reading '+filename+'.' & print
            RESTORE, filename

        endif else begin ; KECK data

            objdata = ['holmbergii','holmbergix','i342','ngc2366', 'ngc1560', $
                       'ngc2903','ngc2976','ngc3109','sextansb']
            check = where(strupcase(objdata) eq strupcase(objname),count)
            if count gt 0L then date = '22dec97/' else date = '23dec97/'

            pushd, paths[2]+date+strlowcase(objname)
            datapath = paths[2]+date+strlowcase(objname)

            if keyword_set(halo) then filename = datapath+'/'+objname+'_halo_error_fit.dat' else $
            if keyword_set(core) then filename = datapath+'/'+objname+'_core_error_fit.dat' else $
              filename = datapath+'/'+objname+'_error_fit.dat'

            print & print, 'Reading '+filename+'.' & print
            RESTORE, filename

        endelse

        popd


; Form function in terms of magnitude and mag errors --------------------------

sig = FLTARR(N_ELEMENTS(flux))
si = FLTARR(N_ELEMENTS(flux), N_ELEMENTS(a))
FOR i = 0, N_ELEMENTS(flux)-1L DO BEGIN
   FOR j = 0, N_ELEMENTS(a)-1L DO BEGIN
      si[i, j] = a[j]*((flux[i])^FLOAT(j))
   ENDFOR
   sig[i] = (2.5/ALOG(10.))*(SQRT(TOTAL(si[i, *]))/flux[i])
ENDFOR

; Optional Plot ---------------------------------------------------------------

IF KEYWORD_SET(plot) THEN BEGIN
   COLORTABLE1
   WINDOW, 6
   PLOT, mag, sig
ENDIF

  return, sig
end

