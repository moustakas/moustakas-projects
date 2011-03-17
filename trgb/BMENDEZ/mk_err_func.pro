pro MK_ERR_FUNC, objname, a, m, hst=hst, halo=halo, core=core, $
                 cutfactor = cutfactor, cutvalue = cutvalue, $
                 order = order, bin = bin, dmag = dmag, save = save

        if not keyword_set(order) then order = 2.0
        if not keyword_set(bin) then bin = 100.
        if not keyword_set(cutfactor) then cutfactor = 0.005
        if not keyword_set(dmag) then dmag = 0.05
        
; read in the data ----------------------------------------------------------
trgb_readata, objname, datapath, data, infobase, hst=hst, halo=halo, core=core, ccut=ccut

imags = TRANSPOSE(data[3,*])
merr = TRANSPOSE(data[4,*])

; Convert to Data Numbers -----------------------------------------------------
fmag = 10.^((30-imags)/2.5)
fsig = (merr/2.5)*ALOG(10.)*fmag

; Make binned luminosity function ---------------------------------------------

fhist = HISTOGRAM(fmag,BINSIZE=bin,MIN=MIN(fmag))
b = N_ELEMENTS(fhist)
f = 100*FINDGEN(b) + MIN(fmag)
m = 30 - 2.5*ALOG10(f)
fbin = (dmag/2.5)*ALOG(10.)*f
;plot, f, fbin
;STOP
;print, '# of Bins=', b

; Form rough error function ---------------------------------------------------
err = FLTARR(b)
std = FLTARR(b)
IF KEYWORD_SET(cutvalue) THEN BEGIN
   GETELEMENT_VECTOR, f, cutvalue, cutoff
   FOR k = 0L, cutoff-1L DO BEGIN
     nearerr = WHERE((fmag GE (f[k])-fbin[k]) AND $
                     (fmag LE (f[k])+fbin[k]), errcount)
     IF (errcount GT 1L) THEN err[k] = MEDIAN(fsig[nearerr]) ELSE err[k] = 0.
     IF (errcount GT 1L) THEN std[k] = STDEV(fsig[nearerr]) ELSE std[k] = 1.e10
   ENDFOR
ENDIF ELSE BEGIN
   FOR j = 0L, b-1L DO BEGIN
     nearerr = WHERE((fmag GE (f[j])-fbin[j]) AND $
                     (fmag LE (f[j])+fbin[j]), errcount)
     IF (errcount GT 1L) THEN err[j] = MEDIAN(fsig[nearerr]) ELSE err[j] = 0.
     IF (errcount GT 1L) THEN std[j] = STDEV(fsig[nearerr]) ELSE std[j] = 1.e10
   ENDFOR
ENDELSE

; Fit a smooth function to the magnitude errors -------------------------------
IF NOT KEYWORD_SET(cutvalue) THEN $ 
       cutoff = MIN(WHERE(fhist LT cutfactor*MAX(fhist)))
a = POLY_FIT(f[0:cutoff], err[0:cutoff]^2., order, yfit, /double)
y = FLTARR(b)
;print, a
yi = FLTARR(b, N_ELEMENTS(a))
FOR i = 0L, b-1L DO BEGIN
   FOR j = 0, N_ELEMENTS(a)-1L DO BEGIN
      yi[i, j] = a[j]*((f[i])^FLOAT(j))
   ENDFOR
   y[i] = TOTAL(yi[i, *])
ENDFOR

; convert into function in terms of magnitude and mag errors ------------------

sigm = FLTARR(b)
si = FLTARR(b, N_ELEMENTS(a))
FOR i = 0L, b-1L DO BEGIN
   FOR j = 0, N_ELEMENTS(a)-1L DO BEGIN
      si[i, j] = a[j]*((f[i])^FLOAT(j))
   ENDFOR
   sigm[i] = (2.5/ALOG(10.))*(SQRT(TOTAL(yi[i, *]))/f[i])
ENDFOR

; Plots -----------------------------------------------------------------------

COLORTABLE1

WINDOW, 0
PLOT, imags, merr, psym = 3, XTITLE = 'Magnitude', YTITLE = 'Magnitude Error'
OPLOT, m, sigm, color = 3, thick = 4
OPLOT, m, (2.5/ALOG(10.))*(err/f), color = 6, thick = 2

WINDOW, 5
PLOT, f, fhist, psym = 10, xrange = [MIN(f),f[cutoff]] 
OPLOT, [f[cutoff],f[cutoff]], [!y.crange[0],!y.crange[1]],$
       line = 2, thick = 2, color = 2
legend, ['0.5% of max at'+STRN(f[cutoff])], box = 0, /clear

WINDOW, 10
PLOT, fmag, fsig^2, psym = 3, $
  XTITLE = 'Flux', YTITLE = 'Flux Error', $
  xrange = [MIN(f), MAX(f)], $
  /ylog, /xlog
OPLOT, f, y, color = 3, thick = 4
OPLOT, f, err^2, color = 6, thick = 2

; Write critical info to an IDL save set --------------------------------------

IF KEYWORD_SET(save) THEN BEGIN

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

            print & print, 'Writing '+filename+'.' & print
            SAVE, bin, order, a, filename = filename

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

            print & print, 'Writing '+filename+'.' & print
            SAVE, bin, order, a, filename = filename

        endelse

        popd

ENDIF

  return
end
