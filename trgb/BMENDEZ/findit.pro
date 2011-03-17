;+
; NAME:
;   findit
;
;
; PURPOSE:
;   to run likelihood analysis on a starlist, looking for TRGB in I
;   band.  Data based either on keck run or on HST study.
;
; CALLING SEQUENCE:
;    findit, objname, hst=hst, halo=halo, minmag=minmag,maxmag=maxmag,colorcut=colorcut
;
;
; INPUTS:
;   objname -- name of galaxy
;
;
; KEYWORD PARAMETERS:
;   hst -- select if hst object, otherwise keck
;   halo -- select if desire to read halo of object only
;   minmag -- minimum magnitude (I=20 default)
;   maxmag -- maximum magnitude (I=25 default)
;   colorcut -- select if desirous of color-cut
;
; MODIFICATION HISTORY:
;  jm, md  23jun00
;-
pro findit, objname, hst=hst, halo=halo, core = core, ccut = ccut, $
              minmag=minmag,maxmag=maxmag, binsize = binsize

  IF NOT KEYWORD_SET(binsize) THEN binsize = 0.1

; read in starlist 
TRGB_READATA, objname, datapath, data, info, $
            halo=halo, core = core, hst=hst, ccut = ccut

; test of likelihood analysis
imag = data[3,*]
print, N_ELEMENTS(imag), ' stars input for galaxy'

;if (N_ELEMENTS(data[*,0]) eq 8 and KEYWORD_SET(colorcut)) then begin
;     color = data[7,*] - info.e_v_i ; color excess correction
;     colorstart = 0.5
;     colorend = 2.0
;     rgbstars = where((color gt colorstart) and (color lt colorend),rgbcount)
;     imag = imag[rgbstars]   ;select subset of stars
;     print, n_elements(imag), ' stars satisfy colorcut'
;endif


if not Keyword_set(minmag) then minmag=20.
if not Keyword_set(maxmag) then maxmag=25.
 

EDGE_DETECT, objname, imag, lnlike, error, max=maxmag, min=minmag, $
             hst=hst, halo=halo, core=core, ccut = ccut

mag=findgen(50)/10. +minmag
;step=findgen(41)/20.
step = findgen(20)/20. + 0.1
beta = findgen(6)/10. +.5

dummy=max(lnlike,i)
; N.B. lnlike has dimension (50,20,6)
beta=i/(50*20)
betaf=beta/10. +.5
cut=fix((i-beta*(50*20))/50.)
;cutf=cut/20.
cutf = cut/20. + 0.1
;magcut=(i-cut*50)/10. +minmag ; find parameters of best fit
magcut=(i mod 50)/10. + minmag  
print, 'best fit paramters: m_cut = ',magcut, ' cutf= ',cutf, ' beta',betaf

gg = FORM_GGS(error, magcut-minmag, cutf, betaf)  ;regenerate best fit
mmag=findgen(500)/100. +minmag
hist=histogram(imag,bin=.1,min=minmag,max=maxmag)

; Plots
COLORTABLE1
WINDOW, 0, xs = 500, ys = 700

PLOT,mag,(hist > 1),psym=10, /ylog, $
  title = info.truename, ytitle = 'log number', $
  position = [0.1,0.4,0.95,0.95]
OPLOT,mmag,gg*hist(25)/gg(250),lines=3, color = 7, thick = 2


llnlike=lnlike[*,*,beta]
CONTOUR, llnlike,mag,step, levels=[-3,-2,-1], $
  xtitle='magnitude', ytitle='log cut amplitude', $
  position = [0.1,0.1,0.95,0.4], /NOERASE


end

