pro findit, objname, hst=hst, halo=halo, core=core, minmag=minmag,maxmag=maxmag, $
     colorcut=colorcut, ccut=ccut
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


        trgb_readata, objname, datapath, data, infobase, hst=hst, halo=halo, core=core, ccut=ccut

; read in starlist 
; 
; test of likelihood analysis
;
 
imag = data[3,*] - infobase.a_i ; apply extinction correction
print, n_elements(imag), ' stars input for galaxy'
if (n_elements(data[*,0]) eq 8 and keyword_set(colorcut)) then begin
     color = data[7,*] - infobase.e_v_i ; color excess correction
     colorstart = 0.5
     colorend = 2.0
     rgbstars = where((color gt colorstart) and (color lt colorend),rgbcount)
     imag = imag[rgbstars]   ;select subset of stars
     print, n_elements(imag), ' stars satisfy colorcut'
endif


;minmag=20.
;maxmag=25.
if not Keyword_set(minmag) then minmag=20.
if not Keyword_set(maxmag) then maxmag=25.
 

edge_detect, imag, lnlike, max=maxmag, min=minmag


 
!p.multi=[0,1,2]
mag=findgen(50)/10. +minmag
step=findgen(41)/20.
contour, lnlike,mag,step, levels=[-3,-2,-1],title=objname, $
  xtitle='cut magnitude', ytitle='log cut amplitude'
dummy=max(lnlike,i)

cut=fix(i/50.)
ccut=cut/20.
magcut=(i-cut*50)/10. +minmag ; find parameters of best fit
print, 'best fit paramters: m_cut = ',magcut, ' ccut= ',ccut

gg=form_gg(magcut-minmag,ccut)  ;regenerate best fit
mmag=findgen(500)/100. +minmag
hist=histogram(imag,bin=.1,min=minmag,max=maxmag)

plot,mag,alog10(hist),psym=10,  xtitle='magnitude',ytitle='log number'
oplot,mmag,alog10(gg)+alog10(hist(25))-alog10(gg(250)),lines=1

!p.multi=0  ;reset plotter

;save,imag,lnlike,minmag, maxmag, mag,step,gg,mmag, file= filename + '.idlsv'

;plot, mag, lnlike, yr=[-3, 0]

stop
end

