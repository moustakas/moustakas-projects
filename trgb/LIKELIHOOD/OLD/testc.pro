pro edge_detect, lf, cutarray, lnlike, alpha=alpha, beta=beta, gamma=gamma

        nstars = long(total(lf.hist))	   ; number of stars
        ncut = n_elements(cutarray)	   ; number of steps in the TRGB discontinuity
        lnlike = fltarr(lf.nhist,ncut)     ; 2D likelihood array

; bin the magnitudes and determine the indices of the 

        indx = fix((lf.mags-lf.minmag)*(1./lf.binsize))
        indx = indx[where((indx ge 0) and (indx lt lf.nhist),nmbins)]
        print, nmbins

; ----------------------------------------------------------------------
; C method
; ----------------------------------------------------------------------

        t1 = systime(/sec)
        doit = call_external('/deep1/ioannis/trgb/idl/LIKELIHOOD/mlike.so', $
                             'mlike_idl', long(lf.nhist), float(lf.mag), $
                             float(cutarray), float(alpha), float(beta), $
                             float(gamma), long(nstars), long(ncut), $
                             long(indx), long(nmbins), lnlike)
        print, 'C time = '+strn(systime(/sec)-t1)+'.'
        
; ----------------------------------------------------------------------
; IDL method
; ----------------------------------------------------------------------

;	marray = fltarr(lf.nhist,lf.nhist) ; shifts the origin to the TRGB
;	for k = 0L, lf.nhist-1L do marray[*,k] = lf.mag-lf.mag[k]

;        t2 = systime(/sec)
;        for j = 0L, ncut-1L do begin
;            cutfactor = cutarray[j]
;            g = form_g(marray,cutfactor,alpha=alpha,beta=beta,gamma=gamma)
;            lnlike[*,j] = - nstars * alog10(total(g,1)) + total(alog10(g[indx,*]),1)
;            print, lnlike[378,j], (alog10(total(g,1)))[378], (total(alog10(g[indx,*]),1)[378]
;        endfor
;        print, 'IDL time = '+strn(systime(/sec)-t2)+'.'

; ----------------------------------------------------------------------

        lnlike = lnlike - max(lnlike)	; normalize to zero
        print, lnlike[*,10]
stop
return
end

pro testc, objname, iter=iter, binsize=binsize, minmag=minmag, maxmag=maxmag, $
                alpha=alpha, beta=beta, gamma=gamma, hst=hst, halo=halo, core=core, $
                ccut=ccut, help=help
;+
; NAME:
;   TRGB_MLIKE
;
; PURPOSE:
;
; INPUTS:
;
; KEYWORD PARAMETERS:
;
; MODIFICATION HISTORY:
;-


	!x.ticklen=0.05
        p_thick = !p.thick & chars = !p.charsize
        charth = !p.charthick & xth = !x.thick & yth = !y.thick
        !p.thick = 1.5 & !p.charsize = 1.5 & !p.charthick = 1.5 & !x.thick=2 & !y.thick=2

	if not keyword_set(iter) then iter = 500L
        if not keyword_set(binsize) then binsize = 0.01

        if not keyword_set(alpha) then alpha = 0.33 ; slope of the RGB
	if not keyword_set(beta) then beta = 0.8    ; slope of the TRGB jump
	if not keyword_set(gamma) then gamma = 5.0  ; slope of the AGB/foreground

        if keyword_set(help) then begin
            print
            print, 'Syntax : trgb_mlike, objname, iter=iter, binsize=binsize, '
            print, 'minmag=minmag, maxmag=maxmag, hst=hst, halo=halo, core=core '
            print, 'ccut=ccut, help=help'
            print
            return
        endif

        objname = 'SextansB'
        ccut = 1
        binsize = 0.05
        
        colortable1

; read in the data

	trgb_readata, objname, datapath, data, infobase, hst=hst, halo=halo, core=core, ccut=ccut

        data[3,*] = data[3,*] - infobase.a_i	; apply the extinction correction
        nstars = n_elements(data[0,*])

        if not keyword_set(minmag) then minmag = 20. ; min(data[3,*])
        if not keyword_set(maxmag) then maxmag = 25. ; max(data[3,*])
 
; generate a luminosity function

        trgb_lfunction, data[3,*], data[4,*], lf, binsize=binsize, minmag=minmag, maxmag=maxmag

        print, 'There are '+strn(total(lf.hist),form='(I6)')+' out of '+strn(nstars,form='(I6)')+ $
          ' stars between '+strn(minmag,form='(F5.2)')+' and '+strn(maxmag,form='(F5.2)')+' mags.'
        print

; find the TRGB        

        cutarray = findgen(11) ; log_10 of the discontinuity (1-0.01 mag): [0,2] mag in 0.05 mag increments
        ncut = n_elements(cutarray)
        
        edge_detect, lf, cutarray, lnlike, alpha=alpha, beta=beta, gamma=gamma

; regenerate the best fit

        lnmax = max(lnlike,maxindx)

        x1 = maxindx mod lf.nhist
        y1 = maxindx/lf.nhist

        trgb = lf.mag[x1]
        trgbwidth = cutarray[y1]

        print & print, 'TRGB: '+strn(trgb,forma='(F6.2)')
        print, 'TRGB Width: '+strn(trgbwidth,forma='(F5.2)')

        gbest = form_g(lf.mag-trgb,trgbwidth,alpha=alpha,beta=beta,gamma=gamma)
        gnorm = nstars/total(gbest) ; normalization constant
        print, gnorm
        
; plot the results

        trgb_lfunction, data[3,*], data[4,*], lfplot, binsize=0.1, minmag=minmag, maxmag=maxmag

        window, 0, xs=300, ys=300
        plot, lfplot.mag, alog10(lfplot.hist), ps=10, xtit='I Magnitude', ytit='Log Number Density', $
          color=7, xsty=3, ysty=3, position=[0.15,0.52,0.95,0.92]
        oplot, lf.mag, alog10(gnorm*gbest), line=0, color=5, thick=1.5
        contour, lnlike, lf.mag, cutarray, levels=[-3,-2,-1], title=infobase.object, $
          ytitle='Log Step Amplitude', position=[0.15,0.12,0.95,0.52], color=3, $
          xtickname=replicate(' ',10), xsty=3, ysty=3, /noerase
        
        !x.ticklen=0
        !p.thick = p_thick & !p.charsize=chars & !p.charthick=charth
        !x.thick = xth & !y.thick=yth

stop
return
end
