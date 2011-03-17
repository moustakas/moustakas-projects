pro rgbline, trgb, lf, logaxis, logline, logcoef
; jm00aug7ucb
; overplot a line onto the logarithmic luminosity function which
; indicates the slope of the red-giant branch

        getelement_vector, lf.mag, trgb, x1
        getelement_vector, alog10(lf.hist), max(alog10(lf.hist)), x2

        logcoef = linfit(lf.mag[x1:x2],alog10(lf.hist[x1:x2]))
        logaxis = findgen((lf.mag[x2]-trgb)*100.)/100.+trgb

        logline = poly(logaxis,logcoef)

return
end

pro makeplot, lf, resplin, respsiglin, resplog, respsiglog, $
              goodlin, goodlog, trgblin, trgblinerr, trgblog, trgblogerr, $
              infobase, nstars, smf=smf
; jm00aug8ucb
; generate a luminosity function


; create a line indicating the slope of the RGB

	rgbline, trgblog, lf, logaxis, logline, logcoef

; ----------------------------------------------------------------------
; linear luminosity function
        plot, [min(lf.mag[goodlin]),max(lf.mag[goodlin])], [0,max(lf.hist[goodlin])], /nodata, $
          line=0, color=1, yminor=3, xminor=3, $
          xsty=3, ysty=3, ytitle='Number Density', yrange = [0,max(lf.hist)], $
          xtickname=replicate(' ',10), position=[0.12,0.32,0.504,0.92]
        if keyword_set(smf) then $
          oplot, lf.mag[goodlin], lf.hist[goodlin], color=7 else $
          oplot, lf.mag[goodlin], lf.hist[goodlin], color=7, ps=10
        xyouts, [0.52,0.52], [0.94,0.94], infobase.truename+' Luminosity Function', $
          /normal, charthick=2, charsize=2.2, color=1, align=0.5
; ----------------------------------------------------------------------
; linear TRGB
        oplot, [trgblin,trgblin], [!y.crange[0],!y.crange[1]], line=2, thick=3, color=2
; ----------------------------------------------------------------------
; legend
        legend, ['Stars = '+strn(nstars), 'Binsize = '+strn(lf.binsize,format='(F4.2)'), $
                 'TRGB = '+strn(trgblin,format='(F6.2)')+$
                 ' '+strn(trgblinerr,format='(F4.2)')], $
          charsize = 1.2, textcolor=[1,1,1], box=0, /clear ;, /right, /bottom
; ----------------------------------------------------------------------
; linear response
        plot, [min(lf.mag[goodlin]),max(lf.mag[goodlin])], $
          [min(resplin[goodlin]),max(resplin[goodlin])], /nodata, color=1, $
          line=0, ps=10, xsty=3, ysty=3, position=[0.12,0.12,0.504,0.32], /noerase, $
          xtit='I Magnitude', ytit='Linear Response', yminor=2, xminor=3
        oplot, lf.mag[goodlin], resplin[goodlin], color=5
;       oploterror, lf.mag[goodlin], resplin[goodlin], respsiglin[goodlin], errcolor=5
        oplot, [trgblin,trgblin], [!y.crange[0],!y.crange[1]], line=2, thick=3, color=2
; ----------------------------------------------------------------------

; ----------------------------------------------------------------------
; log luminosity function
        plot, [min(lf.mag[goodlog]),max(lf.mag[goodlog])], $
          [min(alog10(lf.hist[goodlog])),max(alog10(lf.hist[goodlog]))*1.1], /nodata, $
          line=0, color=1, /noerase, yminor=3, xminor=3, $
          xsty=3, ysty=7, yrange = [min(alog10(lf.hist[goodlog])),max(alog10(lf.hist))*1.1], $
          position=[0.504,0.32,0.908,0.92], $
          ytickname = replicate(' ',10), xtickname=replicate(' ',10)
        axis, yaxis=1, color=1, yminor=3, xminor=3, ysty=7
        if keyword_set(smf) then $
          oplot, lf.mag[goodlog], alog10(lf.hist[goodlog]), color=7 else $
          oplot, lf.mag[goodlog], alog10(lf.hist[goodlog]), color=7, ps=10
        xyouts, [0.988,0.988], [0.62,0.62], 'Log Number Density', $
          /normal, charthick=2, charsize=2.2, color=1, align=0.5, orient=90
; ----------------------------------------------------------------------
; log TRGB
        oplot, [trgblog,trgblog], [!y.crange[0],!y.crange[1]], line=2, thick=3, color=2
; ----------------------------------------------------------------------
; RGB line
        oplot, logaxis, logline, line=3, color=5, thick=3
; ----------------------------------------------------------------------
; legend
        legend, ['Stars = '+strn(nstars), 'Binsize = '+strn(lf.binsize,format='(F4.2)'), $
                 'RGB Slope = '+strn(logcoef[1],format='(F5.2)'), $
        'TRGB = '+strn(trgblog,format='(F6.2)')+$
          ' '+strn(trgblogerr,format='(F4.2)')], $
          charsize = 1.2, textcolor=[1,1,1,1], box=0, /clear ;, /right, /bottom
; ----------------------------------------------------------------------
; log response
        plot, [min(lf.mag[goodlog]),max(lf.mag[goodlog])], $
          [min(resplog[goodlog]),max(resplog[goodlog])], line=0, color=1, ps=10, $
          xsty=3, ysty=7, position=[0.504,0.12,0.908,0.32], /noerase, /nodata, $
          ytickname=replicate(' ',10), yminor=2, xminor=3, xtit='I Magnitude'
        axis, yaxis=1, color=1, yminor=3, xminor=3, ysty=7
        oplot, lf.mag[goodlog], resplog[goodlog], color=5
;       oploterror, lf.mag[goodlog], resplog[goodlog], respsiglog[goodlog], errcolor=5
        xyouts, [0.988,0.988], [0.22,0.22], 'Log Response', $
          /normal, charthick=2, charsize=2.2, color=1, align=0.5, orient=90
        oplot, [trgblog,trgblog], [!y.crange[0],!y.crange[1]], line=2, thick=3, color=2
; ----------------------------------------------------------------------

return
end

pro trgb_edge_old, objname, iter=iter, binsize=binsize, minmag=minmag, maxmag=maxmag, $
                   smf=smf, hst=hst, halo=halo, core=core, ccut=ccut, help=help
;+
; NAME:
;	TRGB_EDGE
;
; PURPOSE:
;	Detect the magnitude of the tip of the red-giant branch (or
;	maximum discontinuity) from a stellar luminosity function. 
;
; INPUTS:
;	objname : galaxy name
;
; OPTIONAL INPUTS:
;	iter	: number of iterations when bootstrapping
;
; KEYWORD PARAMETERS:
;	smf	   : apply the Sakai, Madore & Freedman edge-detection
;		     method 
;	hst  	   : keyword specifying an HST object
;	halo 	   : use the objname_starlist_halo.dat photometry file
;	             created in TRGB_REGIONS
;	ccut	   : use the objname_starlist_ccut.dat photometry file
;		     created in TRGB_COLORCUT
;	help	   : print out all the syntax of the program (all the
;		     keywords) 
;
; OUTPUTS:
;
; COMMON BLOCKS:
;
; EXAMPLE:
;	trgb_edge, 'SextansB', /halo
;		OR
;	trgb_edge, 'ugc07577', /hst, /postscript
;
; PROCEDURES USED:
;	COLORTABLE1, TRGB_READATA, TRGB_LFUNCTION, PLOT_LFUNCTION,
;	PLOT_RESPONSE, BOOT_RESPONSE, TRGB_PRINTOSCREEN, SREAD(),
;	SWRITE() 
;
; MODIFICATION HISTORY:
;	John Moustakas, 2000 June 20, UCB
;	cleaned up, jm00aug4ucb
;	
;-

	!x.ticklen=0.05
        p_thick = !p.thick & chars = !p.charsize
        charth = !p.charthick & xth = !x.thick & yth = !y.thick
        !p.thick = 2 & !p.charsize=1.8 & !p.charthick=2 & !x.thick=2 & !y.thick=2

	if not keyword_set(iter) then iter = 500L
        if not keyword_set(binsize) then binsize = 0.1
        if not keyword_set(minmag) then minmag = 18.0
        if not keyword_set(maxmag) then maxmag = 27.0

        if keyword_set(help) then begin
            print
            print, 'Syntax : trgb_edge, objname, iter=iter, binsize=binsize, '
            print, 'minmag=minmag, maxmag=maxmag, smf=smf, hst=hst, halo=halo, '
            print, 'core=core, ccut=ccut, help=help'
            print
            return
        endif

	colortable1

; read in the data

        trgb_readata, objname, datapath, data, infobase, hst=hst, halo=halo, core=core, ccut=ccut

        data[3,*] = data[3,*] - infobase.a_i	; apply the extinction correction
        nstars = n_elements(data[0,*])

; generate a luminosity function

        trgb_lfunction, data[3,*], data[4,*], lf, binsize=binsize, minmag=minmag, maxmag=maxmag

; ------------------------------------------------------------------------------------
; literature method: generate a smoothed LF as defined in Sakai, Madore, Freedman 1996
; ------------------------------------------------------------------------------------

        if keyword_set(smf) then begin

            del = 1.0           ; smooth factor
            smf = del
            for j = 0L, lf.nhist-1L do begin
                near = where(abs(lf.mags-lf.mag[j]) lt 5.*lf.merr/del,count) ; 5-sigma
                if count gt 0L then $
                  lf.hist[j] = total(1./(sqrt(2.*!pi)*del*lf.merr[near])* $
                                     exp(-0.5*((lf.mags[near]-lf.mag[j])/lf.merr[near]*del)^2))
            endfor
            
        endif

; ----------------------------------------------------------------------

; plot the LF

        window, 0, xs=500, ys=500	; linear scaling
        plot_lfunction, infobase.truename, lf, smf=smf

        window, 1, xs=500, ys=500	; log scaling
        plot_lfunction, infobase.truename, lf, /log, smf=smf

; ----------------------------------------------------------------------
; determine the response or the first derivative of the LF
; ----------------------------------------------------------------------

;       kernal = [-1,0,1]
;       kernal = [-1,-2,0,2,1]
;       resplin = convol(float(hist),kernal)
;       resplog = convol(alog10(hist),kernal)

        resplin = shift(float(lf.hist),-1) - shift(float(lf.hist),1)   ; linear
        resplog = shift(alog10(lf.hist),-1) - shift(alog10(lf.hist),1) ; log

; quantify the variation in the response by bootstrap resampling

        iter1 = long(0.1*iter)	> 10 ; number of iterations in the first pass

        range = where(lf.mag)
        lftemp = create_struct(lf, 'response_lin', resplin, 'response_log', resplog)

        boot_response, iter1, lftemp, range, respsiglin, smf=smf       ; linear
        boot_response, iter1, lftemp, range, respsiglog, /log, smf=smf ; log

        window, 2, xs=500, ys=500 	; plot it
        plot_response, infobase.truename, resplin, respsiglin, $
          resplog, respsiglog, lf.mag

; trap errors

        goodlin = where(((finite(resplin) and finite(respsiglin)) eq 1L) and (respsiglin ne float(0)))
        goodlog = where(((finite(resplog) and finite(respsiglog)) eq 1L) and (respsiglog ne float(0)))

; find the maximum response

        maxlin = max(resplin[goodlin]/respsiglin[goodlin],linindx)
        linmag = (lf.mag[goodlin])[linindx]

        maxlog = max(resplog[goodlog]/respsiglog[goodlog],logindx)
        logmag = (lf.mag[goodlog])[logindx]

; ----------------------------------------------------------------------------
; refine the TRGB detection by creating a nearly continuous LF (small
; binsize).  for the SMF method, use the same binsize as before
; ----------------------------------------------------------------------------

        if keyword_set(smf) then nbinsize = binsize else nbinsize = 0.02
        trgb_lfunction, data[3,*], data[4,*], lf2, binsize=nbinsize, minmag=18., maxmag=27.

        if keyword_set(smf) then begin

            for j = 0L, lf2.nhist-1L do begin
                near = where(abs(lf2.mags-lf2.mag[j]) lt 5.*lf2.merr/del,count) ; 5-sigma
                if count gt 0L then $
                  lf2.hist[j] = total(1./(sqrt(2.*!pi)*del*lf2.merr[near])* $
                                     exp(-0.5*((lf2.mags[near]-lf2.mag[j])/lf2.merr[near]*del)^2))
            endfor
            
        endif

; search within +/- magrange

        magrange = 0.25	; magnitudes (can be chosen)
 
; ------------------------------------------------------------
; linear
; ------------------------------------------------------------

        linrange = where((lf2.mag gt linmag-magrange) and (lf2.mag lt linmag+magrange)) ; linear

        temp = lf2.mag[linrange]
        print, 'Isolating the linear response between '+strn(min(temp),form='(F5.2)')+' and '+$
          strn(max(temp),form='(F5.2)')+' mag:' & print

        lf2temp = create_struct(lf2, 'response_lin', resplin, 'response_log', resplog)
        boot_response, iter, lf2temp, linrange, sig, trgbmaglin, smf=smf ; bootstrap resample

        trgblin = mean(trgbmaglin)
        trgblinerr = stdev(trgbmaglin)

; print the results to the screen

        trgb_printoscreen, trgblin, trgblinerr, median(trgbmaglin), infobase, distlin, vpeclin

; ------------------------------------------------------------
; log
; ------------------------------------------------------------

        logrange = where((lf2.mag gt logmag-magrange) and (lf2.mag lt logmag+magrange)) ; log

        temp = lf2.mag[logrange]
        print, 'Isolating the log response between '+strn(min(temp),form='(F5.2)')+' and '+$
          strn(max(temp),form='(F5.2)')+' mag:' & print

        boot_response, iter, lf2temp, logrange, sig, trgbmaglog, /log, smf=smf

        trgblog = mean(trgbmaglog)
        trgblogerr = stdev(trgbmaglog)

        trgb_printoscreen, trgblog, trgblogerr, median(trgbmaglog), infobase, distlog, vpeclog, /log

; diagnostic plot

;        window, 10, xs=450, ys=450
;        plothist, trgbmaglin, bin=0.02, col=5, $
;          title='Distribution of Errors', xtit='I Magnitude', ytit='Number Density'
;        plothist, trgbmaglog, bin=0.02, col=7, /overplot

; write the results to the data structure

	print & print, 'Updating '+objname+"'s data structure . . . "

        if keyword_set(hst) then $
          infofile = '/deep1/ioannis/trgb/hst_object.dat' else $
          infofile = '/deep1/ioannis/trgb/keck_object.dat'

        info = sread(infofile)
        indx = where(strupcase(info.object) eq strupcase(objname))
        info[indx].distance = distlog
        info[indx].pec_vel = vpeclog
        swrite, info, infofile

        window, 3, xs=800, ys=800, ypos=30
        makeplot, lf, resplin, respsiglin, resplog, respsiglog, $
          goodlin, goodlog, trgblin, trgblinerr, trgblog, trgblogerr, $
          infobase, nstars, smf=smf

        okay = ''
        print & print, 'Type p to create a postscript, '
        read, okay, prompt='file or press ENTER to exit: '

        if strupcase(okay) eq 'P' then begin
            ps_open, datapath+'/'+objname+'_trgb', /ps_fonts
            device, /times
            makeplot, lf, resplin, respsiglin, resplog, respsiglog, $
              goodlin, goodlog, trgblin, trgblinerr, trgblog, trgblogerr, $
              infobase, nstars
            ps_close
        endif

        while !d.window ne -1L do wdelete, !d.window ; delete all windows
        
        !x.ticklen=0
        !p.thick = p_thick & !p.charsize=chars & !p.charthick=charth
        !x.thick = xth & !y.thick=yth

return
end
