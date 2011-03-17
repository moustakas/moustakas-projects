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
          line=0, yminor=3, xminor=3, $
          xsty=3, ysty=3, ytitle='Number Density', yrange = [0,max(lf.hist)], $
          xtickname=replicate(' ',10), position=[0.12,0.32,0.504,0.92];, color=1
        if keyword_set(smf) then $
          oplot, lf.mag[goodlin], lf.hist[goodlin], color=7 else $
          oplot, lf.mag[goodlin], lf.hist[goodlin], color=7, ps=10
        xyouts, [0.52,0.52], [0.94,0.94], infobase.truename+' Luminosity Function', $
          /normal, charthick=2, charsize=2.2, align=0.5;, color=1
; ----------------------------------------------------------------------
; linear TRGB
        oplot, [trgblin,trgblin], [!y.crange[0],!y.crange[1]], line=2, thick=3;, color=2
; ----------------------------------------------------------------------
; legend
        legend, ['Stars = '+strn(nstars), 'Binsize = '+strn(lf.binsize,format='(F4.2)'), $
                 'TRGB = '+strn(trgblin,format='(F6.2)')+$
                 ' '+strn(trgblinerr,format='(F4.2)')], $
          charsize = 1.2, textcolor=[1,1,1], box=0, /clear ;, /right, /bottom
; ----------------------------------------------------------------------
; linear response
        plot, [min(lf.mag[goodlin]),max(lf.mag[goodlin])], $
          [min(resplin[goodlin]),max(resplin[goodlin])], /nodata, $
          line=0, ps=10, xsty=3, ysty=3, position=[0.12,0.12,0.504,0.32], /noerase, $
          xtit='I Magnitude', ytit='Linear Response', yminor=2, xminor=3;, color=1
        oplot, lf.mag[goodlin], resplin[goodlin];, color=5
;       oploterror, lf.mag[goodlin], resplin[goodlin], respsiglin[goodlin], errcolor=5
        oplot, [trgblin,trgblin], [!y.crange[0],!y.crange[1]], line=2, thick=3;, color=2
; ----------------------------------------------------------------------

; ----------------------------------------------------------------------
; log luminosity function
        plot, [min(lf.mag[goodlog]),max(lf.mag[goodlog])], $
          [min(alog10(lf.hist[goodlog])),max(alog10(lf.hist[goodlog]))*1.1], /nodata, $
          line=0, /noerase, yminor=3, xminor=3, $
          xsty=3, ysty=7, yrange = [min(alog10(lf.hist[goodlog])),max(alog10(lf.hist))*1.1], $
          position=[0.504,0.32,0.908,0.92], $
          ytickname = replicate(' ',10), xtickname=replicate(' ',10);, color=1
        axis, yaxis=1, yminor=3, xminor=3, ysty=7;, color=1
        if keyword_set(smf) then $
          oplot, lf.mag[goodlog], alog10(lf.hist[goodlog]), color=7 else $
          oplot, lf.mag[goodlog], alog10(lf.hist[goodlog]), color=7, ps=10
        xyouts, [0.988,0.988], [0.62,0.62], 'Log Number Density', $
          /normal, charthick=2, charsize=2.2, align=0.5, orient=90;, color=1
; ----------------------------------------------------------------------
; log TRGB
        oplot, [trgblog,trgblog], [!y.crange[0],!y.crange[1]], line=2, thick=3;, color=2
; ----------------------------------------------------------------------
; RGB line
        oplot, logaxis, logline, line=3, thick=3;, color=5
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
          [min(resplog[goodlog]),max(resplog[goodlog])], line=0, ps=10, $
          xsty=3, ysty=7, position=[0.504,0.12,0.908,0.32], /noerase, /nodata, $
          ytickname=replicate(' ',10), yminor=2, xminor=3, xtit='I Magnitude'
        axis, yaxis=1, yminor=3, xminor=3, ysty=7;, color=1
        oplot, lf.mag[goodlog], resplog[goodlog];, color=5
;       oploterror, lf.mag[goodlog], resplog[goodlog], respsiglog[goodlog], errcolor=5
        xyouts, [0.988,0.988], [0.22,0.22], 'Log Response', $
          /normal, charthick=2, charsize=2.2, align=0.5, orient=90;, color=1
        oplot, [trgblog,trgblog], [!y.crange[0],!y.crange[1]], line=2, thick=3;, color=2
; ----------------------------------------------------------------------

return
end

pro trgb_edge, objname, iter=iter, binsize=binsize, minmag=minmag, maxmag=maxmag, $
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
; COMMENTS:
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
        if not keyword_set(binsize) then binsize = 0.02

        if keyword_set(help) then begin
            print
            print, 'Syntax : trgb_edge, objname, iter=iter, binsize=binsize, '
            print, 'minmag=minmag, maxmag=maxmag, smf=smf, hst=hst, halo=halo, '
            print, 'core=core, ccut=ccut, help=help'
            print
            return
        endif

;        colortable1

; read in the data

        trgb_readata, objname, datapath, data, infobase, hst=hst, halo=halo, core=core, ccut=ccut

        data[3,*] = data[3,*] - infobase.a_i	; apply the extinction correction
        nstars = n_elements(data[0,*])

        if not keyword_set(minmag) then minmag = min(data[3,*])
        if not keyword_set(maxmag) then maxmag = max(data[3,*])

; generate a luminosity function

        trgb_lfunction, data[3,*], data[4,*], lf, binsize=binsize, minmag=minmag, maxmag=maxmag

        t1 = systime(/seconds)
        if keyword_set(smf) then smf_lfunction, lf, del=1.0
        print, 'Time: ', systime(/seconds)-t1
        
; plot the LF

        window, 0, xs=500, ys=500		; linear scaling
        plot_lfunction, infobase.truename, lf, smf=smf

        window, 1, xs=500, ys=500, ypos=30	; log scaling
        plot_lfunction, infobase.truename, lf, /log, smf=smf

; response (first derivative)

;       kernal = [-1,0,1]
;       kernal = [-1,-2,0,2,1]
;       resplin = convol(float(hist),kernal)
;       resplog = convol(alog10(hist),kernal)

        resplin = shift(float(lf.hist),-1) - shift(float(lf.hist),1)   ; linear
        resplog = shift(alog10(float(lf.hist)),-1) - shift(alog10(float(lf.hist)),1) ; log

; trap errors

        goodlin = where(finite(resplin) eq 1L,nlin)
        goodlog = where(((finite(resplog) and finite(alog10(float(lf.hist)))) eq 1L) $
                        and alog10(float(lf.hist)) gt float(0),nlog)

        if nlog lt nlin then indxarr = goodlog else indxarr = goodlin

; smooth the histograms

        linsmooth = smooth2(float(lf.hist[indxarr]),4)
        logsmooth = smooth2(alog10((float(lf.hist))[indxarr]),4)

; calculate the noise in each bin
        
        linnozero = where(linsmooth ne float(0))
        lognozero = where(logsmooth ne float(0))

        

; calculate the residuals

        linresid = abs(lf.hist[indxarr]-linsmooth)
        logresid = abs(alog10((float(lf.hist))[indxarr])-logsmooth)

        linnozero = where(linresid ne float(0))
        lognozero = where(logresid ne float(0))

; select the most conservative index array (fewest number of "good" elements)

        linweights = resplin[indxarr[linnozero]]*resplog[indxarr[linnozero]] ;/linresid[linnozero]
        logweights = resplin[indxarr[lognozero]]*resplog[indxarr[lognozero]] ;/logresid[lognozero]

;       linweights = 1./linresid[linnozero]
;       logweights = 1./logresid[lognozero]

        window, 10, xs=450, ys=450
        !p.multi = [0,1,2]
        plot, lf.mag[indxarr[linnozero]], resplin[indxarr[linnozero]]*linweights, xsty=3, ysty=3 ;, color=3
;       plot, lf.mag[indxarr], linresid, xsty=3, ysty=3;, color=3
;       plot, lf.mag[indxarr], logweights, xsty=3, ysty=3 ;, color=7
        plot, lf.mag[indxarr[lognozero]], resplog[indxarr[lognozero]]*logweights, xsty=3, ysty=3 ;, color=7
;       plot, lf.mag[indxarr], logresid, xsty=3, ysty=3;, color=7
        !p.multi = 0
        
        maxlin = max(resplin[indxarr[linnozero]]*linweights,linindx)
        maxlog = max(resplog[indxarr[lognozero]]*logweights,logindx)

;        maxlin = max(resplin[indxarr[linnozero]],linindx)
;        maxlog = max(resplog[indxarr[lognozero]],logindx)

        trgblin = (lf.mag[indxarr[linnozero]])[linindx]
        trgblog = (lf.mag[indxarr[lognozero]])[logindx]

; ------------------------------------------------------------
; linear
; ------------------------------------------------------------

        lf2 = create_struct(lf, 'response_lin', resplin, 'response_log', resplog)
        linrange = where(lf2.hist)

        boot_response, iter, lf2, linrange, sig, trgbmaglin, smf=smf ; bootstrap resample
        trgblinerr = stdev(trgbmaglin)

; print the results to the screen

        trgb_printoscreen, trgblin, trgblinerr, infobase, distlin, vpeclin

; ------------------------------------------------------------
; log
; ------------------------------------------------------------
        
        logrange = where(lf2.hist)

        boot_response, iter, lf2, logrange, sig, trgbmaglog, /log, smf=smf
        trgblogerr = stdev(trgbmaglog)

        trgb_printoscreen, trgblog, trgblogerr, infobase, distlog, vpeclog, /log

; diagnostic plot

;        window, 1, xs=450, ys=450
;        plothist, trgbmaglin, bin=0.02
;        plothist, trgbmaglog, bin=0.02, /overplot


stop
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

        okay = 'N'
        print
        read, okay, prompt='Generate a postscript file (Y/[N])? '

        if strupcase(okay) eq 'Y' then begin
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
