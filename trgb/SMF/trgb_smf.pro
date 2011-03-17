pro rgbline, trgb, lf, logaxis, logline, logcoef
; jm00aug7ucb
; overplot a line onto the logarithmic luminosity function which
; indicates the slope of the red-giant branch

        getelement_vector, lf.mag, trgb, x1
        getelement_vector, alog10(lf.hist), max(alog10(lf.hist)), x2

        logcoef = linfit(lf.mag[x1:x2],alog10(lf.hist[x1:x2]))

;       logaxis = findgen((lf.mag[x2]-trgb)*100.)/100.+trgb
        logaxis = findgen((lf.maxmag-lf.minmag)*100.)/100.+lf.minmag
        
        logline = poly(logaxis,logcoef)

return
end

pro makeplot, lf, linsmooth, logsmooth, resplin, resplog, $
              noise, trgblin, trgblinerr, trgblog, trgblogerr, infobase, nstars, smf=smf
; jm00aug8ucb

; create a line indicating the slope of the RGB

	rgbline, trgblog, lf, logaxis, logline, logcoef

; ----------------------------------------------------------------------
; linear luminosity function
        plot, [min(lf.mag),max(lf.mag)], [0,max(lf.hist)], /nodata, $
          line=0, yminor=3, xminor=3, $
          xsty=3, ysty=3, ytitle='N(m)', yrange = [0,max(lf.hist)], $
          xtickname=replicate(' ',10), position=[0.12,0.32,0.504,0.92], color=1
        if keyword_set(smf) then $
          oplot, lf.mag, lf.hist, color=7 else $
          oplot, lf.mag, lf.hist, color=7, ps=10
        oplot, lf.mag, linsmooth, color=3
        xyouts, [0.52,0.52], [0.94,0.94], infobase.truename+' Luminosity Function', $
          /normal, charthick=2, charsize=2.2, align=0.5, color=1
; ----------------------------------------------------------------------
; linear TRGB
        oplot, [trgblin,trgblin], [!y.crange[0],!y.crange[1]], line=2, thick=3, color=2
; ----------------------------------------------------------------------
; RGB line
        oplot, logaxis, 10.^logline, line=3, thick=3, color=5
; ----------------------------------------------------------------------
; legend
        legend, ['Stars '+strn(nstars), 'Binsize '+strn(lf.binsize,format='(F4.2)'), $
                 'TRGB '+strn(trgblin,format='(F6.2)')+$
                 ' '+strn(trgblinerr,format='(F4.2)')], $
          charsize = 1.2, textcolor=[1,1,1], box=0, /clear ;, /right, /bottom
; ----------------------------------------------------------------------
; linear response
        plot, [min(lf.mag),max(lf.mag)], $
          [min(resplin/noise),max(resplin/noise)], /nodata, $
          line=0, ps=10, xsty=3, ysty=3, position=[0.12,0.12,0.504,0.32], /noerase, $
          xtit='I', ytit='Linear Response', yminor=2, xminor=3, color=1
        oplot, lf.mag, resplin/noise, color=5
;       oploterror, lf.mag, resplin, respsiglin, errcolor=5
        oplot, [trgblin,trgblin], [!y.crange[0],!y.crange[1]], line=2, thick=3, color=2
; ----------------------------------------------------------------------

; ----------------------------------------------------------------------
; log luminosity function
        plot, [min(lf.mag),max(lf.mag)], $
          [min(alog10(lf.hist)),max(alog10(lf.hist))*1.1], /nodata, $
          line=0, /noerase, yminor=3, xminor=3, $
          xsty=3, ysty=7, yrange = [min(alog10(lf.hist)),max(alog10(lf.hist))*1.1], $
          position=[0.504,0.32,0.908,0.92], $
          ytickname = replicate(' ',10), xtickname=replicate(' ',10), color=1
        axis, yaxis=1, yminor=3, xminor=3, ysty=7, color=1
        if keyword_set(smf) then $
          oplot, lf.mag, alog10(lf.hist), color=7 else $
          oplot, lf.mag, alog10(lf.hist), color=7, ps=10
        oplot, lf.mag, logsmooth, color=3
        xyouts, [0.988,0.988], [0.62,0.62], 'log N(m)', $
          /normal, charthick=2, charsize=2.2, align=0.5, orient=90, color=1
; ----------------------------------------------------------------------
; log TRGB
        oplot, [trgblog,trgblog], [!y.crange[0],!y.crange[1]], line=2, thick=3, color=2
; ----------------------------------------------------------------------
; RGB line
        oplot, logaxis, logline, line=3, thick=3, color=5
; ----------------------------------------------------------------------
; legend
        legend, ['Stars '+strn(nstars), 'Binsize '+strn(lf.binsize,format='(F4.2)'), $
                 'TRGB '+strn(trgblog,format='(F6.2)')+$
                 ' '+strn(trgblogerr,format='(F4.2)'), $
                 'RGB Slope '+strn(logcoef[1],format='(F5.2)')], $
                 charsize = 1.2, textcolor=[1,1,1,1], box=0, /clear ;, /right, /bottom
; ----------------------------------------------------------------------
; log response
        plot, [min(lf.mag),max(lf.mag)], color=1, $
          [min(resplog*noise),max(resplog*noise)], line=0, ps=10, $
          xsty=3, ysty=7, position=[0.504,0.12,0.908,0.32], /noerase, /nodata, $
          ytickname=replicate(' ',10), yminor=2, xminor=3, xtit='I '
        axis, yaxis=1, yminor=3, xminor=3, ysty=7, color=1
        oplot, lf.mag, resplog*noise, color=5
;       oploterror, lf.mag, resplog, respsiglog, errcolor=5
        xyouts, [0.988,0.988], [0.22,0.22], 'Log Response', $
          /normal, charthick=2, charsize=2.2, align=0.5, orient=90, color=1
        oplot, [trgblog,trgblog], [!y.crange[0],!y.crange[1]], line=2, thick=3, color=2
; ----------------------------------------------------------------------

return
end

pro trgb_smf, objname, iter=iter, binsize=binsize, minmag=minmag, maxmag=maxmag, $
               hst=hst, halo=halo, core=core, ccut=ccut, help=help
;+
; NAME:
;	TRGB_SMF
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
;	binsize : luminosity function bin size
;	minmag  : minimum luminosity function magnitude
;	maxmag  : maximum luminosity function magnitude
;	hst  	: keyword specifying an HST object
;	halo 	: use the objname_starlist_halo.dat photometry file
;	          created in TRGB_REGIONS
;	core 	: use the objname_starlist_core.dat photometry file
;	          created in TRGB_REGIONS
;	ccut	: use the objname_starlist_ccut.dat photometry file
;		  created in TRGB_COLORCUT
;	help	: print out all the syntax of the program (all the
;		  keywords) 
;
; OUTPUTS:
;	A text file of the results of the bootstrap can be written to
;	/deep1/ioannis/trgb/results with the suffix "_edge.txt".  A
;	postscript of the luminosity function with the detected tip
;	can be written to /deep1/ioannis/trgb/plots with the suffix
;	"_edge.ps".  
;
; PROCEDURE:
;	The TRGB is detected by searching for the maximum in the
;	Poisson-weighted first derivative of both the linear and the
;	logarithmic luminosity function.  The starlist is then
;	bootstrap resampled to determine the error in the detection.
;	An iterative 3*sigma-clipping is utilized to reject outlier
;	detections (DJS_ITERSTAT).
;
; COMMENTS:
;
; PROCEDURES USED:
;	COLORTABLE1, TRGB_READATA, TRGB_LFUNCTION, BOOT_ARRAYS,
;	BOOT_RESPONSE, DJS_ITERSTAT, TRGB_DATAPATH(), PLOTSYM, LEGEND,
;	PS_OPEN, PS_CLOSE
;
; MODIFICATION HISTORY:
;	John Moustakas, 2000 June 20, UCB
;
;-

	!x.ticklen=0.03
        p_thick = !p.thick & chars = !p.charsize
        charth = !p.charthick & xth = !x.thick & yth = !y.thick
        !p.thick = 2 & !p.charsize=1.8 & !p.charthick=2 & !x.thick=2 & !y.thick=2

        if keyword_set(help) then begin
            print
            print, 'Syntax : trgb_edge, objname, iter=iter, binsize=binsize, '
            print, 'minmag=minmag, maxmag=maxmag, hst=hst, halo=halo, '
            print, 'core=core, ccut=ccut, help=help'
            print
            return
        endif

        colortable1

; read in the data

        trgb_readata, objname, datapath, data, infobase, hst=hst, halo=halo, core=core, ccut=ccut

        data[3,*] = data[3,*] - infobase.a_i	; apply the extinction correction
        nstars = n_elements(data[0,*])

	if not keyword_set(iter) then iter = 200L
        if not keyword_set(binsize) then binsize = 0.02
        if not keyword_set(minmag) then minmag = min(data[3,*])
        if not keyword_set(maxmag) then maxmag = max(data[3,*])

; generate the Sakai, Madore, Freedman luminosity function

        trgb_lfunction, data[3,*], data[4,*], lf, binsize=binsize, minmag=minmag, maxmag=maxmag

        doit = call_external('/deep1/ioannis/trgb/SMF/smf_lf.so', $
                             'smf_wrap_idl', float(lf.mags), float(lf.merr), $
                             float(lf.mag), long(lf.nhist), long(nstars), $
                             smfhist, resplin, resplog)

        lf.hist = long(smfhist)

; print up some useful information

        print, 'Object    : ', infobase.truename
        print, 'Binsize   : ', strn(binsize,format='(F6.4)')
        print, 'Minmag    : ', strn(minmag,format='(F5.2)')
        print, 'Maxmag    : ', strn(maxmag,format='(F5.2)')
        print, 'Nstars    : ', strn(total(lf.hist),format='(I6)')+'/'+strn(nstars,format='(I6)')
        print, 'Iterations: ', strn(iter,format='(I5)')
        print, '-------------------------------'
        print        
        
; constrain the endpoints

        resplin[0L:1L] = 0. & resplin[lf.nhist-2L:lf.nhist-1L] = 0.
        resplog[0L:1L] = 0. & resplog[lf.nhist-2L:lf.nhist-1L] = 0.
        
; search for the maximum between lf.minmag and the peak of the LF

        getelement_vector, lf.mag, lf.minmag, x1
        getelement_vector, lf.hist, max(lf.hist), x2

        maxlin = max(resplin[x1:x2]/noise[x1:x2],linindx)
        maxlog = max(resplog[x1:x2]*noise[x1:x2],logindx)

        trgblin = (lf.mag[x1:x2])[linindx]
        trgblog = (lf.mag[x1:x2])[logindx]

        print, 'The TRGB was detected at '+strn(trgblin,format='(F6.2)')+' mag'+$
          ' in the linear and '+strn(trgblog,format='(F6.2)')+ ' mag in the log.' & print

        print, 'Generating the bootstrap star lists . . . '
        t1 = systime(/seconds)
	boot_arrays, iter, lf, bootmags, booterrs
        t2 = systime(/seconds)-t1
        print, '. . . which took '+strn(t2,format='(F6.2)')+' seconds.' & print

; bootstrap resample the linear luminosity function

        boot_response, iter, lf, bootmags, booterrs, x1, x2, trgbmaglin

; calculate the standard deviation by iterative sigma clipping

        djs_iterstat, trgbmaglin, sigrej=3., maxiter=10, mean=trgblinmean, median=trgblinmed, $
          sigma=trgblinerr, mask=linmask
        goodlin = where(linmask eq 0B, glincount)
        print, 'A total of '+strn(iter-glincount)+' linear re-samplings were rejected at 3-sigma.' & print
        
; bootstrap resample the log luminosity function
        
        boot_response, iter, lf, bootmags, booterrs, x1, x2, trgbmaglog, /log

        djs_iterstat, trgbmaglog, sigrej=3., maxiter=10, mean=trgblogmean, median=trgblogmed, $
          sigma=trgblogerr, mask=logmask
        goodlog = where(logmask eq 0B, glogcount)
        print, 'A total of '+strn(iter-glogcount)+' logarithmic re-samplings were rejected at 3-sigma.' & print

; plot up a diagnostic plot

        window, 2, xs=400, ys=400
        plotsym, 0, 1., /fill
        plot, lf.mag, histogram(trgbmaglin,bin=lf.binsize), color=7, xsty=1, charsize=1.2, ps=10, $
          xrange=[min(trgbmaglin),max(trgbmaglin)], ytit='Number', xtit='I', ysty=3, $
          title='Distribution of Results', line=0
        oplot, [trgblin,trgblin], [!y.crange[0],!y.crange[1]], line=2, color=2, thick=5
        oplot, lf.mag, histogram(trgbmaglog,bin=lf.binsize), color=3, line=0, ps=10
        oplot, [trgblog,trgblog], [!y.crange[0],!y.crange[1]], line=2, color=5, thick=5
        legend, ['Linear Dist.','Linear TRGB','Log Dist.', 'Log TRGB'], $
          lines=[0,2,0,2], color=[7,2,3,5], charsize = 1.2, textcolor=[1,1,1,1], $
          box=0, /clear, thick=2

; print the results to the screen

        print, 'Linear TRGB : ', strn(trgblin,format='(F7.4)')
        print, 'Linear Error: ', strn(trgblinerr,format='(F7.4)')
        print, 'Log TRGB    : ', strn(trgblog,format='(F7.4)')
        print, 'Log Error   : ', strn(trgblogerr,format='(F7.4)')

        window, 3, xs=800, ys=800, ypos=30
        makeplot, lf, linsmooth, logsmooth, resplin, resplog, noise, $
          trgblin, trgblinerr, trgblog, trgblogerr, infobase, nstars

        path = trgb_datapath()
        plotpath = path[3]          ; plot subdirectory

        okay = 'N' & print
        read, okay, prompt='Generate a postscript plot (Y/[N])? '
        if strupcase(okay) eq 'Y' then begin
            ps_open, plotpath+objname+'_edge', /ps_fonts
            device, /times
            makeplot, lf, linsmooth, logsmooth, resplin, resplog, noise, $
              trgblin, trgblinerr, trgblog, trgblogerr, infobase, nstars
            ps_close
        endif

        respath = path[4]          ; results subdirectory
        resfile = respath+objname+'_edge.txt'

        okay = 'N' & print
        read, okay, prompt='Write the results to a data file (Y/[N])? '
        if strupcase(okay) eq 'Y' then begin
            
            print & print, 'Writing '+resfile+' . . .' & print

            openw, lun1, resfile, /get_lun
            printf, lun1, '# Object    : ', infobase.truename
            printf, lun1, '# Date      : ', systime()
            printf, lun1, '# -------------------------------'
            printf, lun1, '# INPUT:'
            printf, lun1, '# -------------------------------'
            if keyword_set(halo) then printf, lun1, '# Keywords Used : Halo'
            if keyword_set(core) then printf, lun1, '# Keywords Used : Core'
            if keyword_set(ccut) then printf, lun1, '# Keywords Used : Ccut'
            if ((not keyword_set(halo)) and (not keyword_set(core)) and $
                (not keyword_set(ccut))) then printf, lun1, '# Keywords Used : None'
            printf, lun1, '# Binsize   : ', strn(binsize,format='(F6.4)')
            printf, lun1, '# Minmag    : ', strn(minmag,format='(F5.2)')
            printf, lun1, '# Maxmag    : ', strn(maxmag,format='(F5.2)')
            printf, lun1, '# Nstars    : ', strn(total(lf.hist),format='(I6)')+'/'+strn(nstars,format='(I6)')
            printf, lun1, '# Iterations: ', strn(iter,format='(I5)')
            printf, lun1, '# -------------------------------'
            printf, lun1, '# OUTPUT: '
            printf, lun1, '# -------------------------------'
            printf, lun1, '# Linear TRGB : ', strn(trgblin,format='(F7.4)')
            printf, lun1, '# Linear Error: ', strn(trgblinerr,format='(F7.4)')
            printf, lun1, '# Log TRGB    : ', strn(trgblog,format='(F7.4)')
            printf, lun1, '# Log Error   : ', strn(trgblogerr,format='(F7.4)')
            printf, lun1, '# -------------------------------'
            printf, lun1, '# BOOTSTRAP DATA: '
            printf, lun1, '# -------------------------------'
            printf, lun1, '# '
            printf, lun1, '#  Linear TRGB   Log TRGB  '
            printf, lun1, '# '
            for k = 0L, iter-1L do printf, lun1, trgbmaglin[k], trgbmaglog[k], format = '(4x,F7.4,6x,F7.4)'
            free_lun, lun1        

        endif                           

        while !d.window ne -1L do wdelete, !d.window ; delete all windows
        
        !x.ticklen=0
        !p.thick = p_thick & !p.charsize=chars & !p.charthick=charth
        !x.thick = xth & !y.thick=yth

return
end
