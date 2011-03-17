; determine the one-sigma error of the bootstrapped TRGB detections
; and generate a plot of the errors
pro smferror, infobase, lf, iter, trgblin, trgblinboot, trgblog, trgblogboot, linerr, logerr

; ----------------------------------------------------------------------
; calculate the error using the cumulative distribution function
; ----------------------------------------------------------------------

	x1 = long(0.1827*iter) ; 18.27%
        x2 = long(0.50*iter)   ; 50%
	linsrt = sort(trgblinboot)
	logsrt = sort(trgblogboot)

        linerr = (trgblinboot[linsrt])[x2] - (trgblinboot[linsrt])[x1]
        logerr = (trgblogboot[logsrt])[x2] - (trgblogboot[logsrt])[x1]
        
        window, 0, xs=450, ys=450
        plothist, trgblinboot, bin=0.01, color=7, $ ; xrange=[trgblin-0.5,trgblin+0.5]
          ytit='Number', xtit='I', ysty=3, title='Distribution of Results', line=0, xsty=3
        plothist, trgblogboot, color=3, /overplot, bin=0.01, line=2

        path = trgb_datapath()
        ps_open, path[3]+'errors/'+infobase.object+'_smf_error', /ps_fonts
        device, /times, /inches
        plothist, trgblinboot, bin=0.01, color=7, $ ; xrange=[trgblin-0.5,trgblin+0.5]
          ytit='Number', xtit='I', ysty=3, title=infobase.truename, line=0, xsty=3
        oplot, [trgblin,trgblin], [!y.crange[0],!y.crange[1]], line=0, color=2, thick=3
        plothist, trgblogboot, color=3, /overplot, bin=0.01, line=2
        oplot, [trgblog,trgblog], [!y.crange[0],!y.crange[1]], line=2, color=5, thick=3
        legend, ['Linear','Log'], lines=[0,2], box=0, charsize=1.5, color=1
        ps_close
        
return
end

pro trgb_edge_smf, objname, iter=iter, binsize=binsize, minmag=minmag, maxmag=maxmag, $
                   hst=hst, halo=halo, core=core, ccut=ccut

; read in the data

        trgb_readata, objname, datapath, data, infobase, hst=hst, halo=halo, core=core, ccut=ccut
        nstars = n_elements(data[0,*])

	if not keyword_set(iter) then iter = 500L
        if not keyword_set(binsize) then binsize = 0.03
        if not keyword_set(minmag) then minmag = min(data[3,*])
        if not keyword_set(maxmag) then maxmag = max(data[3,*])

        trgb_lfunction, data[3,*], data[4,*], lf, binsize=binsize, minmag=minmag, maxmag=maxmag

        error = trgb_error_func(objname,datapath,lf.mag,minmag=minmag,$ ; smoothed error function
                              maxmag=maxmag,hst=hst,halo=halo,core=core)
;       error = 3.*error

        smf_lfunction, lf, error, phi, resplin, resplog
        smf_edge, lf, phi, resplin, resplog, trgblin, trgblog ; detect the TRGB

        print
        print, 'Linear TRGB: '+strn(trgblin,format='(F6.2)')+' mag'
        print, 'Log TRGB   : '+strn(trgblog,format='(F6.2)')+' mag'
        print
        
; plot it

        plot_smf, lf, phi, resplin, resplog, trgblin, trgblinerr, $
          trgblog, trgblogerr, infobase, nstars


        plot_smflf, infobase.truename, lf, phi, /log
        oplot, [trgblin,trgblin], [!y.crange[0],!y.crange[1]], line=0, color=3, thick=2
        oplot, [trgblog,trgblog], [!y.crange[0],!y.crange[1]], line=2, color=5, thick=2

        print, 'Generating the bootstrap star lists . . . '
        t1 = systime(/seconds)
	boot_arrays, iter, lf, bootmags, booterrs
        t2 = systime(/seconds)-t1
        print, '. . . which took '+strn(t2,format='(F6.2)')+' seconds.' & print

; bootstrap resample

        trgblinboot = fltarr(iter)
        trgblogboot = fltarr(iter)

        for k = 0L, iter-1L do begin

            trgb_lfunction, bootmags[*,k], booterrs[*,k], lfboot, $
              binsize=binsize, minmag=minmag, maxmag=maxmag
            if k mod 10 eq 0L then print, 'Iteration '+strn(k)+'.'

            smf_lfunction, lfboot, error, phiboot, resplin, resplog
            smf_edge, lfboot, phiboot, resplin, resplog, dum1, dum2

            trgblinboot[k] = dum1
            trgblogboot[k] = dum2
            
        endfor

        smferror, infobase, lf, iter, trgblin, trgblinboot, trgblog, trgblogboot, trgblinerr, trgblogerr

; print the results to the screen

        print, 'Linear TRGB : ', strn(trgblin,format='(F7.4)')
        print, 'Linear Error: ', strn(trgblinerr,format='(F7.4)')
        print, 'Log TRGB    : ', strn(trgblog,format='(F7.4)')
        print, 'Log Error   : ', strn(trgblogerr,format='(F7.4)')

        path = trgb_datapath()
        plotpath = path[3]          ; plot subdirectory
        respath = path[4]          ; results subdirectory

; generate a postscript plot of the results

        plot_smf, lf, phi, resplin, resplog, trgblin, trgblinerr, $
          trgblog, trgblogerr, infobase, nstars

        plotfile = plotpath+objname+'_smf'
        ps_open, plotfile, /ps_fonts
        device, /times, /inches
        plot_smf, lf, phi, resplin, resplog, trgblin, trgblinerr, $
          trgblog, trgblogerr, infobase, nstars
        ps_close
        
; write the results out 

        resfile = respath+objname+'_smf.txt'

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
        for k = 0L, iter-1L do printf, lun1, trgblinboot[k], trgblogboot[k], format = '(4x,F7.4,6x,F7.4)'
        free_lun, lun1        

stop        

return
end
