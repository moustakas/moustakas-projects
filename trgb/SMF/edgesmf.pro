pro edgesmf, objname, iter=iter, binsize=binsize, minmag=minmag, maxmag=maxmag, $
             hst=hst, halo=halo, core=core, ccut=ccut

        path = trgb_datapath()
        plotpath = path[3]+'SMF/' ; plot subdirectory
        respath = path[4]+'SMF/'  ; results subdirectory

; read in the data

        trgb_readata, objname, datapath, data, infobase, hst=hst, halo=halo, core=core, ccut=ccut
        a_i = infobase.a_i ; extinction
        nstars = n_elements(data[0,*])

	if not keyword_set(iter) then iter = 500L
        if not keyword_set(binsize) then binsize = 0.03
        if not keyword_set(minmag) then minmag = min(data[3,*])

        trgb_lfunction, data[3,*], data[4,*], dum, binsize=binsize, minmag=minmag
        getelement_vector, dum.hist, max(dum.hist), lfmax               ; LF peak
        if not keyword_set(maxmag) then maxmag = dum.mag[lfmax]

        trgb_lfunction, data[3,*], data[4,*], lf, binsize=binsize, minmag=minmag, maxmag=maxmag

        print, 'Number of stars: ', long(total(lf.hist))
        
        error = trgb_error_func(objname,datapath,lf.mag,minmag=minmag,$ ; smoothed error function
                              maxmag=maxmag,hst=hst,halo=halo,core=core)
        if (strupcase(objname) eq 'SEXTANSB') or (strupcase(objname) eq 'HOLMBERGIX') then error = 3.0*error

        print & print, 'Smoothing the luminosity function . . . '
        smf_lfunction, lf, error, phi, resplin, resplog, noise
        smf_edge, lf, phi, resplin, resplog, noise, trgblin, trgblog ; detect the TRGB

        print
        print, 'Linear TRGB: '+strn(trgblin,format='(F6.2)')+' mag'
        print, 'Log TRGB   : '+strn(trgblog,format='(F6.2)')+' mag'
        print
        
; plot it

        plot_smf, objname, lf, phi, resplin, resplog, trgblin, trgblog, infobase, nstars
        plot_smf, objname, lf, phi, resplin, resplog, trgblin, trgblog, infobase, nstars, /ps

stop

        bootsmf, iter, lf, error, trgblinboot, trgblogboot ; bootstrap resample

        smferror, infobase, lf, iter, trgblin, trgblinboot, trgblog, $ ; calculate the error
          trgblogboot, trgblinerr, trgblogerr

        dmod2d, trgblin, trgblinerr, a_i, dmodlin, sigdmodlin, dlin, sigdlin
        dmod2d, trgblog, trgblogerr, a_i, dmodlog, sigdmodlog, dlog, sigdlog

; print the results to the screen

        print, 'Linear: '
        print, '   TRGB    : ', strn(trgblin,format='(F7.4)')
        print, '   Error   : ', strn(trgblinerr,format='(F7.4)')
        print, '   Distance: ', strn(dlin,format='(F7.4)')
        print, '   Error   : ', strn(sigdlin,format='(F7.4)')
        print, 'Log: '
        print, '   TRGB    : ', strn(trgblog,format='(F7.4)')
        print, '   Error   : ', strn(trgblogerr,format='(F7.4)')
        print, '   Distance: ', strn(dlog,format='(F7.4)')
        print, '   Error   : ', strn(sigdlog,format='(F7.4)')

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
        printf, lun1, '# Linear: '
        printf, lun1, '#    TRGB             : ', strn(trgblin,format='(F7.4)')
        printf, lun1, '#    Error            : ', strn(trgblinerr,format='(F7.4)')
        printf, lun1, '#    Distance Modulus : ', strn(dmodlin,format='(F7.4)')
        printf, lun1, '#    Error            : ', strn(sigdmodlin,format='(F7.4)')
        printf, lun1, '#    Distance         : ', strn(dlin,format='(F7.4)')
        printf, lun1, '#    Error            : ', strn(sigdlin,format='(F7.4)')
        printf, lun1, '# Log: '
        printf, lun1, '#    TRGB             : ', strn(trgblog,format='(F7.4)')
        printf, lun1, '#    Error            : ', strn(trgblogerr,format='(F7.4)')
        printf, lun1, '#    Distance Modulus : ', strn(dmodlog,format='(F7.4)')
        printf, lun1, '#    Error            : ', strn(sigdmodlog,format='(F7.4)')
        printf, lun1, '#    Distance         : ', strn(dlog,format='(F7.4)')
        printf, lun1, '#    Error            : ', strn(sigdlog,format='(F7.4)')
        printf, lun1, '# -------------------------------'
        printf, lun1, '# BOOTSTRAP DATA: '
        printf, lun1, '# -------------------------------'
        printf, lun1, '# '
        printf, lun1, '#  Linear TRGB   Log TRGB  '
        printf, lun1, '# '
        for k = 0L, iter-1L do printf, lun1, trgblinboot[k], trgblogboot[k], format = '(4x,F7.4,6x,F7.4)'
        free_lun, lun1        

return
end
