pro likerror, infobase, lf, iter, trgbboot, widthboot, betaboot, trgberr, widtherr, betaerr

; calculate the error using the cumulative distribution function

	x1 = long(0.1827*iter) ; 18.27%
        x2 = long(0.50*iter)   ; 50%

	trgbsrt = sort(trgbboot)
        widthsrt = sort(widthboot)
        betastr = sort(betaboot)

        trgberr = (trgbboot[trgbsrt])[x2] - (trgbboot[trgbsrt])[x1]
        widtherr = (widthboot[widthsrt])[x2] - (widthboot[widthsrt])[x1]
        betaerr = (betaboot[betasrt])[x2] - (betaboot[betasrt])[x1]
        
        window, 0, xs=450, ys=450
        plothist, trgbboot, bin=0.01, color=7, $ ; xrange=[trgblin-0.5,trgblin+0.5]
          ytit='Number', xtit='I', ysty=3, title='Distribution of Results', line=0, xsty=3

        path = trgb_datapath()
        ps_open, path[3]+'errors/'+infobase.object+'_lhood_error', /ps_fonts
        device, /times, /inches
        plothist, trgbboot, bin=0.01, color=7, $ ; xrange=[trgblin-0.5,trgblin+0.5]
          ytit='Number', xtit='I', ysty=3, title=infobase.truename, line=0, xsty=3
        oplot, [trgblin,trgblin], [!y.crange[0],!y.crange[1]], line=0, color=2, thick=3
        ps_close
        
return
end


pro lhood, objname, iter=iter, binsize=binsize, minmag=minmag, maxmag=maxmag, $
           alpha=alpha, hst=hst, halo=halo, core=core, ccut=ccut
;+
; NAME:
;	TRGB_MLIKE
;
; PURPOSE:
;	Detect the tip of the red-giant branch.
;
; INPUTS:
;	objname : galaxy name
;
; OPTIONAL INPUTS:
;	iter  : number of bootstrap iterations
;
; KEYWORD PARAMETERS:
;
; COMMENTS:
;	Minmag, maxmag and binsize need to be picked intelligently.
;
; PROCEDURES USED:
;	COLORTABLE1, TRGB_READATA, TRGB_LFUNCTION, BOOT_ARRAYS,
;	TRGB_EDGE_DETECT, FORM_G(), DJS_ITERSTAT, STRN(),
;	TRGB_DATAPATH(), PS_OPEN, PS_CLOSE
;
; MODIFICATION HISTORY:
;-

; read in the data

	trgb_readata, objname, datapath, data, infobase, hst=hst, halo=halo, core=core, ccut=ccut

;       data[3,*] = data[3,*] - infobase.a_i	; apply the extinction correction
        nstars = n_elements(data[0,*])

	if not keyword_set(iter) then iter = 100L
	if not keyword_set(binsize) then binsize = 0.01
        if not keyword_set(minmag) then minmag = 20.
        if not keyword_set(alpha) then alpha = 0.30 ; slope of the RGB

; generate a luminosity function

        trgb_lfunction, data[3,*], data[4,*], dum, binsize=binsize, minmag=minmag
        getelement_vector, dum.hist, max(dum.hist), lfmax		; LF peak
        if not keyword_set(maxmag) then maxmag = dum.mag[lfmax]

        trgb_lfunction, data[3,*], data[4,*], lf, binsize=binsize, minmag=minmag, maxmag=maxmag

; print up some useful information

        print, 'Object    : ', infobase.truename
        print, 'Binsize   : ', strn(binsize,format='(F6.4)')             
        print, 'Minmag    : ', strn(minmag,format='(F5.2)')
        print, 'Maxmag    : ', strn(maxmag,format='(F5.2)')
        print, 'Nstars    : ', strn(total(lf.hist),format='(I6)')+'/'+strn(nstars,format='(I6)')
        print, 'Iterations: ', strn(iter,format='(I5)')
        print, 'alpha     : ', strn(alpha,format='(F5.2)')
        print, '-------------------------------'
        print        
        
; create the TRGB width array, cutfactor, and the bright end LF parameter, beta

        cutarray = findgen(3)/20.0+0.05
;       cutarray = findgen(10)/10.+0.1 ; discontinuity amplitude
        betarray = [0.4,0.45,0.5]
;       betarray = findgen(10)/10.+0.1 ; bright end LF slope
;       cutarray = [0.76]
;       betarray = [0.42]
        nbeta = n_elements(betarray)

; generate the smoothed error function

        error = trgb_error_func(objname,datapath,lf.mag,minmag=minmag,$
                                maxmag=maxmag,hst=hst,halo=halo,core=core)
        
        ewidth = round(error*1.44/lf.binsize) ; cube root of 3 to get correct width sigma
        eplus = (indgen(lf.nhist)+ewidth) < (lf.nhist-1L)
        a1 = where(eplus eq lf.nhist-1L,x2)
        eminus = (indgen(lf.nhist)-ewidth) > 0L
        eminus[a1[0]:a1[x2-1L]] = eminus[a1[0]]
        nsmooth = eplus-eminus+1L

; pack everything into a structure

        error_data = {error: float(error), ewidth: long(ewidth), $
                      eplus: long(eplus), eminus: long(eminus), $
                      nsmooth: long(nsmooth), a1: long(a1), x2: long(x2)}
        
; find the TRGB with the original starlist
        
        print, 'Computing the likelihood matrix of the original starlist . . . '
        t1 = systime(/seconds)
        trgb_edge_detect, lf, cutarray, betarray, error_data, lnlike, alpha=alpha
        t2 = systime(/seconds)-t1

        lnmax = max(lnlike,maxindx)

        x1 = maxindx mod lf.nhist
        z1 = maxindx/(lf.nhist*lf.nhist)
        y1 = maxindx/(lf.nhist*nbeta)

        trgb = lf.mag[x1]
        trgbwidth = cutarray[y1]
        beta = betarray[z1]
        
        print
        print, 'TRGB     : '+strn(trgb,format='(F6.2)')+' mag'
        print, 'Amplitude: '+strn(10^(-trgbwidth),format='(F6.2)')
        print, 'Beta     : '+strn(beta,format='(F6.4)')
        print

; plot the luminosity function and the initial best fit

        plot_mlike, lnlike, lf, error_data, infobase, trgb, trgbwidth, beta, nbeta, alpha=alpha

        stop

;; bootstrap resample
;
;        trgbboot = fltarr(iter)
;        widthboot = fltarr(iter)
;        betaboot = fltarr(iter)
;
;        boot_arrays, iter, lf, bootmags, booterrs
;
;        print, 'The bootstrap should take about '+strn(t2*iter/60.,forma='(F6.2)')+' minutes.' & print
;        for k = 0L, iter-1L do begin
;
;            trgb_lfunction, bootmags[*,k], booterrs[*,k], lfboot, $
;              binsize=binsize, minmag=minmag, maxmag=maxmag
;            print, 'Iteration '+strn(k+1)+' . . . computing the likelihood matrix.'
;            trgb_edge_detect, lfboot, cutarray, betarray, error_data, likeboot, alpha=alpha
;
;            maxboot = max(likeboot,indxboot)
;
;            x1boot = indxboot mod lfboot.nhist
;            z1boot = indxboot/(lfboot.nhist*lfboot.nhist)
;            y1boot = indxboot/(lfboot.nhist*nbeta)
;
;            trgbboot[k] = lfboot.mag[x1boot]
;            widthboot[k] = cutarray[y1boot]            
;            betaboot[k] = betarray[z1boot]
;            
;        endfor
;
;; calculate the error 
;        
;        likerror, infobase, lf, iter, trgbboot, widthboot, betaboot, trgberr, widtherr, betaerr
;
;; KS test
;
;;       ksone, lf.mags, 'ks_function', kstat, signif, /plot
;
;; print the solution
;
;        print, 'TRGB      : ', strn(trgb,format='(F7.4)')
;        print, 'Error     : ', strn(trgberr,format='(F7.4)')
;        print, 'TRGB Amp. : ', strn(trgbwidth,format='(F7.4)')
;        print, 'Error     : ', strn(widtherr,format='(F7.4)')
;        print, 'Beta      : ', strn(beta,format='(F5.2)')
;        print, 'Error     : ', strn(betaerr,format='(F5.2)')
;;       print, 'KS        : ', strn(signif,format='F5.2)')
;        
;; write the salient results 
;
;        path = trgb_datapath()
;        plotpath = path[3]          ; plot subdirectory
;        respath = path[4] ; results subdirectory
;
;; write all the results to a structure
;
;        strucfile = respath+infobase.object+'_lhood.dat'
;        print & print, 'Writing '+strucfile+' . . .' & print
;
;        mlike = {truename: string(infobase.truename), mags: float(reform(lf.mags)), $
;                 merr: float(reform(lf.merr)), minmag: float(minmag), maxmag: float(maxmag), $
;                 mag: float(lf.mag), cutarray: float(cutarray), betarray: float(betarray), $
;                 lnlike: float(lnlike), trgb: float(trgb), trgberr: float(trgberr), $
;                 trgbwidth: float(trgbwidth), widtherr: float(widtherr), alpha: float(alpha), $
;                 beta: float(beta), betaerr: float(betaerr), nstars: long(nstars), $
;                 iter: long(iter), maxindx: long(maxindx), date: string(systime())}                 
;        
;        swrite, mlike, strucfile
;        
;        resfile = respath+infobase.object+'_lhood.txt'
;        print & print, 'Writing '+resfile+' . . .'
;        
;        openw, lun1, resfile, /get_lun
;        printf, lun1, '# Object    : ', infobase.truename
;        printf, lun1, '# Date      : ', systime()
;        printf, lun1, '# -------------------------------'
;        printf, lun1, '# INPUT:'
;        printf, lun1, '# -------------------------------'
;        if keyword_set(halo) then printf, lun1, '# Keywords Used : Halo'
;        if keyword_set(core) then printf, lun1, '# Keywords Used : Core'
;        if keyword_set(ccut) then printf, lun1, '# Keywords Used : Ccut'
;        if ((not keyword_set(halo)) and (not keyword_set(core)) and $
;            (not keyword_set(ccut))) then printf, lun1, '# Keywords Used : None'
;        printf, lun1, '# Binsize   : ', strn(binsize,format='(F6.4)')
;        printf, lun1, '# Minmag    : ', strn(minmag,format='(F5.2)')
;        printf, lun1, '# Maxmag    : ', strn(maxmag,format='(F5.2)')
;        printf, lun1, '# Nstars    : ', strn(total(lf.hist),format='(I6)')+'/'+strn(nstars,format='(I6)')
;        printf, lun1, '# Iterations: ', strn(iter,format='(I5)')
;        printf, lun1, '# alpha     : ', strn(alpha,format='(F5.2)')
;        printf, lun1, '# -------------------------------'
;        printf, lun1, '# OUTPUT: '
;        printf, lun1, '# -------------------------------'
;        printf, lun1, '# TRGB      : ', strn(trgb,format='(F7.4)')
;        printf, lun1, '# Error     : ', strn(trgberr,format='(F7.4)')
;        printf, lun1, '# TRGB Amp. : ', strn(trgbwidth,format='(F7.4)')
;        printf, lun1, '# Error     : ', strn(widtherr,format='(F7.4)')
;        printf, lun1, '# Beta      : ', strn(beta,format='(F5.2)')
;        printf, lun1, '# Error     : ', strn(betaerr,format='(F5.2)')
;;       printf, lun1, '# KS        : ', strn(signif,format='(F5.2)')
;        printf, lun1, '# -------------------------------'
;        printf, lun1, '# BOOTSTRAP DATA: '
;        printf, lun1, '# -------------------------------'
;        printf, lun1, '# '
;        printf, lun1, '#  TRGB   TRGB Width   Beta '
;        printf, lun1, '# '
;        for k = 0L, iter-1L do printf, lun1, trgbboot[k], widthboot[k], betaboot[k], $
;          format = '(2x,F7.4,2x,F7.4,4x,F5.2)'
;        free_lun, lun1        


return
end


