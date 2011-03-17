pro deep2_ifitspec, zcat, plinefit, suffix=suffix1, debug=debug, $
  write=write, postscript=postscript
; jm07jun04nyu - based on AGES_IFITSPEC_UNFLUXED
; jm07sep27nyu - reorganized to fit with DEEP2_SPECFIT

    if (n_elements(zcat) eq 0L) then begin
       splog, 'ZCAT should be passed via DEEP2_SPECFIT.'
       return
    endif
    ngalaxy = n_elements(zcat)
    
    if (n_elements(plinefit) ne 0L) then delvarx, plinefit ; NOTE!

    get_juldate, jd
    mjdstr = string(long(jd-2400000L),format='(I5)')

    if (n_elements(suffix1) eq 0L) then suffix = mjdstr else suffix = mjdstr+'_'+suffix1

    psfile = suffix+'_specfit.ps'
    logfile = suffix+'_specfit.log'
    specdatafile = suffix+'_specdata.fits'

    if keyword_set(write) then begin
       splog, filename=logfile
       splog, 'Log file '+logfile+' opened '+systime()
       debug = 0L
    endif
    if keyword_set(postscript) then begin
       splog, 'Opening postscript file '+psfile
       dfpsplot, psfile, /color, /landscape
       debug = 0L
    endif

    if keyword_set(debug) then begin
       im_window, 0, xratio=0.55, yratio=0.6
       postthick1 = 1.0
       postthick2 = 2.0
       postthick3 = 1.0
    endif else begin
       postthick1 = 1.0
       postthick2 = 3.0
       postthick3 = 1.0
    endelse

    splog, 'IDL version: ' + string(!version,format='(99(A," "))')
    spawn, 'uname -a', uname
    splog, 'UNAME: '+uname[0]
    
    splog, 'Current datapath is ', cwd()
    
; set some datapaths, constants, and fitting parameters    
    
    datapath = deep2_path(/dr3)
    specfitpath = deep2_path(/specfit)
    
    light = 2.99792458D5        
    fwhm2sig = 2D0*sqrt(2.0*alog(2.0))
    specres = 0.56*fwhm2sig     ; [FWHM, from Willmer et al.]

    maskwidth = 5.0 ; 50.0 this should stay smallish if using b-spline iterfit
    qsowidth = 50.0 ; 200.0
;   medwidth = 150.0 & telluric = 0
    medwidth = 250.0 & telluric = 1
    vmaxshift = 100.0  ; [km/s]
    sigmax = 500.0     ; [km/s]
    qso_sigmax = 1E4   ; [km/s]
    
    nmaxiter = 1L
    
; read the emission line parameters and constraints

    if n_elements(linefile) eq 0L then linefile = 'elinelist.dat' ; 'elinelist_subset.dat'

    linepars = read_linepars(linepath=specfitpath,linefile=linefile)
    nline = n_elements(linepars)

    zindex = linepars.zindex
    windex = linepars.windex

; loop on each spectrum    

    stime0 = systime(1)
;   for igal = 55, 55 do begin
    for igal = 0L, ngalaxy-1L do begin

       galaxy = strtrim(zcat[igal].galaxy,2)
;      print, format='("Fitting galaxy ",I0,"/",I0,".",A10,$)', igal+1, ngalaxy, string(13b)
       splog, '##########'
       splog, format='("Fitting ",A0,", galaxy ",I0,"/",I0,".",A10,$)', $
         galaxy, igal+1, ngalaxy, string(13b)
       
       spec1 = mrdfits(datapath+strtrim(zcat[igal].file,2),1,/silent)
       spec2 = mrdfits(datapath+strtrim(zcat[igal].file,2),2,/silent)

       flux = [spec1.spec,spec2.spec]
       ivar = [spec1.ivar,spec2.ivar]
       wave = [spec1.lambda,spec2.lambda]
       blue = [spec1.spec*0+1,spec2.spec*0]
       red = [spec1.spec*0,spec2.spec*0+1]

       good = where(ivar gt 0.0,ngood)
       wave = wave[good]
       flux = flux[good]
       ivar = ivar[good]
       blue = blue[good]
       red  = red[good]
       ferr = 1.0/sqrt(ivar)
       npix = n_elements(wave)

       z = zcat[igal].z
       ziter = z
;      splog, 'Galaxy '+galaxy+', z = '+strtrim(string(ziter,format='(G0.0)'),2)

       for iter = 0L, nmaxiter do begin ; iterate the line-fitting

          restwave = wave/(1.0+ziter)
          restflux = flux*(1.0+ziter)
          restivar = ivar/(1.0+ziter)^2.0
          restferr = ferr*(1.0+ziter)
          restspecres = specres/(1.0+ziter)

          snr = median(restflux/restferr)
          
          lineres = replicate(restspecres,nline)          ; instr. resolution [FWHM, Angstrom, rest]
          elineres = light*lineres/linepars.wave/fwhm2sig ; instr. width [sigma, km/s, rest]

; mask emission lines; replace pixels around broad QSO emission lines
; with median values, otherwise the continuum smoothing doesn't
; work very well; setting COSMIC=1 frequently rejects the "red leak" 
          
          restmask = emission_mask(wave,spectrum=flux,good=good,bad=bad,z=ziter,$
            width=maskwidth,bluewidth=0.5*maskwidth,qsowidth=qsowidth,$
            sigrej_cosmic=3.0,qmask=qmask,tellmask=tellmask,qsomask=1,telluric=telluric,$
            bluemask=0,cosmic=0,absmask=1)

          invvar = restivar * restmask; * qmask * tellmask
          sset = bspline_iterfit(restwave,restflux,invvar=invvar,$
            bkspace=medwidth,nord=3.0,lower=0.5,upper=0.5) ;,/silent)
          continuum = bspline_valu(restwave,sset)

          espectrum = restflux - continuum
          espectrum_invvar = restivar ; could include continuum subtraction error here
          
          normconst = max(espectrum) ; normalize
          espectrum = espectrum/normconst
          espectrum_invvar = espectrum_invvar*normconst^2.0

; fit the emission lines; set the inverse variance of pixels in the
; telluric bands to zero

          igood = where(espectrum_invvar gt 0.0)
          for iline = 0L, nline-1L do begin
             if (strmatch(windex[iline],'*qso*',/fold)) then $
               wavewindow = fwhm2sig*linepars[iline].wave*5E3/light else $ ; QSO line
                 wavewindow = fwhm2sig*linepars[iline].wave*5E2/light ; other emission line
             if (linepars[iline].wave gt (min(restwave[igood])+2*lineres[iline])) and $
               (linepars[iline].wave lt (max(restwave[igood])-2*lineres[iline])) then begin
                eindx1 = where((restwave gt linepars[iline].wave-wavewindow) and $
                  (restwave lt linepars[iline].wave+wavewindow),neindx1)
;               djs_plot, restwave[eindx1], espectrum[eindx1], ps=10
;               djs_oplot, linepars[iline].wave*[1,1], !y.crange, color='red'
;               print, wavewindow & cc = get_kbrd(1)
             endif else neindx1 = 0L
             if (neindx1 ne 0L) then if (n_elements(eindx) eq 0L) then $
               eindx = eindx1 else eindx = cmset_op(eindx,'OR',eindx1)
          endfor

          if (n_elements(eindx) eq 0L) then eindx = lindgen(npix)

          if (iter eq 0L) then begin
             sigguess = replicate(100.0/light/alog(10.0),nline) ; [log-Angstrom]
             fvalue = replicate(1.0,nline)
          endif

;         splog, 'Fitting the emission-line spectrum.'
          t0 = systime(1)
          linefit = ilinefit(espectrum[eindx],restwave[eindx],linepars.wave,elineres,$
            invvar=espectrum_invvar[eindx],linename=linepars.line,zindex=zindex,$
;           invvar=espectrum_invvar[eindx]*tellmask[eindx],linename=linepars.line,zindex=zindex,$
            windex=windex,findex=findex,fvalue=fvalue,zguess=0.0D,sigmax=sigmax,qso_sigmax=qso_sigmax,$
            sigguess=sigguess,specfit=specfit,vmaxshift=vmaxshift,silent=0,_extra=extra,$
            linefit_chi2=linefit_chi2,linefit_dof=linefit_dof,linefit_niter=linefit_niter,$
            linefit_status=linefit_status)
          time = string(systime(1)-t0,format='(G0.0)')
;         splog, 'CPU time for line-fitting '+string(systime(1)-t0,format='(G0.0)')+' seconds.'

; update the initial parameters for the next iteration and then
; continue; (remember, LINEFIT gets altered by ICOMBINE_BLENDS);
; actually (jm07mar10nyu) do *not* pass SIGGUESS on the second
; iteration because it can sometimes mess up MPFIT (especially if the
; parameter is outside of its given constraints; so we sacrifice a bit
; of speed-of-convergence for more stability) 
          
;         sigguess = linefit.linesigma/light/alog(10.0)
;         fvalue = linefit.linearea/linefit.linewave/alog(10.0)

;         goodlines = where((linefit.linearea_err gt 0.0),ngoodlines,comp=badlines,ncomp=nbadlines)
          goodlines = where((linefit.linearea gt 0.0) and (linefit.linearea_err gt 0.0),$
            ngoodlines,comp=badlines,ncomp=nbadlines)

          if (nbadlines ne 0L) then linefit[badlines].linez = ziter
          if (ngoodlines ne 0L) then linefit[goodlines].linez = linefit[goodlines].linez + ziter
          
; put back the normalization constant into the emission-line fits and
; the emission-line fluxes

          espectrum = espectrum*normconst
          espectrum_invvar = espectrum_invvar/normconst^2.0
          
          speclinefit = espectrum*0.0
          speclinefit[eindx] = specfit*normconst
          
          goodfit = where(linefit.linearea_err gt 0.0,ngoodfit)
          if (ngoodfit ne 0L) then begin
             linefit[goodfit].linearea = linefit[goodfit].linearea*normconst         ; [erg/s/cm2]
             linefit[goodfit].linearea_err = linefit[goodfit].linearea_err*normconst ; [erg/s/cm2]
          endif

; box fluxes

          goodfit = where(linefit.linebox_err gt 0.0,ngoodfit)
          if (ngoodfit ne 0L) then begin
             linefit[goodfit].linebox = linefit[goodfit].linebox*normconst         ; [erg/s/cm2]
             linefit[goodfit].linebox_err = linefit[goodfit].linebox_err*normconst ; [erg/s/cm2]
          endif

; combine the [O II] doublet (not generalized!); go back to the
; spectrum and compute the LINEBOX measurement properly

          ilinefit = icombine_blends(linefit)

          oii = where(strmatch(linefit.linename,'*OII_3727_1*',/fold) or strmatch(linefit.linename,'*OII_3727_2*',/fold),noii)
          if (noii ne 0L) then linefit = [linefit[oii],ilinefit]

          oii = where(strmatch(linefit.linename,'OII_3727',/fold),noii)
          if (noii ne 0L) then if (linefit[oii].linearea_err ne -2.0) then begin

             dlam = restwave[1]-restwave[0]
             sigwidth = linefit[oii].linesigma_total/light*linefit[oii].linewave
             indx = where((restwave ge linefit[oii].linewave - 3.0*sigwidth) and $ ; +/- 3-sigma
               (restwave le linefit[oii].linewave + 3.0*sigwidth),nindx)

             lineflux = espectrum[indx]
             lineivar = espectrum_invvar[indx]*tellmask[indx]
             denom = total(lineivar+(lineivar eq 0.0))

             linefit[oii].linebox = dlam*nindx*total(lineivar*lineflux)/denom
             linefit[oii].linebox_err = dlam*float(nindx) / sqrt(denom)

;            print, linefit[oii].linearea, linefit[oii].linebox, total(lineflux)*dlam, $
;              linefit[oii].linebox_err, linefit[oii].linearea_err

          endif
          
; measure the local continuum

;         inrange = where((linefit.linewave gt min(restwave)) and (linefit.linewave lt max(restwave)),ninrange)
;         if (ninrange ne 0L) then begin
;            linefit[inrange].linecontlevel = interpol(continuum,restwave,linefit[inrange].linewave)
;            linefit[inrange].linecontlevel_err = 0.0
;         endif
          
; compute upper limits and measure the local continuum for all the
; lines in range; this way we can replace the measured line-fluxed and
; EWs with the upper limits in post-analysis as necessary; also
; construct the EWs 

;         splog, 'Computing upper limits.'
          ulimit = where(linefit.linearea_err gt -2.0,nulimit)
          if (nulimit ne 0L) then begin
             ulinefit = linefit[ulimit]
             ulinefit = iupper_limits(restwave,restflux-speclinefit,ulinefit,$
               snrcut=1.0,glosigma=5.0,ghisigma=20.0,scale=1.0,telluric=1,$
               debug=0)
             linefit[ulimit].linecontlevel     = ulinefit.linecontlevel
             linefit[ulimit].linecontlevel_err = ulinefit.linecontlevel_err
             linefit[ulimit].linelimit         = ulinefit.linearea
          endif
          
          linefit = fill_ew(linefit)

; if there are at least three well-fitted lines then update the
; redshift and iterate; also update the redshift

          if (ngoodlines ge 2L) then begin
             zline = djs_median(linefit[goodlines].linez) 
             splog, 'Updating the redshift: z = '+strtrim(string(ziter,format='(G0.0)'),2)+$
               ' --> '+strtrim(string(zline,format='(G0.0)'),2)
             ziter = zline
          endif else begin
             zline = ziter
             continue
          endelse

;         splog, igal, iter, ngoodlines, z, ziter, ' Time = '+time+' s'
;         if (iter eq nmaxiter) then $
;           struct_print, struct_trimtags(linefit,select=['*name','linez',$
;           '*area*','*box*','*cont*','*area_ew*','*limit*']) ; ,'*sigma*'])

       endfor 

; parse the linefit structure and append; jm08mar31 - do not prepend ZCAT

       plinefit1 = parse_ilinefit(linefit)
;      plinefit1 = struct_addtags(zcat[igal],struct_addtags({linefit_chi2: linefit_chi2,$
       plinefit1 = struct_addtags(struct_addtags({linefit_chi2: linefit_chi2,$
         linefit_dof: linefit_dof, linefit_niter: linefit_niter, linefit_status: $
         linefit_status, continuum_snr: float(snr)},temporary(plinefit1)))

       if (n_elements(plinefit) eq 0L) then plinefit = plinefit1 else $
         plinefit = [plinefit,plinefit1]
       
; debugging plot       

       if (keyword_set(debug) or keyword_set(postscript)) then begin

          espectrum = restflux - continuum
          
          cstats = im_stats(restflux,sigrej=3.0)
          xrange = minmax(restwave)
          yrange = cstats.median_rej + [-4.0,5.0]*cstats.sigma_rej
          yrange[0] = yrange[0]<min(continuum)
          yrange2 = [-1.0*cstats.sigma_rej,max(speclinefit)]

          nsmooth = 5L
          xmargin = [1.0,0.4]
          
; page 1          
          
;         if keyword_set(debug) then wset, 0

          pagemaker, nx=1, ny=1, xspace=0.0, yspace=0.0, ymargin=[0.7,6.2], $
            xmargin=xmargin, position=pos1, /normal
;         pagemaker, nx=1, ny=1, xspace=0.0, yspace=0.0, ymargin=[0.7,1.1], $
;           xmargin=xmargin, position=pos1, /normal

          tmasked = where(tellmask eq 0B,ntmasked)
          title = galaxy+', z = '+strtrim(string(z,format='(F12.4)'),2)+', '+$
            'Q = '+string(zcat[igal].zquality,format='(I0)')
          
          djs_plot, [0], [0], /nodata, xsty=1, ysty=1, position=pos1, $
            xthick=postthick1, ythick=postthick1, charthick=postthick1, xtitle='Rest Wavelength (\AA)', $
            charsize=1.3, ps=10, xrange=xrange, yrange=yrange, ytitle='Flux (ADU)', title=title
          djs_oplot, restwave[where(blue)], smooth(restflux[where(blue)],nsmooth), thick=postthick3, color='blue' ; color=fsc_color('royal blue',100)
          djs_oplot, restwave[where(red)], smooth(restflux[where(red)],nsmooth), thick=postthick3, color='red'
          if (ntmasked ne 0L) then djs_oplot, restwave[tmasked], restflux[tmasked], $
            color='dark green', psym=4, symsize=0.4
          djs_oplot, restwave, continuum, line=0, thick=postthick2;, color='grey'
          djs_oplot, sset.fullbkpt[where(sset.bkmask)], bspline_valu(sset.fullbkpt[where(sset.bkmask)],sset), $
            psym=7, thick=postthick2*2.0, sym=2.0, color='purple'

;         pagemaker, nx=1, ny=2, xspace=0.0, yspace=0.0, ymargin=[0.7,1.1], $
;           xmargin=xmargin, position=pos1, /normal
;
;         tmasked = where(tellmask eq 0B,ntmasked)
;         title = galaxy+', z = '+strtrim(string(z,format='(F12.4)'),2)+', '+$
;           'Q = '+string(zcat[igal].zquality,format='(I0)')
;         
;         djs_plot, restwave, smooth(restflux,nsmooth), xsty=3, ysty=3, position=pos1[*,0], $
;           xtickname=replicate(' ',10), xthick=postthick1, ythick=postthick1, charthick=postthick1, $
;           charsize=1.3, ps=10, yrange=yrange, ytitle='Counts', title=title, color='grey'
;         if (ntmasked ne 0L) then djs_oplot, restwave[tmasked], restflux[tmasked], color='dark green', ps=4
;         djs_oplot, restwave, smooth(continuum,nsmooth), color='red', thick=postthick2
;         djs_plot, restwave, smooth(restflux-continuum,nsmooth), /noerase, xsty=3, ysty=3, position=pos1[*,1], $
;           xtitle='Rest Wavelength', xthick=postthick1, ythick=postthick1, charthick=postthick1, charsize=1.5, $
;           ps=10, yrange=yrange, ytitle='Counts', color='grey'
;         djs_oplot, restwave, smooth(speclinefit,nsmooth), color='blue', thick=postthick3, ps=10
;         djs_oplot, !x.crange, [0,0], color='red', thick=postthick2

; page 2          
          
;         if keyword_set(debug) then wset, 1

          nx = 5.0
          pagemaker, nx=nx, ny=ceil(n_elements(linefit)/nx), position=pos2, /normal, $
            xmargin=xmargin, ymargin=[5.7,0.3], yspace=0.0, xspace=0.1
;         pagemaker, nx=nx, ny=ceil(n_elements(linefit)/nx), position=pos2, /normal, $
;           xmargin=xmargin, ymargin=[0.7,0.3], yspace=0.0, xspace=0.1

          elinewave = restwave*(1.0+z)/(1.0+zline) ; emission line rest wavelength vector
          box = 40.0
          tags = tag_names(plinefit1[0])
          for k = 0L, n_elements(linefit)-1L do begin

             noerase = 1L
;            if k eq 0L then noerase = 0L else noerase = 1L
             line = inice_linename(plinefit1.linename[k])

; emission-line wavelength
             
             match = where(plinefit1.linename[k]+'_WAVE' eq tags)
             match2 = where(plinefit1.linename[k]+'_LINEZ' eq tags)
             match3 = where(plinefit1.linename[k]+'_CHI2' eq tags)
             linewave = plinefit1.(match)

             linematch = where(plinefit1.linename[k] eq tags)
             errflag = plinefit1.(linematch)[1]
             
             if errflag eq -2.0 then begin

                upper = 'Not Measured'
                label = [line,upper]

                djs_plot, [0], [0], /nodata, xsty=3, ysty=3, $
                  xthick=postthick1, ythick=postthick1, thick=postthick2, noerase=noerase, ps=10, $
                  xtickname=replicate(' ',10), ytickname=replicate(' ',10), $
                  position=pos2[*,k], yrange=yrange

             endif else begin
                
; total sigma line-width       

                match = where(plinefit1.linename[k]+'_SIGMA_TOTAL' eq tags)
                lineres = linewave*plinefit1.(match)/light ; [Angstrom]

                leftbox  = linewave-((10.0*lineres)<30.0)
                rightbox = linewave+((10.0*lineres)<30.0)
                get_element, elinewave, [leftbox,rightbox], xx
;               print, linewave, leftbox, rightbox, 10*lineres

                leftbox_zoom  = linewave-3.0*lineres
                rightbox_zoom = linewave+3.0*lineres
                get_element, elinewave, [leftbox_zoom,rightbox_zoom], xx_zoom

;               estats = im_stats(espectrum[xx_zoom[0]:xx_zoom[1]],sigrej=4.0)
;               yrange = estats.median_rej + [-3.0,5.0]*estats.sigma_rej
                yrange = fltarr(2)
                yrange[1] = max(espectrum[xx_zoom[0]:xx_zoom[1]])*2.1
                yrange[0] = -0.2*yrange[1]
                xrange = minmax(elinewave[xx[0]:xx[1]])

                djs_plot, [0], [0], /nodata, xsty=3, ysty=3, xrange=xrange, $
                  xthick=postthick1, ythick=postthick1, noerase=noerase, ps=10, $
                  xtickname=replicate(' ',10), ytickname=replicate(' ',10), $
                  position=pos2[*,k], yrange=yrange, color=fsc_color('dark grey',102)
                blueindx = where(blue[xx[0]:xx[1]]) & redindx = where(red[xx[0]:xx[1]])
                if (blueindx[0] ne -1L) then djs_oplot, (elinewave[xx[0]:xx[1]])[blueindx], $
                  (espectrum[xx[0]:xx[1]])[blueindx], thick=postthick3, ps=10, $
                  color='blue'  ; color='fsc_color('royal blue',100)
                if (redindx[0] ne -1L) then djs_oplot, (elinewave[xx[0]:xx[1]])[redindx], $
                  (espectrum[xx[0]:xx[1]])[redindx], thick=postthick3, ps=10, color='red'
                djs_oplot, elinewave[xx[0]:xx[1]], speclinefit[xx[0]:xx[1]], thick=postthick2, ps=10;, color='red'
                djs_oplot, linewave*[1,1], !y.crange, line=2, thick=postthick2;, color='navy'
                djs_oplot, !x.crange, [0,0], line=0, thick=postthick2;, color='navy'

                case errflag of
                   -3.0: label = [line,'Upper Limit']
                   else: begin
                      snr = 'S/N = '+string((plinefit1.(linematch))[0]/(plinefit1.(linematch))[1],format='(I0)')
                      chi2 = '\chi^{2} = '+strtrim(string(plinefit1.(match3),format='(F12.1)'),2)
                      label = [line,snr+', '+chi2]
                   endelse
                endcase
                
             endelse

             legend, textoidl(label), /left, /top, box=0, charsize=0.8, margin=0, $
               charthick=postthick1, clear=keyword_set(postscript)

          endfor
          
; title    
;         refpos = reform(pos1[*,0])
;         xpos = (refpos[2]-refpos[0])/2.0+refpos[0]
;         ypos = refpos[3]*1.02
;         xyouts, xpos, ypos, title, /normal, charsize=1.5, charthick=postthick1, align=0.5

; old code          
          
;          pagemaker, nx=2, ny=2, xspace=0.2, yspace=0.2, ymargin=[0.4,0.4], $
;            xmargin=[1.1,1.1], position=pos, /normal
;          wbox = 8.0
;; [O II] doublet
;          xrange = 3727+wbox*[-1,1]
;          get_element, restwave, xrange, xx
;          if (xx[0] ge npix-1L) or (xx[1] le 0L) then begin
;             djs_plot, [0], [0], /nodata, xsty=3, ysty=3, position=pos[*,0], $
;               xtitle='', xtickname=replicate(' ',10), ytickname=replicate(' ',10), $
;               xthick=postthick1, ythick=postthick1, charthick=postthick1, charsize=1.5, $
;               ps=10, xrange=xrange, ytitle='', yrange=[0,1], xtickinterval=50
;             legend, 'No measurement', /left, /bottom, box=0, charsize=1.5, charthick=postthick1
;          endif else begin
;             yrange2[1] = max(espectrum[xx[0]:xx[1]])*1.1
;             djs_plot, restwave, espectrum, xsty=3, ysty=3, position=pos[*,0], $
;               xtitle='', xtickname=replicate(' ',10), xthick=postthick1, ythick=postthick1, charthick=postthick1, charsize=1.5, $
;               ps=10, xrange=xrange, yrange=yrange2, xtickinterval=50, ytitle='Counts'
;             djs_oplot, !x.crange, [0,0], thick=postthick1
;             djs_oplot, restwave, speclinefit, color='blue', thick=postthick3, ps=10
;             if (ntmasked ne 0L) then djs_oplot, restwave[tmasked], speclinefit[tmasked], color='dark green', ps=4
;          endelse
;          legend, '[O II]', /left, /top, box=0, charsize=1.5, charthick=postthick1
;; H-beta
;          xrange = 4861+wbox*[-1,1]
;          get_element, restwave, xrange, xx
;          if (xx[0] ge npix-1L) or (xx[1] le 0L) then begin
;             djs_plot, [0], [0], /nodata, /noerase, xsty=3, ysty=11, position=pos[*,1], $
;               xtitle='', xtickname=replicate(' ',10), ytickname=replicate(' ',10), $
;               xthick=postthick1, ythick=postthick1, charthick=postthick1, charsize=1.5, $
;               ps=10, xrange=xrange, ytitle='', yrange=[0,1], xtickinterval=50
;             axis, /yaxis, ytitle='', ythick=postthick1, ysty=3, charsize=1.5, charthick=postthick1, $
;               ytickname=replicate(' ',10)
;             legend, 'No measurement', /left, /bottom, box=0, charsize=1.5, charthick=postthick1
;          endif else begin
;             yrange2[1] = max(espectrum[xx[0]:xx[1]])*1.1
;             djs_plot, restwave, espectrum, /noerase, xsty=3, ysty=11, position=pos[*,1], $
;               xtitle='', xtickname=replicate(' ',10), xthick=postthick1, ythick=postthick1, charthick=postthick1, charsize=1.5, $
;               ps=10, xrange=xrange, ytickname=replicate(' ',10), yrange=yrange2, xtickinterval=50
;             axis, /yaxis, ytitle='Counts', ythick=postthick1, ysty=3, charsize=1.5, charthick=postthick1
;             djs_oplot, !x.crange, [0,0], thick=postthick1
;             djs_oplot, restwave, speclinefit, color='blue', thick=postthick3, ps=10
;             if (ntmasked ne 0L) then djs_oplot, restwave[tmasked], speclinefit[tmasked], color='dark green', ps=4
;          endelse
;          legend, textoidl('H\beta'), /left, /top, box=0, charsize=1.5, charthick=postthick1
;; [O III] 4959
;          xrange = 4959+wbox*[-1,1]
;          get_element, restwave, xrange, xx
;          if (xx[0] ge npix-1L) or (xx[1] le 0L) then begin
;             djs_plot, [0], [0], /nodata, /noerase, xsty=3, ysty=3, position=pos[*,2], $
;               xtitle='', xtickname=replicate(' ',10), ytickname=replicate(' ',10), $
;               xthick=postthick1, ythick=postthick1, charthick=postthick1, charsize=1.5, $
;               ps=10, xrange=xrange, ytitle='', yrange=[0,1], xtickinterval=50
;             legend, 'No measurement', /left, /bottom, box=0, charsize=1.5, charthick=postthick1
;          endif else begin
;             yrange2[1] = max(espectrum[xx[0]:xx[1]])*1.1
;             djs_plot, restwave, espectrum, /noerase, xsty=3, ysty=3, position=pos[*,2], ytitle='Counts', $
;               xtitle='', xtickname=replicate(' ',10), xthick=postthick1, ythick=postthick1, charthick=postthick1, charsize=1.5, $
;               ps=10, xrange=xrange, yrange=yrange2, xtickinterval=50
;             djs_oplot, !x.crange, [0,0], thick=postthick1
;             djs_oplot, restwave, speclinefit, color='blue', thick=postthick3, ps=10
;             if (ntmasked ne 0L) then djs_oplot, restwave[tmasked], speclinefit[tmasked], color='dark green', ps=4
;          endelse
;          legend, textoidl('[O III] \lambda4959'), /left, /top, box=0, charsize=1.5, charthick=postthick1
;; [O III] 5007
;          xrange = 5007+wbox*[-1,1]
;          get_element, restwave, xrange, xx
;          if (xx[0] ge npix-1L) or (xx[1] le 0L) then begin
;             djs_plot, [0], [0], /nodata, /noerase, xsty=3, ysty=11, position=pos[*,3], $
;               xtitle='', xtickname=replicate(' ',10), ytickname=replicate(' ',10), $
;               xthick=postthick1, ythick=postthick1, charthick=postthick1, charsize=1.5, $
;               ps=10, xrange=xrange, ytitle='', yrange=[0,1], xtickinterval=50
;             axis, /yaxis, ytitle='', ythick=postthick1, ysty=3, charsize=1.5, charthick=postthick1, $
;               ytickname=replicate(' ',10)
;             legend, 'No measurement', /left, /bottom, box=0, charsize=1.5, charthick=postthick1
;          endif else begin
;             yrange2[1] = max(espectrum[xx[0]:xx[1]])*1.1
;             djs_plot, restwave, espectrum, /noerase, xsty=3, ysty=11, position=pos[*,3], $
;               xtitle='', xtickname=replicate(' ',10), xthick=postthick1, ythick=postthick1, charthick=postthick1, charsize=1.5, $
;               ps=10, xrange=xrange, ytickname=replicate(' ',10), yrange=yrange2, xtickinterval=50
;             axis, /yaxis, ytitle='Counts', ythick=postthick1, ysty=3, charsize=1.5, charthick=postthick1
;             djs_oplot, !x.crange, [0,0], thick=postthick1
;             djs_oplot, restwave, speclinefit, color='blue', thick=postthick3, ps=10
;             if (ntmasked ne 0L) then djs_oplot, restwave[tmasked], speclinefit[tmasked], color='dark green', ps=4
;          endelse
;          legend, textoidl('[O III] \lambda5007'), /left, /top, box=0, charsize=1.5, charthick=postthick1

          if keyword_set(debug) then begin
             splog, 'Press any key to continue'
             cc = get_kbrd(1)
          endif

       endif

    endfor
    splog, 'Total time = '+string((systime(1)-stime0)/60.0,format='(G0.0)')+' minutes.'

    if keyword_set(write) then begin
       splog, 'Writing '+specdatafile
       mwrfits, plinefit, specdatafile, /create
       spawn, 'gzip -f '+specdatafile, /sh
    endif
    if keyword_set(postscript) then begin
       dfpsclose
       spawn, 'gzip -f '+psfile, /sh
    endif
    
    splog, /close

return
end
    
