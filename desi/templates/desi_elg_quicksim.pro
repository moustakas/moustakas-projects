pro desi_elg_quicksim, dosim=dosim, qaplot=qaplot, debug=debug
; jm14jun23siena - quick simulation of [OII] S/N

    version = desi_deep2_template_version()
    templatepath = getenv('IM_PROJECTS_DIR')+'/desi/templates/'

    simfile = templatepath+'elg_quicksim.fits'

    fluxscale = 1D17
    light = im_light(/kms)
    zref = 1.189
    dlref = dluminosity(zref,/cm)
    
    infile = getenv('DESIMODEL')+'/data/spectra/sn-spec-elg-z1.189.dat'
    desi_quicksim, model='elg', infile=infile, simdat=refsim
    get_element, refsim.wavelength, 3727.0*(1+zref)+20.0*[-1,1], indx
    refoii_snr = sqrt(total(refsim.snvec[indx[0]:indx[1]]^2))
    splog, 8D-17, 70.0, refoii_snr
       
; do the simulation    
    if keyword_set(dosim) then begin
       templatefile = templatepath+version+'/deep2elg_templates_obs_v1.2.fits.gz'
       splog, 'Reading templates'
       templates = mrdfits(templatefile,0,hdr)
       tempinfo = mrdfits(templatefile,1)
       wave = make_wave(hdr)
       ntemp = n_elements(tempinfo)

       zobj = tempinfo.z
       fluxfactor = fluxscale*(1+zobj)/(1+zref)*(dluminosity(zobj,/cm)/dlref)^2.0
       
       sim = struct_addtags(struct_trimtags(tempinfo,$
         select=['objno','sigma_kms','oii_3727','oii_3727_ew']),$
         replicate({oii_snr: -1.0},ntemp))

;      for ii = 0, 100 do begin
       for ii = 0, ntemp-1 do begin
          print, format='("Template=",I4.4,"/",I4.4,A5,$)', ii, ntemp, string(13b)

          desi_quicksim, model='elg', wave*(1+zref)/(1+zobj[ii]), $
            templates[*,ii]*fluxfactor[ii], simdat=sim1

          wavelim = 3728.4837*(1+zref)*(1+10*[-1,1]*sim[ii].sigma_kms/light)
          get_element, sim1.wavelength, wavelim, indx
          sim[ii].oii_snr = sqrt(total(sim1.snvec[indx[0]:indx[1]]^2))
          splog, sim[ii].oii_3727, sim[ii].sigma_kms, sim[ii].oii_snr
          
          if keyword_set(debug) then begin
             djs_plot, wave*(1+zref)/(1+zobj[ii]), templates[*,ii]*fluxfactor[ii], $
               xrange=wavelim*[0.99,1.01], xsty=1, ysty=1
             djs_oplot, sim1.wavelength, sim1.flux, color='blue'
             cc = get_kbrd(1)
             
;            ploterror, sim1.wavelength, sim1.flux, 1.0/sqrt(sim1.invvar), $
;              /trad, xrange=3727.0*(1+zref)+1000.0*[-1,1], xsty=1, ysty=1
             djs_plot, sim1.wavelength, sim1.snvec, $
               xrange=3727.0*(1+zref)+1000.0*[-1,1], xsty=1, ysty=1
             djs_oplot, sim1.wavelength[indx[0]:indx[1]], $
               sim1.snvec[indx[0]:indx[1]], color='orange'
             cc = get_kbrd(1)
          endif
       endfor
       im_mwrfits, sim, simfile, /clobber
    endif

; build the QAplot    
    if keyword_set(qaplot) then begin

       sim = mrdfits(simfile+'.gz',1)
       these = where(sim.oii_snr gt 0)
       sim = sim[these]
       
       psfile = templatepath+'qa_elg_quicksim.eps'
       im_plotconfig, 0, pos, psfile=psfile, height=4.5, $
         xmargin=[1.5,0.4], width=6.6
       ss = texsyms()

;       djs_plot, [0], [0], /nodata, position=pos, $
;         xsty=1, ysty=1, /xlog, /ylog, xtitle='[O II] \lambda3727 ('+$
;         '10^{-17} erg s^{-1} cm^{-2})', ytitle=ss.sqrt+'(S/N_{[O II]})^2}', $
;         xrange=[0.1,200], yrange=[1,200]

       hogg_scatterplot, alog10(fluxscale*sim.oii_3727), alog10(sim.oii_snr), $
         xrange=alog10([0.4,200]), yrange=alog10([1,100]), /internal, $
         levels=[0.25,0.5,0.75,0.95], /nogrey, $
         position=pos, xsty=1, ysty=1, xtitle=textoidl('log_{10} ([O II] \lambda3726,29) ('+$
         '10^{-17} erg s^{-1} cm^{-2})'), ytitle=textoidl('log_{10} '+ss.sqrt+'(S/N_{[O II]})^2'), $
         /outliers, outcolor=cgcolor('grey')
;      djs_oplot, fluxscale*sim.oii_3727, sim.oii_snr, $
;        psym=symcat(6,thick=1), symsize=0.3, color=cgcolor('grey')
       djs_oplot, alog10(fluxscale*[8D-17]), alog10([refoii_snr]), $
         psym=symcat(6,thick=4), color='red', symsize=2

       im_legend, 'spec-elg-z1.189', /left, /top, box=0, psym=6, $
         color='red', margin=0
       
       im_plotconfig, psfile=psfile, /psclose, /pdf
    endif

stop    
    
return
end
    
