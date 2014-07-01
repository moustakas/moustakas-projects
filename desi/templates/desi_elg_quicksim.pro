pro desi_elg_quicksim, dosim=dosim, qaplot=qaplot, debug=debug
; jm14jun23siena - quick simulation of [OII] S/N

    version = desi_deep2_template_version()
    templatepath = getenv('IM_PROJECTS_DIR')+'/desi/templates/'

    simfile = templatepath+'elg_quicksim.fits'

    fluxscale = 1D17
    zref = 1.189
    dl = dluminosity(zref,/cm)
    
    infile = getenv('DESIMODEL')+'/data/spectra/sn-spec-elg-z1.189.dat'
    desi_quicksim, model='elg', infile=infile, simdat=refsim, /silent
    get_element, refsim.wavelength, 3727.0*(1+zref)+20.0*[-1,1], indx
    refoii_snr = sqrt(total(refsim.snvec[indx[0]:indx[1]]^2))
    splog, 8D-17, 70.0, refoii_snr
       
; do the simulation    
    if keyword_set(dosim) then begin
       templatefile = templatepath+version+'/desi_deep2elg_templates_v1.1.fits.gz'
       splog, 'Reading templates'
       templates = mrdfits(templatefile,0,hdr)
       tempinfo = mrdfits(templatefile,1)
       wave = 10D^make_wave(hdr)
       ntemp = n_elements(tempinfo)

       sim = struct_addtags(struct_trimtags(tempinfo,$
         select=['objno','sigma_kms','oii_3727','oii_3727_ew']),$
         replicate({oii_snr: -1.0},ntemp))

stop       
       
       for ii = 0, 1000 do begin
;      for ii = 0, ntemp-1 do begin
          print, format='("Template=",I4.4,"/",I4.4,A5,$)', ii, ntemp, string(13b)

          factor = dluminosity(tempinfo[ii].z,/cm)/dl
          desi_quicksim, model='elg', wave*(1+zref), $
            fluxscale*factor^2*templates[*,ii]/(1+zref), simdat=sim1, /silent
          get_element, sim1.wavelength, 3727.0*(1+zref)+20.0*[-1,1], indx
          sim[ii].oii_snr = sqrt(total(sim1.snvec[indx[0]:indx[1]]^2))
          splog, sim[ii].oii_3727, sim[ii].sigma_kms, sim[ii].oii_snr
          
          if keyword_set(debug) then begin
             djs_plot, wave*(1+zref), fluxscale*templates[*,ii]/(1+zref), $
               xrange=3727.0*(1+zref)+50.0*[-1,1], xsty=1, ysty=1
             djs_oplot, sim1.wavelength, sim1.flux, color='blue'

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

       djs_plot, [0], [0], /nodata, position=pos, $
         xsty=1, ysty=1, /xlog, /ylog, xrange=[0.1,500], $
         yrange=[1,500], xtitle='[O II] \lambda3727 ('+$
         '10^{-17} erg s^{-1} cm^{-2})', ytitle='Total [O II] S/N'
       djs_oplot, fluxscale*sim.oii_3727, sim.oii_snr, $
         psym=symcat(16), symsize=0.8
       djs_oplot, fluxscale*[8D-17], [refoii_snr], psym=7, $
         color='red', symsize=3

       im_legend, 'spec-elg-z1.189', /left, /top, box=0, psym=7, $
         color='red', margin=0
       
       im_plotconfig, psfile=psfile, /psclose, /pdf
    endif

stop    
    
return
end
    
