pro streams_qaplots
; jm13sep17siena - make some QAplots for my CLASH talk

    sample = rsex(streams_path(/propath)+'streams_sample.sex')
    splog, 'IGNORING A2261!!!'
    keep = where(strtrim(sample.shortname,2) ne 'a2261')
    sample = sample[keep]
    ncl = n_elements(sample)

    filt1 = clash_filterlist(short=short1,instr=instr1,$
      weff=weff1,zpt=zpt1,/dropbluest)

; QAplot of my preliminary F160W models
    dimagepath = streams_path(/dimage)+'macs1206/'
    image = mrdfits(dimagepath+'macs1206-f160w.fits.gz',0,/silent)
    model = mrdfits(dimagepath+'macs1206-f160w-model.fits.gz',0,/silent)

    psfile = streams_path()+'qa_macs1206_dimage.ps'
;   im_plotconfig, 11, pos, psfile=psfile, yspace=[0.05,0.05], width=3.5, $
;     height=3.5*[1,1,1], ymargin=[0.2,0.2], xmargin=[0.2,0.2], $
;     ypage=11.0, xpage=3.9

    im_plotconfig, 13, pos, psfile=psfile, xmargin=[0.1,0.1], $
      xspace=[0.01,0.01], width=2.76*[1,1,1], height=2.76, $
      ymargin=[0.2,0.2]
    
    cgimage, image, /keep, stretch=3, clip=3, minvalue=-1.0, $
      maxvalue=5.0, position=pos[*,0], /negative
    cgimage, model, /keep, stretch=3, clip=3, minvalue=-1.0, $
      maxvalue=5.0, position=pos[*,1], /noerase, /negative
    cgimage, image-model, /keep, stretch=3, clip=3, minvalue=-1.0, $
      maxvalue=5.0, position=pos[*,2], /noerase, /negative

    im_plotconfig, psfile=psfile, /psclose, /pdf
    
    

stop
    
; gather everything    
    skypath = streams_path(/skysub)
    ff = file_search(skypath+'skyinfo-*.fits.gz',count=nsky)

; plot the surface brightness profiles
    psfile = streams_path()+'qa_bcg_f110w_f160w.ps'
    im_plotconfig, 0, pos, psfile=psfile, height=5.0
    djs_plot, [0], [0], /nodata, xrange=[0.3,150.0], yrange=[-0.3,1.3], $
      xsty=1, ysty=1, xtitle='Semi-Major Axis (kpc)', $
      ytitle='F110W - F160W (AB mag)', position=pos, /xlog

;   for ic = 0, 0 do begin
    for ic = 0, ncl-1 do begin
       cluster = strtrim(sample[ic].shortname,2)

; read the skyinfo structure to get the filters
       skyinfo = mrdfits(streams_path(/skysub)+'skyinfo-'+$
         cluster+'.fits.gz',1,/silent)
       short = strtrim(skyinfo.band,2)
       nfilt = n_elements(skyinfo)
       
       apfile = streams_path(/bcg)+cluster+'/'+cluster+'-apphot.fits.gz'
       splog, 'Reading '+apfile
       phot = mrdfits(apfile,1,/silent)

       f110w = where(short eq 'f110w')
       f160w = where(short eq 'f160w')
       color = phot.abmag[f110w,*]-phot.abmag[f160w,*]
       djs_oplot, phot.radius, color, line=0, psym=-8
    endfor
    im_plotconfig, psfile=psfile, /psclose, /pdf

stop    
    
; surface brightness limit    
    psfile = streams_path()+'qa_sblimit.ps'
    im_plotconfig, 0, pos, psfile=psfile, height=5.0
    djs_plot, [0], [0], /nodata, xrange=[0.3,1.6], yrange=[27,23.0], $
      xsty=1, ysty=1, xtitle='Filter Wavelength (\mu'+'m)', $
      ytitle='SB Limit (1\sigma, mag arcsec^{-2})', $
      position=pos
    djs_oplot, !x.crange, [0,0], line=0, color='grey'
    for ii = 0, nsky-1 do begin
       sky = mrdfits(ff[ii],1,/silent)
       for jj = 0, n_elements(sky)-1 do begin
          weff = weff1[where(strtrim(sky[jj].band,2) eq strtrim(short1,2))]
          plots, weff/1D4, sky[jj].sblimit, psym=symcat(16), symsize=1.0
;         oploterror, weff/1D4, sky[jj].mode, sky[jj].sigma, psym=symcat(16), symsize=0.5
       endfor
    endfor
    im_plotconfig, psfile=psfile, /psclose, /pdf

; sky mode    
    psfile = streams_path()+'qa_skymode.ps'
    im_plotconfig, 0, pos, psfile=psfile, height=5.0, xmargin=[1.3,0.4], width=6.8
    djs_plot, [0], [0], /nodata, xrange=[0.3,1.6], yrange=[-0.6,0.2], $
      xsty=1, ysty=1, xtitle='Filter Wavelength (\mu'+'m)', $
      ytitle='Sky Mode (10^{-12} erg s^{-1} cm^{-2} Hz^{-1})', $
      position=pos
    djs_oplot, !x.crange, [0,0], line=0, color='grey'
    for ii = 0, nsky-1 do begin
       sky = mrdfits(ff[ii],1,/silent)
       for jj = 0, n_elements(sky)-1 do begin
          weff = weff1[where(strtrim(sky[jj].band,2) eq strtrim(short1,2))]
          plots, weff/1D4, sky[jj].mode, psym=symcat(16), symsize=1.0
;         oploterror, weff/1D4, sky[jj].mode, sky[jj].sigma, psym=symcat(16), symsize=0.5
       endfor
    endfor
    im_plotconfig, psfile=psfile, /psclose, /pdf
    
stop    

return
end
    
