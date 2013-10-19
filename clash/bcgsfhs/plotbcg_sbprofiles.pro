pro plotbcg_sbprofiles

    sample = rsex(streams_path(/propath)+'streams_sample.sex')
    splog, 'IGNORING A2261!!!'
    keep = where(strtrim(sample.shortname,2) ne 'a2261')
    sample = sample[keep]
    ncl = n_elements(sample)

    pixscale = 0.065
    
; plot these bands    
    thesefilt = ['f160w','f814w','f475w']
    color = ['orange','firebrick','forest green']
    nthese = n_elements(thesefilt)
    
;   for ic = 0, 0 do begin
    for ic = 0, ncl-1 do begin
       cluster = strtrim(sample[ic].shortname,2)
       print & splog, cluster, sample[ic].z
       datapath = streams_path(/postman_bcg)+cluster+'/'

       arcsec2kpc = dangular(sample[ic].z,/kpc)/206265D ; [kpc/arcsec]
       
       galphot = mrdfits(datapath+cluster+'-ellipse-image.fits.gz',1,/silent)
       modphot = mrdfits(datapath+cluster+'-ellipse-model.fits.gz',1,/silent)

       reffilt = where(strtrim(galphot.band,2) eq 'f160w') ; reference filter
       
       xrange = [3*pixscale*arcsec2kpc,max(galphot[reffilt].radius_kpc)]
;      xrange = [pixscale*arcsec2kpc,max(galphot.radius_kpc)*1.5]
       yrange = [max(modphot.sb0fit),min(modphot.sb0fit)]
       yrange = [29,17]
       
       djs_plot, [0], [0], /nodata, position=pos, xsty=3, ysty=1, $
        yrange=yrange, xrange=xrange, /xlog, $
         xtitle='Semi-Major Axis (kpc)', ytitle='\mu (mag arcsec^{-2})'
       im_legend, strupcase(cluster), /right, /top, box=0, margin=0
       for ii = 0, nthese-1 do begin
          this = where(thesefilt[ii] eq strtrim(galphot.band,2))

;         notzero = lindgen(n_elements(galphot[this].sb0fit_ivar))
          notzero = where(galphot[this].sb0fit_ivar gt 0)
          djs_oplot, galphot[this].radius_kpc[notzero], galphot[this].sb0fit[notzero], $
            color=cgcolor('medium gray')
;         oploterror, galphot[this].radius_kpc[notzero], galphot[this].sb0fit[notzero], $
;           1.0/sqrt(galphot[this].sb0fit_ivar[notzero]), color=cgcolor('dodger blue'), $
;           psym=symcat(16), symsize=0.4

          good = where(modphot[this].sb0fit lt modphot[this].sblimit,ngood)
          splog, thesefilt[ii], max(modphot[this].radius_kpc[good])
          djs_oplot, modphot[this].radius_kpc[good], modphot[this].sb0fit[good], $
            color=cgcolor(color[ii])
;         djs_oplot, 10^!x.crange, galphot[this].sblimit*[1,1], line=5, $
;           color=cgcolor(color[ii])
       endfor

       cc = get_kbrd(1)
    endfor

stop    
    
return
end
    
