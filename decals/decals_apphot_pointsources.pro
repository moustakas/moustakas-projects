pro decals_apphot_pointsources, sdss=sdss
; jm14nov25siena - do aperture photometry of the non-SDSS point
; sources in each brick

    edrdir = getenv('HOME')+'/tmp/'
;   edrdir = getenv('DECALS_DIR')+'/release/edr/'
    qadir = getenv('HOME')+'/tmp/'
    bricklist = file_search(edrdir+'tractor/tractor-??????.fits',count=nbrick)

    band = ['g','r','z']
    nband = n_elements(band)

    rsat = 17.0 ; (conservative) magnitude above which stars are saturated

    if keyword_set(sdss) then begin
       rbright = 17.0
       rfaint = 26.0
       rbin = 0.5
    endif else begin
       rbright = 21.5
       rfaint = 26.0
       rbin = 0.25
    endelse

    nrad = 14
    radius = findgen(nrad)+1 ; pixels

    for ii = 1L, 1 do begin
;   for ii = 0L, nbrick-1 do begin
       cat = mrdfits(bricklist[ii],1,hdr)
       brick = strtrim(cat[0].brickid,2)

; photometer the SDSS-detected stars
       if keyword_set(sdss) then begin
          these = where(cat.type eq 'S' and cat.brick_primary eq 'T' and $
            cat.decam_flux[1] gt 0 and cat.decam_flux[1] lt 10.0^(-0.4*(rsat-22.5)) and $
            cat.sdss_objid ne 0)
          qafile = qadir+'qa_apphot_'+brick+'_sdss.ps'
       endif else begin
          these = where(cat.type eq 'S' and cat.brick_primary eq 'T' and $
            cat.decam_flux[1] gt 0 and cat.decam_flux[1] lt 10.0^(-0.4*(rsat-22.5)) and $
            cat.sdss_objid eq 0)
          qafile = qadir+'qa_apphot_'+brick+'.ps'
       endelse
       cat = cat[these]

; build the residual image in each band and do photometry 
       for jj = 1, 1 do begin ; just the r-band for now
;      for jj = 0, nband-1 do begin
          image = mrdfits(edrdir+'coadd/image2-'+brick+'-'+band[jj]+'.fits',0,hdr)
          invvar = mrdfits(edrdir+'coadd/image-'+brick+'-'+band[jj]+'.fits',1)
          model = mrdfits(edrdir+'coadd/model-'+brick+'-'+band[jj]+'.fits.gz',0)

          extast, hdr, astr
          pixscale = sqrt(abs(determ(astr.cd)))*3600

          flux = djs_phot(cat.x,cat.y,radius,skyrad,image,$
            invvar,calg='none',flerr=ferr,peakval=peak)
          flux_resid = djs_phot(cat.x,cat.y,radius,skyrad,image-model,$
            invvar,calg='none',flerr=ferr_resid,peakval=peak_resid)

;         djs_plot, -2.5*alog10(cat.decam_flux[1])+22.5, peak, psym=8, ysty=3
          
; compute the variance weighted mean curve-of-growth in bins of
; r-magnitude in both the data and in the model-subtracted data 
          hist = histogram(22.5-2.5*alog10(cat.decam_flux[1]),$
            locations=rbinmag,bin=rbin,reverse=rev,min=rbright+rbin,$
            max=rfaint)
          rbinmag = rbinmag-rbin/2.0
          nbin = n_elements(hist)

          out = {$
            rmag:                   0.0,$
            nstar:                    0,$
            radius:              radius,$
            cog:           fltarr(nrad)-1,$
            cog_err:       fltarr(nrad)-1,$
            cog_resid:     fltarr(nrad)-1,$
            cog_resid_err: fltarr(nrad)-1}
          out = replicate(out,nbin)
          out.rmag = rbinmag
          out.nstar = hist

          for ibin = 0, nbin-1 do begin
             if (rev[ibin+1] gt rev[ibin]) and (hist[ibin] ge 5) then begin
                indx = rev[rev[ibin]:rev[ibin+1]-1]
            
                for irad = 0, nrad-1 do begin
                   good = where(ferr[indx,irad] gt 0.0,ngood)
                   if ngood ge 5 then begin
                      out[ibin].cog[irad] = im_weighted_mean(flux[indx[good],irad],$
                        errors=ferr[indx[good],irad],wmean_err=mnerr)
                      out[ibin].cog_err[irad] = mnerr
                   endif

                   good = where(ferr_resid[indx,irad] gt 0.0,ngood)
                   if ngood ge 5 then begin
                      out[ibin].cog_resid[irad] = im_weighted_mean(flux_resid[indx[good],irad],$
                        errors=ferr_resid[indx[good],irad],wmean_err=mnerr)
                      out[ibin].cog_resid_err[irad] = mnerr
                   endif
                endfor
             endif
          endfor

          im_plotconfig, 30, pos, psfile=qafile, charsize=1.1, $
            width=2.2*[1,1,1], xspace=[0.4,0.4], xmargin=[0.9,0.1], $
            yspace=[0.05,0.05,0.05,0.05,0.05], ymargin=[0.4,1.1], $
            height=[1.5,1.5,1.5,1.5,1.5,1.5]
          for ibin = 0, nbin-1 do begin
             if ibin le 14 then xtickname = replicate(' ',10) else delvarx, xtickname
             yrange = [min(out[ibin].cog)<min(out[ibin].cog_resid),$
               max(out[ibin].cog)>max(out[ibin].cog_resid)]
             if out[ibin].cog[0] eq -1.0 then ytickname = replicate(' ',10) else $
               delvarx, ytickname
             djs_plot, [0], [0], position=pos[*,ibin], xsty=1, ysty=1, $
               xrange=[0,max(radius)], yrange=yrange*1.2, noerase=ibin gt 0, $
               xtickname=xtickname, ytickname=ytickname
             djs_oplot, !x.crange, [0,0], color='grey', line=5
             im_legend, ['<r>='+string(out[ibin].rmag,format='(F6.3)'),$
               'N='+strtrim(out[ibin].nstar,2)], /left, /top, box=0, margin=0
             if out[ibin].cog[0] ne -1.0 then begin
                oploterror, out[ibin].radius, out[ibin].cog, $
                  out[ibin].cog_err
                oploterror, out[ibin].radius, out[ibin].cog_resid, $
                  out[ibin].cog_resid_err, color=cgcolor('orange'), $
                  errcolor=cgcolor('orange')
             endif
             xyouts, 0.05, pos[1,6], align=0.5, orientation=90, $
               'Cumulative Flux (nanomaggies)', /norm, charsize=1.5
             xyouts, djs_mean([pos[0,16],pos[2,16]]), 0.05, align=0.5, $
               'Aperture radius (pixels)', /norm, charsize=1.5
          endfor 
          im_plotconfig, psfile=qafile, /psclose, /pdf

stop
       endfor

    endfor

return
end
