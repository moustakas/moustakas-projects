pro streams_build_psfs
; jm13sep04siena - build the PSF in each band; for now use Molino's PSFs

; specifiy the filters and some other handy info    
    filt = clash_filterlist(short=short,instr=instr,$
      weff=weff,zpt=zpt,/dropbluest)
    nfilt = n_elements(filt)

    outpath = streams_path(/psfs)
    npix = 61
    
    psfpath = getenv('CLASH_ARCHIVE')+'/PSFs/AMolino_PSF_2013/final_PSF_models/'
;   instr1 = repstr(instr,'acs','acswfc')
    for ib = 0, nfilt-1 do begin
       psffile = psfpath+instr[ib]+'_'+short[ib]+'.psfmos.fits'
       if file_test(psffile) eq 0 then begin
; the only missing one is F555W; use a dummy link for now
          splog, 'No PSF for band '+short[ib]+'!'
          this = where(strtrim(short,2) eq 'f475w')
          spawn, 'ln -sf '+short[this]+'-bpsf.fits '+outpath+short[ib]+'-bpsf.fits', /sh
       endif else begin
          bpsf1 = mrdfits(psffile,0,hdr,/silent)
          sd = dsigma(bpsf1)

; embed into a larger postage stamp
          npix1 = (size(bpsf1,/dim))[0]
          bpsf = fltarr(npix,npix)
          bpsf = randomn(seed,npix,npix)*0.5*sd
          
          embed_stamp, bpsf, bpsf1, npix/2L-npix1/2, npix/2L-npix1/2
          bpsf = bpsf/total(bpsf) ; normalize

;         ww = where(bpsf eq 0.0,nww)
;         med = im_median(bpsf1,sigrej=2.0)
;         bpsf[ww] = med+randomn(seed,nww)*0.5*sd
;         djs_plot, bpsf[*,30], /ylog
    
;         sectors_photometry, bpsf, 0.0, 0.0, npix/2, npix/2, $
;           radius, phi, counts, n_sector=1
;         ww = where(finite(radius))
;         mge_fit_sectors, radius[ww], phi[ww], counts[ww], 0.0, ngauss=4, sol=sol
          
          dfit_mult_gauss, bpsf, 1, amp, psfsig, model=model, /quiet
;         splog, 'Hack!!'
          bpsf = model
       
;; taper             
;         xx = replicate(1.0,npix)#findgen(npix)-float(npix/2L)
;         yy = findgen(npix)#replicate(1.0,npix)-float(npix/2L)
;         rr2 = xx^2+yy^2
;    
;         bpsf = bpsf*exp(-0.5*rr2/(psfsig[0]*8.0)^2)
;         bpsf = bpsf/total(bpsf)

; write out
          sxaddpar, hdr, 'PSFSIGMA', float(psfsig[0]), ' Gaussian sigma [pixel]'

          outfile = outpath+short[ib]+'-bpsf.fits'
          splog, 'Writing '+outfile
          mwrfits, float(bpsf), outfile, hdr, /create, /silent
          mwrfits, float(model), outfile, hdr, /silent

;; variable PSF; just use a constant here
;         sz = size(bpsf,/dim)
;         outfile = path+cluster[ic]+'-'+short[ib]+'-vpsf.fits'
;         mwrfits, float(bpsf), outfile, hdr, /create, /silent
;         mwrfits, bpsf*0+1, outfile, /silent
;         mwrfits, fltarr(1,1), outfile, /silent
       endelse
    endfor

return
end
    
