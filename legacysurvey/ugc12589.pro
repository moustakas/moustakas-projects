pro ugc12589
; jm15jun24siena - make a pretty picture of the interacting galaxy UGC12589 for
; the cover of the monthly NOAO Currents issue
    
; first hack common.wcs_for_brick() so pixscale=0.065 then do
;
; python $TRACTOR_DIR/projects/desi/runbrick.py --no-write --stage image_coadds
; --radec 351.26 0.007 --width 2500 --height 2500

    band = ['g','r','z']
    dir = 'coadd/cus/custom-351260p00007/'
    
; then interpolate over the saturated star
    for ii = 0, n_elements(band)-1 do begin
       im = mrdfits(dir+'decals-custom-351260p00007-image-'+band[ii]+'.fits',0,hdr)
       fim = djs_maskinterp(im,im eq 0,iaxis=1)
       sim = gauss_smooth(fim,2.0,/edge_truncate)
       mwrfits, sim, 'ugc12589-'+band[ii]+'.fits', hdr, /create
    endfor

; and then
; python trilogy.py trilogy.in
    
    
return
end
