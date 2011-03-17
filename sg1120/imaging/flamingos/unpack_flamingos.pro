pro unpack_flamingos, unpack=unpack
; jm07jan17nyu

    datapath = flamingos_path(/mosaics)
    outpath = datapath
    
    imagelist = datapath+['as_sg1120a_ks.fits','as_sg1120b_ks.fits']
    weightlist = datapath+['expmap_sg1120a_ks.fits','expmap_sg1120b_ks.fits']
    outimagelist = outpath+file_basename(repstr(imagelist,'as_',''))
    outweightlist = repstr(outimagelist,'.fits','.weight.fits')
    
    nimage = n_elements(imagelist)
    
    for ii = 0L, nimage-1L do begin
       
       image = mrdfits(imagelist[ii],0,hdr,/silent)
       weight = mrdfits(weightlist[ii],0,/silent)

       mkhdr, outhdr, image
       sxaddpar, outhdr, 'OBJECT', 'SG1120'
;      sxaddpar, outhdr, 'EXPTIME', sxpar(hdr,'EXP_TIME'), ' exposure time (seconds)'
       sxaddpar, outhdr, 'EXPTIME', 1.0, ' exposure time (seconds)'
       sxaddpar, outhdr, 'NCOMBINE', sxpar(hdr,'NCOMBINE'), ' number of stacked frames'
;      sxaddpar, outhdr, 'GAIN', sxpar(hdr,'GAIN_1')*sxpar(hdr,'NCOMBINE'), ' total gain'
       sxaddpar, outhdr, 'ZUNIT', 'count/s', ' image units'
       sxaddpar, outhdr, 'FILTER', 'Ks', ' filter bandpass'
       sxaddpar, outhdr, 'ZPT_FAT', sxpar(hdr,'ZPT_FAT'), ' magnitude zero-point'
       sxaddpar, outhdr, 'ZPTNSTAR', sxpar(hdr,'ZPTNSTAR'), ' number of stars used to determine ZPT_FAT'
       sxaddpar, outhdr, 'PHOTFLAG', 1, ' photometrically calibrated field'

       extast, hdr, astr
       putast, outhdr, astr
       sxdelpar, outhdr, 'DATE'
       sxdelpar, outhdr, 'COMMENT'
       sxdelpar, outhdr, 'HISTORY'

       splog, 'Writing '+outimagelist[ii]
       if keyword_set(unpack) then mwrfits, float(image), outimagelist[ii], outhdr, /create
       
       splog, 'Writing '+outweightlist[ii]
;      if keyword_set(unpack) then mwrfits, byte(weight gt 0), outweightlist[ii], hdr, /create
       if keyword_set(unpack) then mwrfits, float(weight), outweightlist[ii], outhdr, /create
       
    endfor

; ---------------------------------------------------------------------------    
    
;    datapath = flamingos_path()+'gonzalez.21mar2006/'
;    outpath = flamingos_path(/mar06)
;
;    imagelist = file_search(datapath+'*.fits',count=nimage)
;;   for ii = 0L, nimage-1L do spawn, 'ln -sfv '+imagelist[ii]+' '+outpath+file_basename(imagelist[ii])
;
;    outimagelist = outpath+'ss_'+strmid(file_basename(imagelist),19,/reverse)
;    outweightlist = repstr(outimagelist,'.fits','.weight.fits')
;
;    nimage = n_elements(imagelist)
;
;    for ii = 0L, nimage-1L do begin
;
;       image = mrdfits(imagelist[ii],0,hdr,/silent)
;
;       sxaddpar, hdr, 'OBJECT', 'SG1120'
;       sxaddpar, hdr, 'EXPTIME', sxpar(hdr,'EXP_TIME'), ' exposure time (seconds)'
;
;       if keyword_set(unpack) then begin
;          splog, 'Writing '+outimagelist[ii]
;          mwrfits, float(image), outimagelist[ii], hdr, /create
;
;          splog, 'Writing '+outweightlist[ii]
;          mwrfits, byte(image ne 0), outweightlist[ii], hdr, /create
;       endif 
;
;    endfor
       
;; ---------------------------------------------------------------------------    
;
;    image = mrdfits(datapath+'as_sg1120a_ks.fits',0,hdr,/silent)
;    weight = mrdfits(datapath+'expmap_sg1120a_ks.fits',0,weighthdr,/silent)
;
;    splog, 'Writing '+outpath+'sra.sg1120a_ks.fits'
;    mwrfits, image, outpath+'sra.sg1120a_ks.fits', hdr, /create
;    mwrfits, weight, outpath+'sra.sg1120a_ks.fits', weighthdr
;
;; ---------------------------------------------------------------------------    
;    
;    image = mrdfits(datapath+'as_sg1120b_ks.fits',0,hdr,/silent)
;    weight = mrdfits(datapath+'expmap_sg1120b_ks.fits',0,weighthdr,/silent)
;
;    splog, 'Writing '+outpath+'sra.sg1120b_ks.fits'
;    mwrfits, image, outpath+'sra.sg1120b_ks.fits', hdr, /create
;    mwrfits, weight, outpath+'sra.sg1120b_ks.fits', weighthdr
    
return
end
