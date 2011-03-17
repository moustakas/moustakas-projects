pro hst_split_gonzalez_drizzle
; jm07aug15nyu - Anthony's multidrizzled images have the data in
;                the first extension and the weight map in the second
;                extension; split these into individual FITS files for
;                use with HST_MOSAICS

    outpath = hst_path()+'gonzalez_drizzle/'
    datapath = outpath+'MEF/'

    flist = [file_search(datapath+'p?.fits'),file_search(datapath+'p10.fits')]
    fcount = n_elements(flist)
    pos = string(lindgen(fcount)+1L,format='(I2.2)')
    outlist = outpath+'hst_f814w_drz_pos'+pos+'.fits'
    
;   outlist = outpath+repstr(file_basename(flist),'p','hst_f814w_drz_pos')
    outweight = repstr(outlist,'.fits','.weight.fits')

    for i = 0L, fcount-1L do begin
       splog, 'Reading '+flist[i]
       hdr = headfits(flist[i],ext=0L)

       sci = mrdfits(flist[i],1,scihdr,/silent)
       sxaddpar, scihdr, 'OBJECT', 'SG1120_POS'+pos[i]
       splog, 'Writing '+outlist[i]
       writefits, outlist[i], 0, hdr
       writefits, outlist[i], temporary(sci), scihdr, /append

       wgt = mrdfits(flist[i],2,wgthdr,/silent)
       sxaddpar, wgthdr, 'OBJECT', 'SG1120_POS'+pos[i]
       splog, 'Writing '+outweight[i]
       writefits, outweight[i], 0, hdr
       writefits, outweight[i], temporary(wgt), wgthdr, /append

       print
    endfor
    
return
end
    
