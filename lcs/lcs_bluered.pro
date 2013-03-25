pro lcs_bluered, clobber=clobber
; jm13mar11siena - select blue and red galaxies

    lcspath = getenv('LCS_DATA')

; coefficients for selecting star-forming galaxies    
    coeff = [-2.10,0.6]

    cl = rsex(getenv('LCS_DIR')+'/lcs_sample.cat')
    ncl = n_elements(cl)
    struct_print, cl
    
    for ii = 0, ncl-1 do begin
       nsa = mrdfits(lcspath+strlowcase(cl[ii].cluster)+'_nsa.fits.gz',1)
       lir = mrdfits(lcspath+strlowcase(cl[ii].cluster)+'_lir.fits.gz',1)
       mass = alog10(nsa.mass)
       nuvmr = nsa.absmag[1]-nsa.absmag[4] ; NUV-r color

       sf = where(nuvmr lt poly(mass,coeff),nsf)

       im_mwrfits, nsa[sf], lcspath+strlowcase(cl[ii].cluster)+'_nsa_blue.fits', clobber=clobber
       im_mwrfits, lir[sf], lcspath+strlowcase(cl[ii].cluster)+'_lir_blue.fits', clobber=clobber
       
;      djs_plot, mass, nuvmr, xr=[7,12], yr=[0,8], xsty=1, ysty=1, psym=8
;      djs_oplot, mass[sf], nuvmr[sf], color='blue', psym=8

    endfor       

; combine the sample into "big" and "small" clusters
    bignsa = [mrdfits(lcspath+'coma_nsa.fits.gz',1),mrdfits(lcspath+'hercules_nsa.fits.gz',1),$
      mrdfits(lcspath+'abell1367_nsa.fits.gz',1),mrdfits(lcspath+'abell2052_nsa.fits.gz',1)]
    biglir = [mrdfits(lcspath+'coma_lir.fits.gz',1),mrdfits(lcspath+'hercules_lir.fits.gz',1),$
      mrdfits(lcspath+'abell1367_lir.fits.gz',1),mrdfits(lcspath+'abell2052_lir.fits.gz',1)]
    im_mwrfits, bignsa, 'bigclusters_nsa.fits', /clobber
    im_mwrfits, biglir, 'bigclusters_lir.fits', /clobber
    
    smallnsa = [mrdfits(lcspath+'mkw11_nsa.fits.gz',1),mrdfits(lcspath+'mkw8_nsa.fits.gz',1),$
      mrdfits(lcspath+'awm4_nsa.fits.gz',1),mrdfits(lcspath+'abell2063_nsa.fits.gz',1),$
      mrdfits(lcspath+'ngc6107_nsa.fits.gz',1)]
    smalllir = [mrdfits(lcspath+'mkw11_lir.fits.gz',1),mrdfits(lcspath+'mkw8_lir.fits.gz',1),$
      mrdfits(lcspath+'awm4_lir.fits.gz',1),mrdfits(lcspath+'abell2063_lir.fits.gz',1),$
      mrdfits(lcspath+'ngc6107_lir.fits.gz',1)]
    im_mwrfits, smallnsa, 'smallclusters_nsa.fits', /clobber
    im_mwrfits, smalllir, 'smallclusters_lir.fits', /clobber
    
return
end
    
