pro ndwfs_sdss_stars
; jm09aug25ucsd - grab the stars in the NDWFS field

    indir = getenv('RESEARCHPATH')+'/data/ndwfs/'

    outfile = indir+'NDWFS_SDSS_stars.fits.gz'
    sdss = sdss_sweep_circle(218.02,34.28,3.0,type='star')
    im_mwrfits, sdss, outfile

;   djs_plot, sdss1.ra, sdss1.dec, ps=3, color='red', xsty=3, ysty=3
;   djs_oplot, ndwfs1.r_alpha_j2000, ndwfs1.r_delta_j2000, ps=3                

stop    
    
return    
end
