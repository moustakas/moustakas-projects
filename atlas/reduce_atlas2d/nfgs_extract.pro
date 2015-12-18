pro nfgs_extract, wfits=wfits, noplot=noplot
; jm02jun28uofa
; jm03sep15uofa - updated
    
;   skyspecfile = 'kpno_sky.fits'
    
    nfgspath = atlas_path(/analysis1d)+'nfgs/spec2d/'
    datapath = atlas_path()
    outpath = atlas_path(/analysis1d)+'nfgs/spec1d/'
    pushd, datapath
    
; ----------------------------------------------------------------
; extract the NFGS galaxies that we observed with matched apertures
; -----------------------------------------------------------------
    
    ispec, 'ngc_5596_drift_027.fits.gz', deredden=1, /heliocor, datapath=nfgspath, aperture=52.8, $
      /meanprofile, /tracespec, outpath=outpath, noplot=noplot, wfits=wfits, /gzip ; 0.88x0.45

    ispec, 'ngc_5993_drift_027.fits.gz', deredden=1, /heliocor, datapath=nfgspath, aperture=69.6, $
      /meanprofile, /tracespec, outpath=outpath, noplot=noplot, wfits=wfits, /gzip ; 1.16x0.45

    ispec, 'ngc_6123_drift_006.fits.gz', deredden=1, /heliocor, datapath=nfgspath, aperture=66.0, $
      /meanprofile, /tracespec, outpath=outpath, noplot=noplot, wfits=wfits, /gzip ; 1.10x0.10

; ----------------------------------------------------------------
; extract galaxies that are in the atlas that are also in the NFGS
; (with mismatched apertures)
; ----------------------------------------------------------------

; I need to obtain this spectrum (122.A12300+4259) from R. Jansen
;   ispec, 'ugc_07690_drift_060.fits.gz', deredden=1, /heliocor, skyspecfile=skyspecfile, aperture=80, $
;     skyaperture=[45,50], outpath=outpath, noplot=noplot, wfits=wfits, /gzip

    ispec, 'ngc_0695_drift_030.fits.gz', deredden=1, /heliocor, refrow=59, aperture=45, $
      /meanprofile, /tracespec, sbox=10, fluxrange=[-0.1,2]*1E-15, $
      outpath=outpath, noplot=noplot, wfits=wfits, /gzip

; S/N = 6 (marginal)
    ispec, 'ugc_04787_drift_120.fits.gz', deredden=1, /heliocor, aperture=50, $
      /meanprofile, /tracespec, outpath=outpath, noplot=noplot, wfits=wfits, /gzip

; S/N = 16 (okay)
    ispec, 'ngc_3104_drift_105.fits.gz', deredden=1, /heliocor, aperture=120, $
      /meanprofile, /tracespec, outpath=outpath, noplot=noplot, wfits=wfits, /gzip

; S/N = 13 (okay)
    ispec, 'ngc_3510_drift_120.fits.gz', deredden=1, /heliocor, aperture=80, $
      /meanprofile, /tracespec, outpath=outpath, noplot=noplot, wfits=wfits, /gzip

; S/N = 11 (marginal)
    ispec, 'ngc_4248_drift_090.fits.gz', deredden=1, /heliocor, aperture=100, $
      /meanprofile, /tracespec, outpath=outpath, noplot=noplot, wfits=wfits, /gzip

; S/N = 10 (okay)
    ispec, 'ngc_4288_drift_090.fits.gz', deredden=1, /heliocor, aperture=100, $
      /meanprofile, /tracespec, outpath=outpath, noplot=noplot, wfits=wfits, /gzip

    ispec, 'ugc_07950_drift_040.fits.gz', deredden=1, /heliocor, aperture=80, $
      /meanprofile, /tracespec, outpath=outpath, noplot=noplot, wfits=wfits, /gzip

; IUE galaxy
    ispec, 'ugc_09560_drift_030.fits.gz', deredden=1, /heliocor, aperture=40, $
      /meanprofile, /tracespec, outpath=outpath, noplot=noplot, wfits=wfits, /gzip

    ispec, 'ngc_7620_drift_040.fits.gz', deredden=1, /heliocor, aperture=60, $
      /meanprofile, /tracespec, traceorder=2, outpath=outpath, noplot=noplot, wfits=wfits, /gzip

; S/N = 12 (okay)    
    ispec, 'ugc_00685_drift_040.fits.gz', deredden=1, /heliocor, aperture=100, $
      /meanprofile, /tracespec, outpath=outpath, noplot=noplot, wfits=wfits, /gzip

    popd
    
return
end
