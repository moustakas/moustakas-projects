pro templates_extract, wfits=wfits, noplot=noplot
; jm02jun25uofa

    datapath = '/home/ioannis/kennicutt/templates/spec2d/'
    outpath = '/home/ioannis/kennicutt/templates/spec1d/'

    pushd, datapath

    ispec, 'bd+023375.fits.gz', deredden=1, /heliocor, /tracespec, traceorder=2, aperture=30, skyaperture=60, outpath=outpath, wfits=wfits, noplot=noplot, /gzip
    ispec, 'bd+302347.fits.gz', deredden=1, /heliocor, /tracespec, traceorder=2, aperture=30, skyaperture=60, outpath=outpath, wfits=0, noplot=noplot, /gzip ; no radial velocity!
    ispec, 'bd+302611.fits.gz', deredden=1, /heliocor, /tracespec, traceorder=2, aperture=30, skyaperture=60, outpath=outpath, wfits=wfits, noplot=noplot, /gzip
    ispec, 'hd102224.fits.gz',  deredden=1, /heliocor, /tracespec, traceorder=2, aperture=30, skyaperture=60, outpath=outpath, wfits=wfits, noplot=noplot, /gzip
    ispec, 'hd104979.fits.gz',  deredden=1, /heliocor, /tracespec, traceorder=2, aperture=30, skyaperture=60, outpath=outpath, wfits=wfits, noplot=noplot, /gzip
    ispec, 'hd108317.fits.gz',  deredden=1, /heliocor, /tracespec, traceorder=2, aperture=30, skyaperture=60, outpath=outpath, wfits=wfits, noplot=noplot, /gzip
    ispec, 'hd110964.fits.gz',  deredden=1, /heliocor, /tracespec, traceorder=2, aperture=30, skyaperture=60, outpath=outpath, wfits=0, noplot=noplot, /gzip ; no radial velocity!
    ispec, 'hd126868.fits.gz',  deredden=1, /heliocor, /tracespec, traceorder=2, aperture=30, skyaperture=60, outpath=outpath, wfits=wfits, noplot=noplot, /gzip
    ispec, 'hd132052.fits.gz',  deredden=1, /heliocor, /tracespec, traceorder=2, aperture=30, skyaperture=60, outpath=outpath, wfits=wfits, noplot=noplot, /gzip
    ispec, 'hd147470.fits.gz',  deredden=1, /heliocor, /tracespec, traceorder=2, aperture=30, skyaperture=60, outpath=outpath, wfits=wfits, noplot=noplot, /gzip
    ispec, 'hd147550.fits.gz',  deredden=1, /heliocor, /tracespec, traceorder=2, aperture=30, skyaperture=60, outpath=outpath, wfits=wfits, noplot=noplot, /gzip
    ispec, 'hd150177.fits.gz',  deredden=1, /heliocor, /tracespec, traceorder=2, aperture=30, skyaperture=60, outpath=outpath, wfits=wfits, noplot=noplot, /gzip
    ispec, 'hd157089.fits.gz',  deredden=1, /heliocor, /tracespec, traceorder=2, aperture=30, skyaperture=60, outpath=outpath, wfits=wfits, noplot=noplot, /gzip
    ispec, 'hd158659.fits.gz',  deredden=1, /heliocor, /tracespec, traceorder=2, aperture=30, skyaperture=60, outpath=outpath, wfits=wfits, noplot=noplot, /gzip
    ispec, 'hd161677.fits.gz',  deredden=1, /heliocor, /tracespec, traceorder=2, aperture=30, skyaperture=60, outpath=outpath, wfits=wfits, noplot=noplot, /gzip
    ispec, 'hd162028.fits.gz',  deredden=1, /heliocor, /tracespec, traceorder=2, aperture=30, skyaperture=60, outpath=outpath, wfits=wfits, noplot=noplot, /gzip
    ispec, 'hd162652.fits.gz',  deredden=1, /heliocor, /tracespec, traceorder=2, aperture=30, skyaperture=60, outpath=outpath, wfits=wfits, noplot=noplot, /gzip
    ispec, 'hd162736.fits.gz',  deredden=1, /heliocor, /tracespec, traceorder=2, aperture=30, skyaperture=60, outpath=outpath, wfits=wfits, noplot=noplot, /gzip
    ispec, 'hd163346.fits.gz',  deredden=1, /heliocor, /tracespec, traceorder=2, aperture=30, skyaperture=60, outpath=outpath, wfits=wfits, noplot=noplot, /gzip
    ispec, 'hd163792.fits.gz',  deredden=1, /heliocor, /tracespec, traceorder=2, aperture=30, skyaperture=60, outpath=outpath, wfits=wfits, noplot=noplot, /gzip
    ispec, 'hd164257.fits.gz',  deredden=1, /heliocor, /tracespec, traceorder=2, aperture=30, skyaperture=60, outpath=outpath, wfits=wfits, noplot=noplot, /gzip
    ispec, 'hd164557.fits.gz',  deredden=1, /heliocor, /tracespec, traceorder=2, aperture=30, skyaperture=60, outpath=outpath, wfits=wfits, noplot=noplot, /gzip
    ispec, 'hd166384.fits.gz',  deredden=1, /heliocor, /tracespec, traceorder=2, aperture=30, skyaperture=60, outpath=outpath, wfits=wfits, noplot=noplot, /gzip
    ispec, 'hd167946.fits.gz',  deredden=1, /heliocor, /tracespec, traceorder=2, aperture=30, skyaperture=60, outpath=outpath, wfits=wfits, noplot=noplot, /gzip
    ispec, 'hd173158.fits.gz',  deredden=1, /heliocor, /tracespec, traceorder=2, aperture=30, skyaperture=60, outpath=outpath, wfits=wfits, noplot=noplot, /gzip
    ispec, 'hd175305.fits.gz',  deredden=1, /heliocor, /tracespec, traceorder=2, aperture=30, skyaperture=60, outpath=outpath, wfits=wfits, noplot=noplot, /gzip
    ispec, 'hd184406.fits.gz',  deredden=1, /heliocor, /tracespec, traceorder=2, aperture=30, skyaperture=60, outpath=outpath, wfits=wfits, noplot=noplot, /gzip
    ispec, 'hd186293.fits.gz',  deredden=1, /heliocor, /tracespec, traceorder=2, aperture=30, skyaperture=60, outpath=outpath, wfits=wfits, noplot=noplot, /gzip
    ispec, 'hd192832.fits.gz',  deredden=1, /heliocor, /tracespec, traceorder=2, aperture=30, skyaperture=60, outpath=outpath, wfits=wfits, noplot=noplot, /gzip
    ispec, 'hd205512.fits.gz',  deredden=1, /heliocor, /tracespec, traceorder=2, aperture=30, skyaperture=60, outpath=outpath, wfits=wfits, noplot=noplot, /gzip
    ispec, 'hd338529.fits.gz',  deredden=1, /heliocor, /tracespec, traceorder=2, aperture=30, skyaperture=60, outpath=outpath, wfits=wfits, noplot=noplot, /gzip
    ispec, 'sao20899.fits.gz',  deredden=1, /heliocor, /tracespec, traceorder=2, aperture=30, skyaperture=60, outpath=outpath, wfits=wfits, noplot=noplot, /gzip
    ispec, 'sao63349.fits.gz',  deredden=1, /heliocor, /tracespec, traceorder=2, aperture=30, skyaperture=60, outpath=outpath, wfits=wfits, noplot=noplot, /gzip
    
    popd

return
end
