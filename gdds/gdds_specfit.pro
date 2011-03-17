pro gdds_specfit
; jm05feb08uofa
    
    datapath = gdds_path(/spec1d)

    pushd, datapath
    speclist = file_search('*.fits.gz')
    popd

    suffix = 'gdds'

; try a lower spectral resolution
    
    specdata = ispeclinefit(speclist,specres=16.0,snrcut=0.0,dustmodel=0,$
      datapath=datapath,linepath=datapath,suffix=suffix,/charlot,/zcrosscor,$
      /postscript,/write,/nologfile,starvdisp=0.0)

stop    
    
    outfile = 'gdds_speclinefit.fits'
    line = parse_ispeclinefit(root='gdds',outfile=outfile,/write,/hbeta,snrcut=1.0)
    
return
end
    

