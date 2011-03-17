pro sc1120_specfit
; jm05jan28uofa

    analysis_path = sc1120_path(/analysis)
    datapath = sc1120_path(/spec1d)

    pushd, datapath
    speclist = file_search('*.fits',count=fcount)
    popd

    suffix = 'sc1120'

; an instrumental resolution of 10.6 and a 100 km/s velocity
; dispersion gives a spectral resolution of ~12 Angstroms at 7080
; Angstroms 
    
    specdata = ispeclinefit(speclist,specres=10.0,snrcut=0.0,dustmodel=0,$
      datapath=datapath,linepath=analysis_path,suffix=suffix,/charlot,/nolog,$
      Zmulti=0,/zcrosscor,/postscript,/write,starvdisp=50.0,vmaxshift=500.0)
    
    sc1120_parse_specfit

stop
    
return
end    
