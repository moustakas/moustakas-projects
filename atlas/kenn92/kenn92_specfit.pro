pro kenn92_specfit, _extra=extra
; jm04jan04uofa
; jm05aug03uofa - updated

    datapath = kenn92_path(/spec1d)
    specfitpath = kenn92_path(/specfit)
    kenn92 = kenn92_read_info()

;   suffix = 'test'
    suffix = 'kenn92'

    keep = where((kenn92.drift_lowresflag eq 0L) and (kenn92.drift_agnflag eq 0L))
    speclist = strtrim(kenn92[keep].drift_file,2)

    specdata = ispeclinefit(speclist,specres=6.0,snrcut=0.0,dustmodel=0,$
      datapath=datapath,linepath=specfitpath,suffix=suffix,Zmulti=Zmulti,$
      /charlot,/zcrosscor,/postscript,/write,vmaxshift=1000.0,/nologfile,$
      starvdisp=100.0,nback=4)

    kenn92_parse_specfit, line, datapath=specfitpath

stop
    
return
end    
