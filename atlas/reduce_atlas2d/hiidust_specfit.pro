pro hiidust_specfit
; jm06mar27uofa
    
    datapath = '/home/ioannis/06mar/extracted/'
    specfitpath = datapath

    suffix = 'hiidust'

    speclist = file_basename(file_search(datapath+'*.ms.fits*'))

    specdata = ispeclinefit(speclist,specres=8.0,snrcut=1.0,dustmodel=0,$
      datapath=datapath,linepath=specfitpath,suffix=suffix,Zmulti=Zmulti,$
      /charlot,zcrosscor=1,/postscript,/write,vmaxshift=vmaxshift,/nologfile,$
      starvdisp=100.0)

    snrcut_linedust = 3.0
    snrcut_abundance = 3.0

    line = parse_ispeclinefit(datapath=datapath,root=root,$
      select_lines=select_lines,syserr=syserr,disttag=disttag,$
      disterrtag=disterrtag,photerrtag=photerrtag,outfile=outfile,$
      /electrondensity,snrcut_linedust=snrcut_linedust,$
      snrcut_abundance=snrcut_abundance,/match,/kauffmann,/write,$
      /odonnell,/nopropagate,_extra=extra)

return
end
    
