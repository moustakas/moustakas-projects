pro kenn92_parse_specfit, line, datapath=datapath, _extra=extra
; jm04apr28uofa

    if n_elements(datapath) eq 0L then datapath = kenn92_path(/specfit)
    analysis_path = kenn92_path(/analysis)
    
    root = 'kenn92'
    outfile = root+'_speclinefit.fits'

; read the ancillary data and prepend

    infofile = 'kenn92_info.fits.gz'
    if (file_test(analysis_path+infofile,/regular) eq 0L) then begin
       splog, 'Unable to find '+analysis_path+infofile+'.'
       return
    endif

    splog, 'Reading '+analysis_path+infofile+'.'
    kenn92 = mrdfits(analysis_path+infofile,1,/silent)
    kenn92 = struct_trimtags(kenn92,except=except)

; parse and write out everything
    
    select_lines = ['OII_3727','H_BETA','OIII_5007','H_ALPHA']
    disttag = 'DISTANCE'
    disterrtag = 'DISTANCE_ERR'

    snrcut_linedust = 1.0
    snrcut_abundance = 1.0

    line = parse_ispeclinefit(datapath=datapath,prepend=kenn92,root=root,$
      trimtags='GALAXY',select_lines=select_lines,$
      disttag=disttag,disterrtag=disterrtag,photerrtag=photerrtag,$
      outfile=outfile,/electrondensity,snrcut_linedust=snrcut_linedust,$
      snrcut_abundance=snrcut_abundance,/match,/kauffmann,/write,$
      /odonnell,/nopropagate,_extra=extra)

return
end
