; RELEGATED! - jm08mar25nyu
stop
pro ediscs_parse_specfit, ancillary, _extra=extra
; jm06sep27nyu
    
    version = ediscs_specfit_version()

    datapath = ediscs_path(/specfit)
    analysis_path = ediscs_path(/analysis)

    ancillary = ediscs_read_ancillary()
;   ancillary = ediscs_read_info()
    
; parse and write out the spectral fitting results

    select_lines = ['OII_3727','H_BETA','OIII_5007','H_ALPHA']
    snrcut_linedust = 1.0
    snrcut_abundance = 1.0

    root = 'ediscs'
    outfile = root+'_speclinefit_'+version+'.fits'
 
    fitspecfile = file_basename(file_search(datapath+'ediscs_specdata_'+version+'.fits.gz',count=fcount))
;   fitspecfile = file_basename(file_search(datapath+'?????_ediscs_specdata.fits.gz',count=fcount))
    if (fcount ne 1L) then begin
       splog, 'FITSPECFILE not found!'
       return
    endif else fitspecfile = fitspecfile[0]

; PREPEND supersedes ANCILLARY    

    t0 = systime(1)
    line = parse_ispeclinefit(fitspecfile=fitspecfile,datapath=datapath,$
      prepend=ancillary,trimtags=['SPECFILE','GALAXY'],$$ ;ancillary=ancillary,$
      root=root,select_lines=select_lines,syserr=syserr,disttag=disttag,disterrtag=disterrtag,$
      photerrtag=photerrtag,outfile=outfile,/electrondensity,snrcut_linedust=snrcut_linedust,$
      snrcut_abundance=snrcut_abundance,/kauffmann,/write,/odonnell,/nopropagate,$
      /noabundance,_extra=extra)
    splog, 'Total time = '+strtrim(string((systime(1)-t0)/60.0,format='(F12.1)'),2)+' minutes.'

return
end
