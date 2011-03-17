pro sc1120_parse_specfit, line, _extra=extra
; jm05jan28uofa

; read the ancillary data
    
    datapath = sc1120_path(/analysis)

    root = 'sc1120'
    outfile = root+'_speclinefit.fits'

    pushd, datapath
    fitspecfile = (file_search('?????_'+root+'_specdata.fits*'))[0]
    popd

    data = mrdfits(datapath+root+'_data.fits.gz',1,/silent)

; read and append the synthesized magnitudes and K-corrections 
;
;   if file_test(datapath+root+'_kcorr_merged.fits.gz',/regular) then begin
;
;      splog, 'Reading and appending '+datapath+root+'_kcorr_merged.fits.gz.'
;      kcorr = mrdfits(datapath+root+'_kcorr_merged.fits.gz',1,/silent)
;      
;      match, data.galaxy, kcorr.galaxy, match1, match2
;      data = struct_addtags(data[match1],struct_trimtags(kcorr[match2],except=['GALAXY']))
;
;   endif

    select_lines = ['OII_3727','H_BETA','OIII_5007','H_ALPHA']

    snrcut_linedust = 5.0
    snrcut_abundance = 5.0

; parse and write out

    line = parse_ispeclinefit(fitspecfile=fitspecfile,datapath=datapath,$
      prepend=data,root=root,/nopropagate,trimtags=['GALAXY','SPECFILE'],$
      /match,select_lines=select_lines,snrcut_abundance=snrcut_abundance,$
      snrcut_linedust=snrcut_linedust,/kauffmann,outfile=outfile,$
      /write,_extra=extra)
    
stop
    
return
end
