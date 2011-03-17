;+
; NAME:
;       SINGS_PARSE_SPECFIT
;
; PURPOSE:
;       Parse the output from fitting the SINGS galaxies.
;
; INPUTS:
;
; OPTIONAL INPUTS:
;
; KEYWORD PARAMETERS:
;       nuclear - parse the nuclear spectral fitting
;       drift20 - parse the 20" scan spectral fitting
;       drift56 - parse the 56" scan spectral fitting
;
; OUTPUTS:
;
; OPTIONAL OUTPUTS:
;       line - output from PARSE_ISPECLINEFIT()
;
; COMMENTS:
;
; EXAMPLES:
;
; MODIFICATION HISTORY:
;       J. Moustakas, 2005 Jul 26, U of A
;       jm08oct20nyu - totally rewritten - do not use
;         PARSE_ISPECLINEFIT() 
;-

pro sings_parse_specfit, nuclear=nuclear, drift20=drift20, drift56=drift56

    version = sings_version(/specfit)
    outpath = sings_path(/specfit)
    specfitpath = outpath+version+'/'

    if (not keyword_set(nuclear)) and (not keyword_set(drift20)) and $
      (not keyword_set(drift56)) then begin
       splog, 'Either NUCLEAR *or* DRIFT20 *or* DRIFT56 keyword must be set!'
       return
    endif

    if (keyword_set(nuclear) and keyword_set(drift20)) or $
      (keyword_set(nuclear) and keyword_set(drift56)) or $
      (keyword_set(drift20) and keyword_set(drift56)) then begin
       splog, 'Only one keyword (NUCLEAR, DRIFT20, or DRIFT56) can be set at the same time!'
       return
    endif

    if keyword_set(nuclear) then begin
       root = 'sings_nuclear'
       except = ['DRIFT20*','DRIFT56*']
    endif

    if keyword_set(drift20) then begin
       root = 'sings_drift20'
       except = ['NUCLEAR*','DRIFT56*']
    endif

    if keyword_set(drift56) then begin
       root = 'sings_drift56'
       except = ['NUCLEAR*','DRIFT20*']
    endif

; read the SINGS info file

    sings = sings_read_info()

; find and read the current SPECDATA file

    allfiles = file_search(specfitpath+'*'+root+'*_specdata.fits.gz')
    specdatafile = allfiles[(reverse(sort(allfiles)))[0]]
    splog, 'Reading '+specdatafile
    specdata = mrdfits(specdatafile,1,/silent)

    match, strtrim(specdata.galaxy,2), strtrim(sings.galaxy,2), m1, m2
    
    out = struct_addtags(struct_trimtags(sings[m2],except=except),$
      struct_trimtags(specdata[m1],except=['galaxy']))
    out = out[sort(strtrim(out.galaxy,2))]

; write out    
    
    outfile = outpath+root+'_speclinefit_'+version+'.fits'
    splog, 'Writing '+outfile
    mwrfits, out, outfile, /create
    spawn, 'gzip -f '+outfile, /sh

; old/obsolete code    
;   select_lines = ['OII_3727','H_BETA','OIII_5007','H_ALPHA']
;   snrcut_linedust = 1.0
;   snrcut_abundance = 1.0
;
;   line = parse_ispeclinefit(datapath=datapath,prepend=sings,root=root,$
;     trimtags='GALAXY',select_lines=select_lines,$
;     outfile=outfile,snrcut_linedust=snrcut_linedust,$
;     snrcut_abundance=snrcut_abundance,/match,/kauffmann,/write,$
;     /odonnell,/nopropagate,/noabundance,_extra=extra)

return
end
