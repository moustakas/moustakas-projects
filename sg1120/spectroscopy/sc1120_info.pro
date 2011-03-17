pro sc1120_info, wfits=wfits
; jm05jan28uofa
; merge the redshift catalog and other information from the headers
; into a FITS data structure; this routine should be run after
; SC1120_WRITE_SIGSPEC 
    
    red, h100=0.7, omega0=0.3, omega_lambda=0.7

    datapath = sc1120_path(/sigspec)
    outpath = sc1120_path(/analysis)

; forage the headers in DATAPATH

    splog, 'Foraging headers in '+datapath
    pushd, datapath & speclist = file_search('*.fits',count=nspec) & popd
    forage = sc1120_forage(speclist,datapath=datapath)

    tags = [$
      'GALAXY',  $
      'SPECFILE',$
      'MASK',    $
      'QUADRANT',$
      'SLIT',    $
      'RA',      $
      'DEC',     $
      'AIRMASS', $
      'EXPTIME', $
      'Z',       $
      'Z_ERR',   $
      'Q',       $
      'YOBJ',    $
      'YSKY']

    info = struct_trimtags(forage,select=tags)

; add tags

    addinfo = {$
      distance:             -999D, $
      distance_err:         -999D}
    addinfo = replicate(addinfo,nspec)

    addinfo.distance     = dluminosity(info.z,/Mpc)
    addinfo.distance_err = 0.0

    info = struct_addtags(info,addinfo)
    
    if keyword_set(wfits) then begin
       infofile = 'sc1120_data.fits'
       splog, 'Writing '+outpath+infofile+'.'
       mwrfits, info, outpath+infofile, /create
       spawn, ['gzip -f '+outpath+infofile], /sh
    endif

return
end
    
