function ages_allpasses, fluxed=fluxed
; jm09nov13ucsd - return all the passes in AGES

    spec1dpath = ages_path(/spec1d)+$
      'fluxed/before_skysubpca/'
    allfiles = file_search(spec1dpath+'spectra_???.fits.gz')
    allpass = strmid(file_basename(allfiles),8,3)

    if keyword_set(fluxed) then begin
       keep = where((allpass ne '106') and $
         (allpass ne '110') and (allpass ne '209') and $
         (allpass ne '310') and (allpass ne '311'))
       allpass = allpass[keep]
    endif
    
return, allpass
end
    
