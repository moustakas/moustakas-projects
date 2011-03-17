function read_ediscs_specfit, specfile, cluster, silent=silent
; jm07apr16nyu - read the best-fitting EDisCS spectra; sort by cluster 

    nspecfile = n_elements(specfile)
    if (nspecfile eq 0L) then begin
       print, 'specfit = read_ediscs_specfit(specfile,_extra=extra)'
       return, -1L
    endif

    datapath = ediscs_path(/specfit)

    specfile = strcompress(specfile,/remove)
    cluster = strcompress(cluster,/remove)
    ucluster = cluster[uniq(cluster,sort(cluster))]
    ncluster = n_elements(ucluster)

; sort through the various SPECFILES and find the most recent one

    specfit = fltarr(2201,5,nspecfile)
    
    for icluster = 0L, ncluster-1L do begin

       specdatafiles = file_basename(file_search(datapath+'?????_ediscs_'+ucluster[icluster]+'_specdata.fits.gz',count=fcount))

       if (fcount eq 0L) then message, 'Problem here!'
       specdatafile = specdatafiles[(reverse(sort(specdatafiles)))[0]]
       specfitfile = repstr(specdatafile,'specdata','specfit')

       specdata = mrdfits(datapath+specdatafile,1,/silent)

       these = where((ucluster[icluster] eq cluster),nthese)
       for ii = 0L, nthese-1L do begin

          if (not keyword_set(silent)) then print, format='("Reading cluster ",I0,"/",I0," and specfile ",I0,"/",I0,".    ",A1,$)', $
            icluster+1, ncluster, ii+1, nthese, string(13b)

          indx = where(strtrim(specfile[these[ii]],2) eq strtrim(specdata.specfile,2),nindx)
          if (nindx eq 0L) then message, 'Problem here!'

          ext_no = indx[0]+2
          specfit1 = mrdfits(datapath+specfitfile,ext_no,/silent)
          specfit[0L:n_elements(reform(specfit1[*,0]))-1L,*,these[ii]] = specfit1

       endfor

    endfor
    if (not keyword_set(silent)) then print
    
return, specfit
end
