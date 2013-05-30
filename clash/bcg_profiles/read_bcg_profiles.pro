function read_bcg_profiles, cluster
; jm13may20siena - read a M. Postman style BCG profile file 

    profilepath = getenv('CLASH_DATA')+'/bcg_profiles/13may20/'
    if n_elements(cluster) eq 0 then begin
       splog, 'Need CLUSTER name.'
       return, -1
    endif

    clash = rsex(getenv('CLASH_DIR')+'/clash_sample.sex')
    filt = clash_filterlist(short=short,weff=weff,zpt=zpt)
    kl = k_lambda(weff,/odon)
    nfilt = n_elements(filt)

    this =  where(strtrim(cluster,2) eq strtrim(clash.cluster_short,2))
    if this[0] eq -1 then begin
       splog, 'Unknown cluster '+cluster
       return, -1
    endif
    arcsec2kpc = dangular(clash[this].z,/kpc)/206265D ; [kpc/arcsec]
    ebv = clash[this].ebv
    
    nrmax = 25
    data_template = {$
      sma:     fltarr(nrmax)-99.0,$
      mu:      fltarr(nrmax)-99.0,$
      mu_err:  fltarr(nrmax)-99.0,$
      ell:     fltarr(nrmax)-99.0,$
      ell_err: fltarr(nrmax)-99.0,$
      pa:      fltarr(nrmax)-99.0,$
      pa_err:  fltarr(nrmax)-99.0}
    allfile = profilepath+cluster+'_mosaic_*_*_'+short+'_*.txt'

    nrow = 9
    check = fltarr(nfilt)+1
    for ii = 0, nfilt-1 do begin
       data1 = data_template
       file1 = file_search(allfile[ii],count=nfile)
       if nfile eq 0 then check[ii] = 0 else begin
          txt = djs_readilines(allfile[ii],nhead=4,head=head)
          scale = double(strmid(head[2],strpos(head[2],'Scale = ')+8)) ; [arcsec/pixel]

          niso = n_elements(txt)
          newtxt = reform(strsplit(strjoin(repstr(repstr(repstr(txt,' ',';'),')',''),$
            '(',''),';'),';',/extract),nrow,niso)
          data1.ell[0:niso-1] = reform(newtxt[3,*])
          data1.ell_err[0:niso-1] = reform(newtxt[4,*])
          data1.pa[0:niso-1] = reform(newtxt[5,*])
          data1.pa_err[0:niso-1] = reform(newtxt[6,*])

; semi-major axis          
          data1.sma[0:niso-1] = reform(newtxt[0,*])*arcsec2kpc ; [kpc]
          
; surface brightness
          mu = reform(newtxt[1,*])
          mu_err = reform(newtxt[2,*])
          data1.mu[0:niso-1] = -2.5*alog10(mu)+5.0*alog10(scale)+$ ; [ABmag/arcsec^2]
            zpt[ii]-kl[ii]*ebv
          data1.mu_err[0:niso-1] = 2.5*mu_err/mu/alog(10)
          
;         splog, short[ii], niso, minmax(data1.sma)
       endelse
       if ii eq 0 then data = data1 else data = [data,data1]
    endfor

; was this cluster fitted?    
    if total(check) eq 0.0 then begin
       splog, 'No SB profiles fitted for cluster '+cluster
       return, -1
    endif
    
return, data
end
