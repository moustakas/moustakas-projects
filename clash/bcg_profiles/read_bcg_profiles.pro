function read_bcg_profiles, cluster, scale=scale, these_filters=these_filters
; jm13may20siena - read a M. Postman style BCG profile file 

;   profilepath = getenv('CLASH_DATA')+'/bcg_profiles/13oct29/'
    profilepath = getenv('CLASH_DATA')+'/bcg_profiles/14jul12/'
    if n_elements(cluster) eq 0 then begin
       splog, 'Need CLUSTER name.'
       return, -1
    endif

    clash = rsex(getenv('CLASH_DIR')+'/clash_sample.sex')

    if n_elements(scale) eq 0 then scale = '065mas' ; other option is '030mas'
    
; read the filters, optionally restricting the list
    filt = clash_filterlist(short=short,weff=weff,zpt=zpt)
    if n_elements(these_filters) ne 0 then begin
       match, short, these_filters, m1, m2
;      match, filt, these_filters, m1, m2
       if n_elements(m2) ne n_elements(these_filters) or m2[0] eq -1 then $
         message, 'Problem matching filters!'
       srt = sort(m2)
       filt = filt[m1[srt]]
       short = short[m1[srt]]
       weff = weff[m1[srt]]
       zpt = zpt[m1[srt]]
    endif
    
    kl = k_lambda(weff,/odon)
    nfilt = n_elements(filt)

    this =  where(strtrim(cluster,2) eq strtrim(clash.shortname,2))
    if this[0] eq -1 then begin
       splog, 'Unknown cluster '+cluster
       return, -1
    endif
    arcsec2kpc = dangular(clash[this].z,/kpc)/206265D ; [kpc/arcsec]
    ebv = clash[this].ebv
    
    nrmax = 25
    data_template = {$
      band:                        '',$
      sma:         fltarr(nrmax)-99.0,$
      mu:          fltarr(nrmax)-99.0,$
      mu_err:      fltarr(nrmax)-99.0,$
      abmag:       fltarr(nrmax)-99.0,$
      abmag_err:   fltarr(nrmax)-99.0,$
      maggies:     fltarr(nrmax),$
      ivarmaggies: fltarr(nrmax),$
      ell:         fltarr(nrmax)-99.0,$
      ell_err:     fltarr(nrmax)-99.0,$
      pa:          fltarr(nrmax)-99.0,$
      pa_err:      fltarr(nrmax)-99.0}
    allfile = profilepath+cluster+'_'+scale+'_'+short+'_drz_????????_bcg.txt'
;   allfile = profilepath+cluster+'_mosaic_*'+scale+'*_*_'+short+'_*.txt'
;   splog, allfile

    nrow = 9
    check = fltarr(nfilt)+1
    for ii = 0, nfilt-1 do begin
       data1 = data_template
       data1.band = short[ii]
       file1 = file_search(allfile[ii],count=nfile)
; multiple files
       if nfile gt 1 then begin
          splog, 'Warning: multiple profile files!'
          niceprint, '   '+file_basename(file1)
          date = strmid(file1,15,8,/reverse)
          this = max(date,indx1)
          file1 = file1[indx1]
          nfile = 1
       endif
       
       if nfile eq 0 then check[ii] = 0 else begin
          txt = djs_readilines(file1,nhead=4,head=head)
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
          counts = reform(newtxt[1,*]) ; flux [electron/s/pixel]
          counts_err = reform(newtxt[2,*])
          
; surface brightness in [ABmag/arcsec^2]
          data1.mu[0:niso-1] = -2.5*alog10(counts)+5.0*alog10(scale)+zpt[ii]-kl[ii]*ebv
          data1.mu_err[0:niso-1] = 2.5*counts_err/counts/alog(10)
          
;         fact = 10.0^(-0.4*(zpt[ii]-kl[ii]*ebv))
;         data1.maggies[0:niso-1] = counts*fact
;         data1.ivarmaggies[0:niso-1] = 1.0/(counts_err*fact)^2
;         splog, short[ii], niso, minmax(data1.sma)
       endelse
       if ii eq 0 then data = data1 else data = [data,data1]
    endfor
    
;; now go back through and compute the apparent magnitude (and maggies)
;; at each position along the semimajor axis by multiplying the SB by
;; the area of each isophote; this is a hack because the position
;; angles and ellipticities of the isophotes change with semi-major
;; axis, but just assume uniform (median) ellipse parameters
;    f160w = where(short eq 'f160w')
;    good = where(data[f160w].pa gt -90 and data[f160w].ell gt -90,ngood)
;    medpa = djs_median(data[f160w].pa[good])
;    medell = djs_median(data[f160w].ell[good])
;
;    sma = data[f160w].sma[good]
;    area = fltarr(ngood)
;    for rr = 0, ngood-1 do begin
;       case rr of
;          0: begin
;             rin = 0
;             rout = data[f160w].sma[good[rr]]+(data[f160w].sma[good[rr+1]]-data[f160w].sma[good[rr]])/2
;          end
;          ngood-1: begin
;             diff = (data[f160w].sma[good[rr]]-data[f160w].sma[good[rr-1]])/2
;             rin = data[f160w].sma[good[rr]]-diff
;             rout = data[f160w].sma[good[rr]]+diff
;          end
;          else: begin
;             rin = data[f160w].sma[good[rr]]-(data[f160w].sma[good[rr]]-data[f160w].sma[good[rr-1]])/2
;             rout = data[f160w].sma[good[rr]]+(data[f160w].sma[good[rr+1]]-data[f160w].sma[good[rr]])/2
;          end
;       endcase
;;      splog, rin, data[f160w].sma[good[rr]], rout
;
;       area[rr] = !pi*medell*(rout^2-rin^2)/arcsec2kpc^2 ; [arcsec^2]
;;      area[rr] = !pi*data[f160w].ell[good[rr]]*(rout^2-rin^2)/arcsec2kpc^2 ; [arcsec^2]
;    endfor
;
;; now get ABmag and maggies in each band by scaling by the area
;    for ii = 0, nfilt-1 do begin
;       good1 = where(data[ii].mu[good] gt -90,ngood1)
;       if ngood1 ne 0 then begin
;          data[ii].abmag[good[good1]] = data[ii].mu[good[good1]] - 5*alog10(area[good1])
;          data[ii].abmag_err[good[good1]] = data[ii].mu[good[good1]]
;          data[ii].maggies[good[good1]] = mag2maggies(data[ii].abmag[good[good1]],$
;            magerr=data[ii].abmag_err[good[good1]],ivarmaggies=ivarmaggies)
;          data[ii].ivarmaggies[good[good1]] = ivarmaggies
;       endif   
;;      niceprint, area, data[ii].abmag[good] & print
;    endfor

; was this cluster fitted?    
    if total(check) eq 0.0 then begin
       splog, 'No SB profiles fitted for cluster '+cluster
       return, -1
    endif
    
return, data
end
