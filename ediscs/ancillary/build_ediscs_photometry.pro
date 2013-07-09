function read_completeness
; parse Bianca's completeness files; see the emails from Bianca
; in PATH for details 

    path = ediscs_path(/catalogs)+'completeness/'
    mycatpath = ediscs_path(/mycatalogs)
    alltab = file_search(path+'*.tab',count=ntab)

    result_template = {galaxy: '', ra: 0.0D, dec: 0.0D, $
      z: 0.0, zflag: 0, spec_weight: 0.0, geom_weight: 0.0}
    for itab = 0L, ntab-1L do begin
;      splog, 'Reading '+alltab[itab]
       readcol, alltab[itab], idnew, wmag, wgeom, /silent, $
         ra, dec, z, format='A,F,F,F,F,A', comment='#'
       nobj = n_elements(idnew)
       result1 = replicate(result_template,nobj)

       result1.galaxy = idnew
       result1.ra = 15.0D*ra
       result1.dec = dec
       result1.zflag = strmatch(z,'*:*') ; 1=uncertain redshift
       result1.z = float(repstr(z,':','')) 
       result1.spec_weight = wmag
       result1.geom_weight = wgeom
       if (n_elements(result) eq 0L) then result = result1 else $
         result = [result1,result]
    endfor

return, result
end

pro build_ediscs_photometry, outphot, clobber=clobber
; jm09jun17ucsd - build the EDisCS photometry structure

    common ediscs_phot, allout, ircat
    
    vv = 'v23'
    catpath = ediscs_path(/catalogs)
    mycatpath = ediscs_path(/mycatalogs)
    
    obsfilters = [['FORS_B','FORS_V','FORS2_R','FORS_I']+$
      '_ccd',['SOFI_J','SOFI_Ks']]+'_atm.par'
    vega2ab = k_vega2ab(filterlist=obsfilters,/kurucz)
    weff = k_lambda_eff(filterlist=obsfilters)
    kl = k_lambda(weff,/odonnell)

; read the aperture corrections table (which, note, is line-matched to
; the photometry table)
    apercorpath = catpath+'apercor/'
    apercorfile = apercorpath+'kronapercorr.v2.3.fits.gz'
    splog, 'Reading '+apercorfile
    allapercor = mrdfits(apercorfile,1)

; parse the photometry       
    bands = ['B','V','R','I','J','K']
    nband = n_elements(bands)

    if (n_elements(allout) eq 0L) then begin
       photpath = catpath+'photometry/'
       all = file_search(photpath+'*.'+vv+'.fits.gz',count=nall)
       for ii = 0, nall-1 do begin
          phot = mrdfits(all[ii],1)
          nobj = n_elements(phot)
; build the output structure
          if (ii eq 0) then begin
             out1 = struct_trimtags(phot[0],select=['ediscsid',$
               'old_ediscsid','ra_hours','dec_deg'])
             out1 = create_struct(out1,mrd_struct('maut'+bands,replicate('-999.0',nband),1))
             out1 = create_struct(out1,mrd_struct('maute'+bands+'_apsim',replicate('-999.0',nband),1))
             out1 = create_struct(out1,mrd_struct('miso'+bands,replicate('-999.0',nband),1))
             out1 = create_struct(out1,mrd_struct('misoe'+bands+'_apsim',replicate('-999.0',nband),1))
             out1 = create_struct(out1,mrd_struct('map'+bands+'_1',replicate('-999.0',nband),1))
             out1 = create_struct(out1,mrd_struct('map'+bands+'_1e_apsim',replicate('-999.0',nband),1))
             out1 = create_struct(out1,mrd_struct('w'+bands,replicate('-999.0',nband),1))
             out1 = create_struct(out1,mrd_struct('flag'+bands,replicate('-999L',nband),1))
             out1 = create_struct(out1,mrd_struct('sexrhl'+bands,replicate('-999.0',nband),1))
;            for jj = 0, nband-1 do out1 = create_struct(out1,$
;              'maut'+bands[jj],0.0,'maute'+bands[jj]+'_apsim',0.0,$
;              'w'+bands[jj],0.0,'flag'+bands[jj],0.0,'sexrhl'+bands[jj],0.0)
          endif
          out = replicate(out1,nobj)
          struct_assign, phot, out, /nozero
          if (ii eq 0) then allout = out else allout = [allout,out]
       endfor
    endif 

    phot = struct_addtags(allapercor,allout)
    ngal = n_elements(phot)
       
; construct the final photometry; the idea here is to scale everything
; to the aperture-corrected I-band auto magnitude using the 1"
; diameter aperture color for crowded objects, (FLAG AND 1) EQ 1, and
; the isophotal aperture color otherwise
    outphot = mrd_struct('phot_'+bands,replicate('-999.0',nband),1)
    outphot = create_struct(outphot,mrd_struct('phot_'+bands+'_err',$
      replicate('-999.0',nband),1))
    outphot = replicate(outphot,ngal)
    outphot = struct_addtags(phot,temporary(outphot))

; isophotal colors and errors
    for kk = 0, nband-1 do begin
       apertag = tag_indx(outphot[0],'map'+bands[kk]+'_1') ; 1" aperture magnitude
       isotag = tag_indx(outphot[0],'miso'+bands[kk])      ; isophotal magnitude
       apertagerr = tag_indx(outphot[0],'map'+bands[kk]+'_1e_apsim')
       isotagerr = tag_indx(outphot[0],'misoe'+bands[kk]+'_apsim')
       
       flagtag = tag_indx(outphot[0],'flag'+bands[kk])
       outtag = tag_indx(outphot[0],'phot_'+bands[kk])
       outtagerr = tag_indx(outphot[0],'phot_'+bands[kk]+'_err')
       
       good = where((outphot.apercor gt -900.0) and $
         (outphot.mapi_1 gt 0.0) and (outphot.mapi_1 lt 90.0) and $
         (outphot.misoi gt 0.0) and (outphot.misoi lt 90.0) and $
         (outphot.mauti gt -900.0) and (outphot.mauti lt 90.0) and $
         (outphot.(apertag) gt 0.0) and (outphot.(apertag) lt 90.0) and $
         (outphot.(isotag) gt 0.0) and (outphot.(isotag) lt 90.0),ngood)
       itot = outphot[good].mauti + outphot[good].apercor ; total I-band magnitude
; blended
       blend = where((outphot[good].flagi and 1) eq 1,comp=clean)
       outphot[good[blend]].(outtag) = itot[blend] + $
         (outphot[good[blend]].(apertag)-outphot[good[blend]].mapi_1)
       outphot[good[blend]].(outtagerr) = outphot[good[blend]].(apertagerr)
; isolated          
       outphot[good[clean]].(outtag) = itot[clean] + $
         (outphot[good[clean]].(isotag)-outphot[good[clean]].misoi)
       outphot[good[clean]].(outtagerr) = outphot[good[clean]].(isotagerr)
       
       toobig = where(outphot[good[clean]].(outtag) gt 90)
       if (toobig[0] ne -1) then stop
    endfor

;   djs_plot, outphot.mautb, outphot.phot_b-outphot.mautb, xr=[18,28], yr=[-0.75,0.75], ps=4, ysty=3, xsty=3
;   djs_plot, outphot.mautv, outphot.phot_v-outphot.mautv, xr=[18,28], yr=[-0.75,0.75], ps=4, ysty=3, xsty=3
;   djs_plot, outphot.mautr, outphot.phot_r-outphot.mautr, xr=[18,28], yr=[-0.75,0.75], ps=4, ysty=3, xsty=3
;   djs_plot, outphot.mauti, outphot.phot_i-outphot.mauti, xr=[18,28], yr=[-0.75,0.75], ps=4, ysty=3, xsty=3
;   djs_plot, outphot.mautj, outphot.phot_j-outphot.mautj, xr=[18,28], yr=[-0.75,0.75], ps=4, ysty=3, xsty=3
;   djs_plot, outphot.mautk, outphot.phot_k-outphot.mautk, xr=[18,28], yr=[-0.75,0.75], ps=4, ysty=3, xsty=3

; add the photometric redshifts and IR photometry
    ircatpath = catpath+'ircatalogs/'
    if n_elements(ircat) eq 0L then begin
       all = file_search(ircatpath+'*mastertable.fits',count=nall)
       ircat = mrdfits(all[0],1)
       for ii = 1, nall-1 do ircat = [ircat,mrdfits(all[ii],1)]
    endif
    
    match2, strtrim(outphot.ediscsid,2), strtrim(ircat.ediscs_id,2), m1, m2
    good = where(m1 ne -1,comp=crap)

    outircat = im_empty_structure(ircat[0],empty_value=-999.0,ncopies=ngal)
    outircat[good] = im_struct_assign(ircat[m1[good]],outircat[good])

    out = struct_addtags(outphot,struct_trimtags(outircat,$
      except=['*iso*','ediscs_id*','ra','dec']))

    tags = tag_names(out)
    newtags = repstr(repstr(tags,'RA_HOURS','RA'),'DEC_DEG','DEC')
    out = im_struct_trimtags(out,select=tags,newtags=newtags)

; add the spectroscopic completeness weights from Bianca; there are
; duplicates, so to be safe, loop
    splog, 'Adding spectroscopic completeness'
    comp = read_completeness()
    keep = where(comp.z gt 0.0)
    comp = comp[keep]

    out = struct_addtags(out,replicate({spec_weight: 1.0, geom_weight: 1.0},ngal))
    for ii = 0, ngal-1 do begin
       cg = strtrim(comp.galaxy,2)
       m1 = where(strtrim(out[ii].ediscsid,2) eq cg or $
         strtrim(out[ii].old_ediscsid,2) eq cg,nm1)
; jm13jul02siena; from visual inspection, it looks like the multiple
; matches have different redshifts but the same spectroscopic weights,
; to just take the zeroth weights 
;      if nm1 gt 1 then begin
;         splog, 'Multiple matches!'
;         struct_print, comp[m1]
;      endif
       if nm1 ne 0 then begin
          out[ii].spec_weight = comp[m1[0]].spec_weight
          out[ii].geom_weight = comp[m1[0]].geom_weight
       endif
    endfor
    
;   match, strtrim(sout.ediscsid,2), strtrim(comp.galaxy,2), m1, m2
;   miss = lindgen(n_elements(comp))
;   remove, m2, miss
;   struct_print, comp[miss]
    
; write out
    outfile = mycatpath+'ediscs_all_photometry.'+vv+'.fits'
    im_mwrfits, out, outfile, clobber=clobber

; now match to the spec1d info file and write out; note:
; spherematching does not work here
    info = mrdfits(mycatpath+'ediscs_spec1d_info.fits.gz',1)
    ngal = n_elements(info)
    match2, strtrim(info.galaxy,2), strtrim(out.ediscsid,2), m1, m2
    
    good = where(m1 ne -1,comp=crap)
    sout = im_empty_structure(out[0],empty_value=-999.0,ncopies=ngal)
    sout[good] = im_struct_assign(out[m1[good]],sout[good])
    
;   m1 = lonarr(ngal)
;   for ii = 0, ngal-1L do m1[ii] = where(strtrim(info[ii].galaxy,2) eq $
;     strtrim(allout.ediscsid,2))
;   good = where(m1 ne -1,comp=crap)
       
;   spherematch, 15.0D*allout.ra_hours, allout.dec_deg, $
;     info.ra, info.dec, 3.0/3600.0, m1, m2
;   match, strtrim(allout.old_ediscsid,2), strtrim(info.galaxy,2), m1, m2

    outfile = mycatpath+'ediscs_photometry.'+vv+'.fits'
    im_mwrfits, sout, outfile, clobber=clobber

; prepend the SPEC1D info file
;      outphot = struct_addtags(info,struct_trimtags(outphot,$
;        except=['galaxy','cluster']))

;; correct for MW extinction and convert Vega --> AB
;       for kk = 0, nband-1 do begin
;          ftag = tag_indx(outphot[0],'phot_'+bands[kk])
;          good = where(outphot.(ftag) gt 0.0)
;          outphot[good].(ftag) = outphot[good].(ftag) + vega2ab[kk] - $
;            kl[kk]*info[good].ebv_mw
;       endfor
    
return
end

;;    if keyword_set(bamford) then begin
;;; read and parse Bamford's "total" and aperture photometry
;;; files, and append; correct for MW extinction
;;       photfile = catpath+'photometry_SED_table_noredcorr_2007-08-30.fits.gz'
;;       photfile_aper = catpath+'photometry_spec_comparison_table_2007-08-30.fits.gz'
;;
;;       splog, 'Reading '+photfile
;;       phot1 = mrdfits(photfile,1)
;;       ngal = n_elements(phot1)
;;       phot = im_struct_trimtags(phot1,select=['B','BERR','V','VERR','R','RERR',$
;;         'I','IERR','J','JERR','K','KERR'],newtags='phot_'+['b','b_err','v','v_err',$
;;         'r','r_err','i','i_err','j','j_err','k','k_err'])
;;       for i = 0L, ngal-1L do for j = 0L, n_tags(phot)-1L do if (phot[i].(j) eq 0.0) or $
;;         (abs(phot[i].(j)) gt 900.0) then phot[i].(j) = -999.0 else begin
;;          if odd(j) then phot[i].(j) = sqrt(phot[i].(j)^2.0 + phot1[i].total_corr_err^2.0) else $
;;            phot[i].(j) = phot[i].(j) + phot1[i].total_corr
;;       endelse
;;       for iband = 0L, n_elements(obsfilters)-1L do begin ; correct for MW extinction
;;          good = where((phot.(2*iband) gt -900.0) and (phot.(2*iband+1) gt -900.0))
;;          phot[good].(2*iband) = phot[good].(2*iband) + vega2ab[iband] - $ ; Vega-->AB
;;            k_lambda(weff[iband],/odonnell,r_v=3.1)*info[good].ebv_mw
;;       endfor
;;
;;       splog, 'Reading '+photfile_aper
;;       phot1_aper = mrdfits(photfile_aper,1)
;;       phot_aper = im_struct_trimtags(phot1_aper,select=['B','BERR','V','VERR','R','RERR',$
;;         'I','IERR','J','JERR','K','KERR'],newtags='phot_aper_'+['b','b_err','v','v_err',$
;;         'r','r_err','i','i_err','j','j_err','k','k_err'])
;;       for i = 0L, ngal-1L do for j = 0L, n_tags(phot_aper)-1L do if (phot_aper[i].(j) eq 0.0) or $
;;         (abs(phot_aper[i].(j)) gt 900.0) then phot_aper[i].(j) = -999.0
;;       for iband = 0L, n_elements(obsfilters)-1L do begin ; correct for MW extinction
;;          good = where((phot_aper.(2*iband) gt -900.0) and (phot_aper.(2*iband+1) gt -900.0))
;;          phot_aper[good].(2*iband) = phot_aper[good].(2*iband) + vega2ab[iband] - $ ; Vega-->AB
;;            k_lambda(weff[iband],/odonnell,r_v=3.1)*info[good].ebv_mw
;;       endfor
;;
;;       outphot = struct_addtags(phot,phot_aper)
;;       outphot = struct_addtags(struct_trimtags(info,$
;;         select=['galaxy','cluster','ra','dec','ebv_mw','z']),outphot)
;;
;;; read Bamford's uberfile to include the effective weights and SE flags 
;;       bamfordfile = catpath+'combined_spec_table_selected_2007-08-22.fits.gz'
;;       splog, 'Reading '+bamfordfile
;;       bamford = mrdfits(bamfordfile,1)
;;       
;;       bam = struct_trimtags(bamford,select=['phot_w*','*phot_flag*'])
;;       outphot = struct_addtags(outphot,bam)
;;
;;       outfile = mycatpath+'ediscs_photometry.bamford.fits'
;;
;;    endif else begin
;;       phot = im_empty_structure(allout[0],empty_value=-999.0,ncopies=ngal)
;;       phot[good] = im_struct_assign(allout[m1[good]],phot[good])
;;
;;       apercor = im_empty_structure(allapercor[0],empty_value=-999.0,ncopies=ngal)
;;       apercor[good] = im_struct_assign(allapercor[m1[good]],apercor[good])
;;       phot = struct_addtags(temporary(phot),apercor)
;;
