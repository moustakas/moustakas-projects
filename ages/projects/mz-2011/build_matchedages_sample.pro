pro build_matchedages_sample, q0=q0, clobber=clobber
; jm11jun07ucsd - build an SDSS sample that matches AGES in each
; redshift bin 

    mzpath = mz_path()
;   matchpath = mz_path(/matchedages)
    zbins = mz_zbins(nz)

; read the data    
    agesancillary = read_mz_sample(/mzhii_ancillary)
    agesispec = read_mz_sample(/mzhii_ispec)
    agesmass = read_mz_sample(/mzhii_mass)
    agesohdust = read_mz_sample(/mzhii_log12oh)

    sdssancillary = read_mz_sample(/mzhii_ancillary,/sdss)
    sdssispec = read_mz_sample(/mzhii_ispec,/sdss)
    sdssmass = read_mz_sample(/mzhii_mass,/sdss)
    sdssohdust = read_mz_sample(/mzhii_log12oh,/sdss)

; match in mass, luminosity, and log[EW(Hb)]
    hmin = [9.0,-23.5,0.1]
    hmax = [12.0,-18.5,2.0]
    binsz = [0.15,0.15,0.05]
    nbins = round((hmax-hmin)/binsz)
    nmatched = 8000 ; number of "mock" galaxies

    if (n_elements(q0) eq 0) then q0 = 1.5 ; from Cool+11
    qz0 = 0.1
    evolsuffix = 'q'+string(q0,format='(F4.2)')+'-qz0.1'

; deal with each calibration separately 
    ncalib = 3
    for ii = 0, ncalib-1 do begin
       case ii of
          0: begin
             t04 = 1 & m91 = 0 & kk04 = 0
             calib = 't04'
          end
          1: begin
             t04 = 0 & m91 = 1 & kk04 = 0
             calib = 'm91'
          end
          2: begin
             t04 = 0 & m91 = 0 & kk04 = 1
             calib = 'kk04'
          end
       endcase

       out_template = {q0: q0, qz0: qz0, z: 0.0, oh: 0.0, $
         oh_err: 0.0, weight: 1.0, mb: 0.0, mass: 0.0, ewhb: 0.0}
       outfile = mzpath+'matchedages_'+calib+'_'+evolsuffix+'.fits'
;      outfile = matchpath+'matchedages_mass_zbin'+string(findgen(nzbin)+1,$
;        format='(I0)')+'_'+evolsuffix+'_'+calib+'.fits'
       if im_file_test(outfile+'.gz',clobber=clobber) then return

       sinfo = mzlz_grab_info(sdssohdust,sdssancillary,sdssmass,$
         t04=t04,m91=m91,kk04=kk04,/nolimit)
       for iz = 0, nz-1 do begin
          ainfo = mzlz_grab_info(agesohdust,agesancillary,agesmass,$
            t04=t04,m91=m91,kk04=kk04,zmin=zbins[iz].zlo,$
            zmax=zbins[iz].zup,/nolimit)

; build an SDSS sample matched in stellar mass, luminosity and EW(Hb) 
          amass = ainfo.mass
          amb = ainfo.mb_ab
          aewhb = alog10(agesispec[ainfo.indx].h_beta_ew[0])

          smass = sinfo.mass
          smb = sinfo.mb_ab - q0*(zbins[iz].zbin-qz0)
          sewhb = alog10(sdssispec[sinfo.indx].h_beta_ew[0])

          hages = hist_nd(transpose([[amass],[amb],[aewhb]]),nbins=nbins,min=hmin,max=hmax,rev=arev)
          hsdss = hist_nd(transpose([[smass],[smb],[sewhb]]),nbins=nbins,min=hmin,max=hmax,rev=srev)
          ntotbins = cmproduct(nbins)
          nperbin = round(nmatched*hages/total(hages))

; for testing, get the indices of the AGES galaxies used 
          delvarx, indx
          for kk = 0, ntotbins-1 do begin
             if (arev[kk+1] gt arev[kk]) and (nperbin[kk] gt 0) then begin
                indx1 = arev[arev[kk]:arev[kk+1]-1]
                if (n_elements(indx) eq 0) then indx = indx1 else $
                  indx = [indx,indx1]
             endif
          endfor

; now select a Monte Carlo subset of SDSS galaxies matching the AGES
; distributions 
          delvarx, sel
          for kk = 0, ntotbins-1 do begin
             if (srev[kk+1] gt srev[kk]) and (nperbin[kk] gt 0) then begin
                inbin = srev[srev[kk]:srev[kk+1]-1]
                sel1 = inbin[round(randomu(seed,nperbin[kk])*n_elements(inbin))]
                if (n_elements(sel) eq 0) then sel = sel1 else $
                  sel = [sel,sel1]
             endif
          endfor

          out1 = replicate(out_template,n_elements(sel))
          out1.z = zbins[iz].zbin          
          out1.oh = sinfo.oh[sel]
          out1.oh_err = sinfo.oh_err[sel]
          out1.mb = smb[sel]
          out1.mass = smass[sel]
          out1.ewhb = sewhb[sel]
          if (n_elements(out) eq 0L) then out = out1 else out = [out,out1]

;         plothist, amb[indx], /peak, bin=0.1               
;         plothist, smb[sel], /peak, bin=0.1, /over, thick=8

;         plothist, alog10(sdssispec[sel].h_beta_ew[0]), /peak, bin=binsz[2]
;         plothist, alog10(agesispec[ainfo.indx].h_beta_ew[0]), /over, bin=binsz[2], color=djs_icolor('red'), /peak
;         cc = get_kbrd(1)
;
;         plothist, sdssancillary[sel].k_ubvrijhk_absmag_00[1], /peak, bin=binsz[1]
;         plothist, ainfo.mb_ab, /over, bin=binsz[1], color=djs_icolor('red'), /peak
;         cc = get_kbrd(1)
;
;         plothist, sdssmass[sel].mass_50, /peak, bin=binsz[0]
;         plothist, ainfo.mass, /over, bin=binsz[0], color=djs_icolor('red'), /peak
;         cc = get_kbrd(1)
       endfor ; redshift loop
       im_mwrfits, out, outfile, clobber=clobber
    endfor    ; calibration loop

return
end
