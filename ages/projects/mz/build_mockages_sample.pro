;+
; NAME:
;   BUILD_MOCKAGES_SAMPLE
;
; PURPOSE:
;   Simulate the AGES sample at different redshifts using the SDSS. 
;
; INPUTS: 
;
; KEYWORD PARAMETERS: 
;   select - select a sample from the SDSS in each redshift bin that
;     passes the AGES flux limits
;   evolve - allow the SDSS 
;
; OUTPUTS: 
;
; COMMENTS:
;
; MODIFICATION HISTORY:
;
; Copyright (C) 2010, John Moustakas
; 
; This program is free software; you can redistribute it and/or modify 
; it under the terms of the GNU General Public License as published by 
; the Free Software Foundation; either version 2 of the License, or
; (at your option) any later version. 
; 
; This program is distributed in the hope that it will be useful, but 
; WITHOUT ANY WARRANTY; without even the implied warranty of
; MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
; General Public License for more details. 
;-

pro build_mockages_sample, select=select, noevolve=noevolve, kcorr=kcorr
; jm09mar19nyu - Simulate AGES at different redshifts using the SDSS 
; jm10oct16ucsd - major rewrite
; echo "build_mockages_sample, /select, /kcorr" | idl > & log.noevolve &
; echo "build_mockages_sample, /select, /kcorr, /evolve" | idl > & log.evolve &

; read the SDSS photometry and SDSS/MZ parent sample
    mockpath = ages_path(/projects)+'mz/mockages/'
;   sdss = read_mz_sample(/sdss,/parent)
    sdss = read_mz_sample(/sdss,/mzhii_ancillary)
    sdssispec = read_mz_sample(/sdss,/mzhii_ispec)
    ohdust = read_mz_sample(/sdss,/mzhii_log12oh)
    agesispec = read_mz_sample(/mzhii_ispec)
    nsdss = n_elements(sdss)

; initialize the parameters of the simulation
    vname = mz_vname()
    if keyword_set(noevolve) then begin
       q0 = 0.0
       qz0 = 0.0
       evolsuffix = 'noevol'
       beta = 0.0
    endif else begin
       q0 = 1.5 ; from Cool+11
       qz0 = 0.1
       beta = 3.5
       evolsuffix = 'q1.50-beta3.50'
    endelse

    zbins = mz_zbins(nzbin)
;   zbins = zbins[0] & nzbin = 1
    struct_print, zbins

; output file names
    select_outfile = mockpath+'mockages_select_zbin'+$
      string(findgen(nzbin)+1,format='(I0)')+'_'+evolsuffix+'.fits'
    mockages_outfile = repstr(select_outfile,'_select','')
    
    mzfilters = mz_filterlist()
    nfilt = n_elements(mzfilters)
    ifaint = mz_ifaint(ibright=ibright,select_filter=select_filter)
    iselect = (where(strtrim(mzfilters,2) eq select_filter))[0]
    hbcut = mz_hbcut() ; H-beta flux cut

; ##################################################
; select the sample in each redshift bin
    if keyword_set(select) then begin
       thesemin = 10000L ; minimum number of mock SDSS galaxies per iteration
       nitermax = 1E5    ; maximum number of iterations
       nmockmax = 5000L  ; desired number of galaxies per redshift bin

; output selection structure
       select_info_template = {$
         object_position:             0L,$ 
         z:                          0.0,$
         mock_maggies:     fltarr(nfilt),$
         mock_ivarmaggies: fltarr(nfilt),$
         h_beta:                     0.0,$
         h_beta_ew:                  0.0}

; loop on each redshift bin
       t1 = systime(1)
       for iz = 0, nzbin-1 do begin
; iterate until we build a sufficiently large sample
          t0 = systime(1)
          iiter = 0L & nmock = 0L & thatsenough = 0
          while (thatsenough eq 0L) do begin
             print, format='("Iter ",I0,"/",I0,", NMOCK = ",I0," ",A5,$)', $
               iiter+1L, nitermax, nmock, string(13b)
             
; choose a random subset of SDSS galaxies (with replacement) with a
; probability equal to 1/Vmax
             vmax = sdss.vmax_evol ; h=0.7, but it doesn't matter
             maxvmax = max(1.0/vmax)
             chances = (1.0/vmax)/maxvmax

; to minimize overhead, build up THESE until it has at least THESEMIN
; elements
             these = where((randomu(seed,nsdss) lt chances),nthese)
             while nthese lt thesemin do begin
                these1 = where((randomu(seed,nsdss) lt chances),nthese1)
                if (nthese1 ne 0L) then these = [these,these1]
                nthese = n_elements(these) ; & print, nthese
             endwhile
             
             if (nthese gt 0L) then begin
; distribute galaxies uniformly in volume
                vinner = (lf_comvol(zbins[iz].zlo))[0] ; h=1
                vouter = (lf_comvol(zbins[iz].zup))[0] ; h=1
                vol = vinner+(vouter-vinner)*randomu(seed,nthese)
                zages = lf_vtoz(vol)
                
; compute the rest-frame AGES magnitudes
;               sdss_nmaggies = sdss[these].k_maggies*1D9           ; nanomaggies
;               sdss_ivarnmaggies = sdss[these].k_ivarmaggies/1D9^2 ; nanomaggies
                sdss_nmaggies = sdss[these].k_maggies[0:4]*1D9           ; nanomaggies
                sdss_ivarnmaggies = sdss[these].k_ivarmaggies[0:4]/1D9^2 ; nanomaggies
                mock_mag = sdss2ages(sdss[these].z,zages,nmgy=sdss_nmaggies,$
                  ivar=sdss_ivarnmaggies,vname=vname,q0=q0,q1=0.0,qz0=qz0,$
                  agesphot_ivar=mock_ivarmag,out_filterlist=mzfilters,$
                  in_filterlist=sdss_filterlist(),chi2=chi2)

; only keep objects with Ibright<I<Ifaint and H(beta)>H(beta)_cut,
; allowing for L(Hb) evolution
;               dm_sdss = lf_distmod(sdss[these].z,omega0=omega0,omegal0=omegal0)
;               dm_ages = lf_distmod(zages,omega0=omega0,omegal0=omegal0)
                hb = sdssispec[these].h_beta[0]*(dluminosity(sdss[these].z)/dluminosity(zages))^2.0*(1.0+zages)^beta
                keep = where((mock_mag[iselect,*] ge ibright) and $
                  (mock_mag[iselect,*] le ifaint) and (hb gt hbcut),nkeep)
                if (nkeep ne 0L) then begin
                   select_info1 = replicate(select_info_template,nkeep)
                   select_info1.object_position = these[keep]
                   mock_maggies = 10.0^(-0.4*mock_mag[*,keep])
                   select_info1.z = zages[keep]
                   select_info1.mock_maggies = mock_maggies
                   select_info1.mock_ivarmaggies = mock_ivarmag[*,keep]/(0.4*alog(10.0)*mock_maggies)^2.0
                   select_info1.h_beta = hb[keep]
                   select_info1.h_beta_ew = sdssispec[these[keep]].h_beta_ew[0]
                   if (n_elements(select_info) eq 0L) then select_info = select_info1 else $
                     select_info = [temporary(select_info),select_info1]
                endif 
             endif 
             iiter++ ; increase the counter
             nmock = n_elements(select_info)
             if (iiter ge nitermax) or (nmock ge nmockmax) then thatsenough = 1
          endwhile 
          splog, 'Number of iterations = ', iiter
          splog, 'Number of mock galaxies = ', nmock
          splog, 'Time for this redshift bin = '+strtrim(string((systime(1)-t0)/60.0,$
            format='(F12.1)'),2)+' minutes'
          im_mwrfits, select_info, select_outfile[iz], /clobber
          delvarx, select_info
       endfor ; zbins
       splog, 'Time to build sample = '+strtrim(string((systime(1)-t1)/60.0,$
         format='(F12.1)'),2)+' minutes'
    endif 
    
; ##################################################
;  compute K-corrections and write out the final big table of global
;  properties, rest-frame quantities, and oxygen abundances
    if keyword_set(kcorr) then begin
       hbbinsz = 0.05
       hbmin = -0.5
       hbmax = 2.0
       nfinal = 7500 ; final number of galaxies per redshift bin
       
       for iz = 0, nzbin-1 do begin
          splog, 'Computing K-corrections for zbin ', iz
          splog, 'Reading '+select_outfile[iz]+'.gz'
          info = mrdfits(select_outfile[iz]+'.gz',1)

; match the EW(Hb) distributions
          these = where((agesispec.z gt zbins[iz].zlo) and $
            (agesispec.z lt zbins[iz].zup),nthese)
          ewhbages = histogram(alog10(agesispec[these].h_beta_ew[0]),$
            min=hbmin,max=hbmax,bin=hbbinsz)
          ewhbsdss = histogram(alog10(info.h_beta_ew),$
            min=hbmin,max=hbmax,bin=hbbinsz,rev=rev)
          nbins = n_elements(ewhbages)

;         nperbin = ewhbages
          nperbin = round(nfinal*ewhbages/total(ewhbages))

          nfinal = long(total(nperbin))
          delvarx, sel
          for ii = 0, nbins-1 do begin
             if (rev[ii+1] gt rev[ii]) and (nperbin[ii] gt 0) then begin
                inbin = rev[rev[ii]:rev[ii+1]-1]
                sel1 = inbin[round(randomu(seed,nperbin[ii])*n_elements(inbin))]
;               if (nperbin[ii] gt 100) then stop
                if (n_elements(sel) eq 0) then sel = sel1 else $
                  sel = [sel,sel1]
             endif
          endfor

;         im_plothist, alog10(info[sel].h_beta_ew), bin=hbbinsz, $
;           histmin=hbmin, histmax=hbmax, t1, t2
;         im_plothist, alog10(agesispec[these].h_beta_ew[0]), bin=hbbinsz, $
;           histmin=hbmin, histmax=hbmax, /over, color='red', j1, j2

stop          
          
          kcorr = mz_kcorrect(info[sel].z,info[sel].mock_maggies,$
            info[sel].mock_ivarmaggies,filterlist=mzfilters)
          out = struct_addtags(info[sel],kcorr)
          out = struct_addtags(temporary(out),ohdust[info[sel].object_position])
          im_mwrfits, out, mockages_outfile[iz], /clobber

stop
       endfor       
    endif
    
return
end
    
