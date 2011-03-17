;+
; NAME:
;   COSMICIMF_LIMITS
;
; PURPOSE:
;   Compute the limiting SFR and stellar mass as a function of
;   redshift.  
;
; INPUTS: 
;
; KEYWORD PARAMETERS: 
;
; OUTPUTS: 
;
; COMMENTS:
;    If I_lim=19.95 is the apparent magnitude limit of the survey, I
;    is the apparent magnitude of any source, and M_g is its absolute
;    magnitude (computed using K-correct), then the limiting absolute
;    g-band magnitude for that source is given simply by: 
;       M_g_lim = M_g - (I-I_lim) 
;
;    The limiting stellar mass, log_M_lim, is similarly given by:
;    log_M_lim = log_M + 0.4(I-I_lim) 
;
; MODIFICATION HISTORY:
;   J. Moustakas, 2010 Mar 18, UCSD
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

pro cosmicimf_limits, lf24=lf24

    cosmicimfpath = ages_path(/projects)+'cosmicimf/'
    parent = read_cosmicimf_sample()
    sfrs = read_cosmicimf_sample(/sfrs)
    mass = read_cosmicimf_sample(/mass)
    zbins = cosmicimf_zbins(nzbins,lf24=lf24)

    quant = 0.5
    I_lim = 19.95
    zbin = 0.02
    zaxis = im_array(0.05,0.75,zbin)
    nzaxis = n_elements(zaxis)

; define the limiting luminosity or mass at the *upper* end of each
; redshift bin!    
    if keyword_set(lf24) then begin
       limits = struct_addtags(zbins,replicate({l24_lim: 0.0},nzbins))

       good = where(sfrs.mips)
       agn = where(sfrs.mips and sfrs.agn)
       stats = im_medxbin(sfrs[good].z,sfrs[good].l24,0.05)
       offset = -0.05
       im_poly_iter, stats.binctr, stats.sigy05+offset, 4, coeff=coeff
       limits.l24_lim = poly(limits.zup,coeff)
       
; simple QAplot
       psfile = cosmicimfpath+'limits_lf24.ps'
       im_plotconfig, 0, pos, psfile=psfile
       djs_plot, sfrs[good].z, sfrs[good].l24, position=pos, $
         psym=6, sym=0.2, ysty=3, xsty=1, xtitle='Redshift', $
         ytitle='log_{10} [L(24 \mu'+'m)/L_{\odot}]', $
         xrange=[0.0,0.8]
       djs_oplot, sfrs[agn].z, sfrs[agn].l24, ps=6, sym=0.4, color='blue'
       djs_oplot, stats.binctr, stats.sigy05+offset, psym=symcat(15), sym=2, color='red'
       for kk = 0, nzbins-1 do begin
          djs_oplot, limits[kk].zlo*[1,1], !y.crange, line=5, color='orange'
          djs_oplot, limits[kk].zup*[1,1], !y.crange, line=5, color='orange'
          djs_oplot, [limits[kk].zlo,limits[kk].zup], limits[kk].l24_lim*[1,1], $
            color='orange', line=5
       endfor
       djs_oplot, zaxis, poly(zaxis,coeff), color='blue'
       im_plotconfig, psfile=psfile, /psclose, /gzip
       struct_print, limits

; write out
       outfile = cosmicimfpath+'limits_lf24.fits'
       im_mwrfits, limits, outfile, /clobber

    endif else begin

; define the stellar mass limits for each IMF
       peg = read_cosmicimf_sample(/pegase)
       imf = strtrim(peg.imf,2)
       nimf = n_elements(peg)
       limits = struct_addtags(zbins,replicate({mass_lim: 0.0},nzbins))

       outfile = cosmicimfpath+'limits_mf.fits'
       psfile = cosmicimfpath+'limits_mf.ps'
       im_plotconfig, 0, pos, psfile=psfile, ymargin=[0.6,1.1]
;      for ii = 0, 0 do begin
       for ii = 0, nimf-1 do begin

; compute the limiting stellar mass on a fine redshift bin          
          masslim = mass.mass[ii] + 0.4*(parent.i_tot-I_lim)
          minmass = zaxis*0.0
          for jj = 0, nzaxis-2 do begin
             these = where((parent.z ge zaxis[jj]) and (parent.z lt zaxis[jj+1]),nthese)
             weight1 = parent[these].final_weight
             masslim1 = masslim[these]
             minmass[jj] = im_quantile(masslim1,weight1,quant=quant)
          endfor

; fit a smooth polynomial
          good = where((zaxis gt 0.05) and (zaxis lt 0.75) and $
            (minmass ne 0.0),ngood)
          im_poly_iter, zaxis[good], minmass[good], 4, coeff=coeff
          
; interpolate the polynomial model at the *upper* redshift bin to get
; the limiting mass in each bin
          limits.mass_lim = poly(limits.zup,coeff)

; simple QAplot
          title = strtrim(repstr(peg[ii].imf,'cosmicimf_','\alpha_{2}='),2)
          djs_plot, parent.z, mass.mass[ii], position=pos, $
            xsty=1, ysty=3, psym=6, sym=0.2, title=title, $
            xrange=[0.0,0.8], xtitle='Redshift', ytitle=cosmicimf_masstitle()
          djs_oplot, zaxis[0:nzaxis-2], minmass[0:nzaxis-2], color='red', psym=-8
          for kk = 0, nzbins-1 do begin
             djs_oplot, limits[kk].zlo*[1,1], !y.crange, line=5, color='orange'
             djs_oplot, limits[kk].zup*[1,1], !y.crange, line=5, color='orange'
             djs_oplot, [limits[kk].zlo,limits[kk].zup], limits[kk].mass_lim*[1,1], $
               color='orange', line=5
          endfor
          djs_oplot, zaxis, poly(zaxis,coeff), color='blue'
          struct_print, limits

          splog, 'Writing '+outfile
          mwrfits, limits, outfile, create=(ii eq 0)
       endfor ; IMF loop
       im_plotconfig, psfile=psfile, /psclose, /gzip
       spawn, 'gzip -f '+outfile, /sh
    endelse

return
end
    

;
;       limits = struct_addtags(zbins,replicate({l1500_lim: 0.0},nzbins))
;       limits_coeff = [8.1,2.5]
;       limits.l1500_lim = poly(zbins.zup,limits_coeff)
;
;; QAplot       
;       zaxis = im_array(0.05,0.75,0.02)
;       djs_plot, sfrs.z, sfrs.l1500, ps=6, sym=0.2, ysty=3, xsty=3
;       oplot, zaxis, poly(zaxis,limits_coeff)
;
