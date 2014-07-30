pro parse_xray_mass, debug=debug, clobber=clobber
; jm14jul28siena - parse the X-ray mass Table 4 from Megan's
; paper to get M(500) and write out a new FITS table
    
    propath = bcgmstar_path(/propath)
    outfile = propath+'xray_mass_parsed.fits'

    bigG = 6.67384D-11*1.99D30/1D3^2/3.086D22 ; [convert from m^3/kg/s^2 to km^2*Mpc/Msun/s2]

    xray = rsex(propath+'xray_mass.sex')
    ncl = n_elements(xray)
    out = struct_addtags(xray,replicate({m500: 0.0, m500_err: 0.0},ncl))

;; check to make sure my M(2500) matches the Table - yes!
;    xx_2500 = xray.r2500/xray.rs
;    m2500 = 4.625*xray.vmax^2*xray.r2500/bigG*(alog(1+xx_2500)/xx_2500-1.0/(1+xx_2500))
;    niceprint, xray.cluster, m2500, xray.m2500
    
    xx_500 = xray.r500/xray.rs
    out.m500 = 4.625*xray.vmax^2*xray.r500/bigG*(alog(1+xx_500)/xx_500-1.0/(1+xx_500))

; use Monte Carlo to get the uncertainty
    nmc = 500
    for ic = 0, ncl-1 do begin
       vmax_mc = xray[ic].vmax+randomn(seed,nmc)*xray[ic].vmax_err/2
       r500_mc = xray[ic].r500+randomn(seed,nmc)*xray[ic].r500_err/2
       rs_mc = xray[ic].rs+randomn(seed,nmc)*xray[ic].rs_err/2
       xx_500_mc = r500_mc/rs_mc
       good = where(xx_500_mc gt 0)
       
       m500_mc = 4.625*vmax_mc[good]^2*r500_mc[good]/bigG*(alog(1+xx_500_mc[good])/$
         xx_500_mc[good]-1.0/(1+xx_500_mc[good]))
;      good = where(finite(m500_mc) and m500_mc gt 0)
       out[ic].m500_err = djsig(m500_mc)

       splog, xray[ic].cluster, out[ic].m500, out[ic].m500_err, out[ic].m500/out[ic].m500_err

       if keyword_set(debug) then begin
          im_plothist, m500_mc[good]/1E14, bin=0.1
          djs_oplot, out[ic].m500*[1,1]/1E14, !y.crange, color='green'
          djs_oplot, (out[ic].m500+out[ic].m500_err)*[1,1]/1E14, !y.crange, color='green', line=5
          djs_oplot, (out[ic].m500-out[ic].m500_err)*[1,1]/1E14, !y.crange, color='green', line=5
          cc = get_kbrd(1)
       endif
    endfor

    im_mwrfits, out, outfile, clobber=clobber
    
return
end
