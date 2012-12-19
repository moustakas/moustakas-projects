pro sfrm_plots
; jm11may19ucsd - look at SFR/M vs M plots and compare with Salim

;;    jj = mrdfits('/Users/ioannis/research/data/isedfit_sfhgrid/sfhgrid04/bc03/charlot/chab_montegrid.fits.gz',1)
;;    these = where(jj.nburst gt 0,nthese)
;;    for ii = 0, nthese-1 do begin
;;       sfh = isedfit_reconstruct_sfh(jj[these[ii]],aburst=aburst)
;;;      ab[0:sfh[ii].nburst-1,ii] = aburst
;;       if (ii eq 0) then begin
;;          tb = jj[these[ii]].tburst[0:jj[these[ii]].nburst-1]
;;          fb = jj[these[ii]].fburst[0:jj[these[ii]].nburst-1]
;;          ab = aburst
;;       endif else begin
;;          tb = [tb,jj[these[ii]].tburst[0:jj[these[ii]].nburst-1]]
;;          fb = [fb,jj[these[ii]].fburst[0:jj[these[ii]].nburst-1]]
;;          ab = [ab,aburst]
;;       endelse
;;    endfor
;;    plot, tb, ab, psym=6, /ylog

    rr = mrdfits('./test2_bc03_chab_charlot_sfhgrid09.fits.gz',1)
;   rr = mrdfits('./test1_bc03_chab_charlot_sfhgrid08.fits.gz',1)
;   oo = mrdfits('./test2_bc03_chab_charlot_sfhgrid04.fits.gz',1)

    
    maxis = range(8.0,12.0,1000)
    
    djs_plot, rr.mass_avg, ((rr.sfr100_avg-rr.mass_avg)+9)>(-4), xsty=3, ysty=3, $
      psym=6, xr=[8.8,12], yr=[-4.2,0.7], sym=0.3
    med = im_medxbin(rr.mass_avg,rr.sfr100_avg-rr.mass_avg+9,0.05,/verb)
    oploterror, med.xbin, med.medy, med.sigy, color=djs_icolor('red'), errcolor=djs_icolor('red')
    djs_oplot, maxis, poly(maxis-10.0,[-9.83,-0.35])+9, line=0, thick=5, color='blue'     ; Salim+07

stop

;   djs_oplot, rr.mass_50, ((rr.sfr100_avg-rr.mass_50)+9)>(-4), psym=6, sym=0.3, color='orange'
;   djs_oplot, rr.mass_50, ((rr.sfr100_avg-rr.mass_50)+9)>(-4), psym=6, sym=0.3, color='cyan'
    djs_oplot, rr.mass_50, ((rr.sfr_avg-rr.mass_50)+9)>(-4), psym=6, sym=0.3, color='green'
;   djs_oplot, maxis, poly(maxis-10.0,[-9.83,-0.35])+9, line=0, thick=5, color='red'     ; Salim+07

    lo = where(maxis le 9.4) & hi = where(maxis ge 9.4)
    djs_oplot, maxis[lo], poly(maxis[lo]-10.0,[-9.65,-0.17])+9, line=0, thick=3, color='red'     ; Salim+07
    djs_oplot, maxis[hi], poly(maxis[hi]-10.0,[-9.87,-0.53])+9, line=0, thick=3, color='red'     ; Salim+07



    djs_plot, rr.mass_50, rr.sfr100-rr.mass_50, psym=6
    
    
    maxis = range(8.0,12.0,50)
    dfpsplot, 'junk.ps', /square, /color
    hogg_scatterplot, rr.mass_50, rr.sfr100_mode-rr.mass_50, psym=6, yr=[-14,-8], $
      xr=[8,11.9], ytickinterval=1, xsty=1, ysty=1, /internal, /outlier
    djs_oplot, maxis, -maxis, line=5, color='red'
    dfpsclose
    

stop    
    
    djs_plot, rr.mass_50, rr.sfr100_mode-rr.sfr100_avg, xsty=3, ysty=3

stop      
      

    med = im_medxbin(rr.mass_50,alog10(rr.chi2),0.1,/verb)
    djs_plot, rr.mass_50, alog10(rr.chi2), psym=6, sym=0.3, xsty=3, ysty=3, xr=[8.8,12]
    oploterror, med.xbin, med.medy, med.sigy, color=djs_icolor('red'), errcolor=djs_icolor('red')

stop    
    
    
    parent = mrdfits('parent.fits.gz',1)
    mpa = mrdfits(sdss_path(/mpa_dr7)+'mpamassoh_dr7_v5_2.fits.gz',1)
;   mpa = mrdfits(sdss_path(/mpa_dr4)+'mpamassoh_dr4_v5_1b.fits.gz',1)
    gd = where(mpa.ra gt 0.0 and mpa.sfr_flag eq 0)
    spherematch, parent.ra, parent.dec, mpa[gd].ra, mpa[gd].dec, $
      1.0/3600.0, m1, m2

    mass = rr[m1].mass_avg
    xx = rr[m1].sfr100_avg
;   xx = rr[m1].sfr100_50
    yy = mpa[gd[m2]].sfr_median
    med = im_medxbin(mass,yy-xx,0.1,/verb)

    djs_plot, mass, yy-xx, psym=6, sym=0.3, xsty=3, ysty=3, xr=[8.5,12], yrange=[-3,3]
    oploterror, med.xbin, med.medy, med.sigy, color=djs_icolor('red'), errcolor=djs_icolor('red')


stop    

return
end
    
