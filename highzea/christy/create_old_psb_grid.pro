
function make_sfh, tau, time = time, tburst=tburst, sf_peak=sf_peak, $
         sf_start = sf_start

if not keyword_set(sf_start) then sf_start = 30.0  ; starting SFR in M_sun/yr
if not keyword_set(sf_peak) then sf_peak = 100.0  ; SFR at peak of burst

dt = 50e6 ; time step in yrs
time = findgen(140) * dt ; time in yrs
time[0] = 1e6

tburst = 5.5e9    ; time burst occurs
sig_burst = 100e6 ; width of gaussian
contin_tau = 5e9  ; decay time of continuous SFR

sfh = sf_start * exp(-time/contin_tau)
burst = sf_peak * exp(-0.5 * ((time - tburst)/sig_burst)^2)

sfh = sfh + burst
postb = where(time ge tburst)
sfh[postb] = sfh[postb[0]] * exp(-(time[postb] - tburst)/tau)

tsf = total(sfh*dt)
sfh = sfh / tsf * 1e11 ; form  a  total mass of 10^11

;plot, time/1e9, sfh, xtitle = 'time (Gyr)', ytitle = 'SFR (M_sun/yr)'

return, sfh

end

pro create_old_psb_grid

sf_start = 30.0
sf_peak = 100.0
modelname = 'psb_models_bc03_z05_cmd.fit'

; Read SDSS filters
fdir = '/home/tremonti/local/idlutils/data/filters/'
readcol, fdir + 'sdss_jun2001_g_atm.dat', gwl, gfl
readcol, fdir + 'sdss_jun2001_r_atm.dat', rwl, rfl
readcol, fdir + 'sdss_jun2001_i_atm.dat', iwl, ifl


; Read Model Grid
m_ssp = mrdfits(hizea_path(/models) + $
        'bc03_pad94_chab_z05_ssp_uvpatch_mmtres.fit', 1)

m_tab = mrdfits(hizea_path(/models) + 'bc03_extras_z050.fit', 1)

minwl = 2000
maxwl = 9000 
normwl = [5200, 5800]

tau_str = ['10', '25', '50', '100', '200', '300', '400', '500']
tau_val = tau_str * 1e6
ntau = n_elements(tau_val)

out_ages = 5.5e9 + [0.0, 25.0, 50.0, 75.0, 100.0, 125.0, 150.0, 175.0, $
                    200.0, 250.0, 300.0, 350.0, 400.0, $
                    500.0, 600.0, 700.0, 800.0, 900.0, 1000.0, $
                    1100.0, 1200.0, 1300.0, 1400.0, 1500.0, 1750, $
                    2000.0, 3000.0, 4000.0, 5000.0] * 1e6 ; in yrs

nage = n_elements(out_ages)

ok = where(m_ssp.wave ge minwl and m_ssp.wave le maxwl, nwave)

m_burst = {burst_age: out_ages-5.5e9, burst_tau: tau_val, $
           wave:double(m_ssp.wave[ok]), flux: dblarr(nwave, nage, ntau), $
           mass: dblarr(nage, ntau), sfr: dblarr(nage, ntau), $
           log_nlyc: dblarr(nage, ntau), $ 
           ha_lum: dblarr(nage, ntau), ha_ew: dblarr(nage, ntau), $
           OII_lum: dblarr(nage, ntau), OII_ew: dblarr(nage, ntau), $   
           burst_norm: dblarr(nage, ntau), $
           plaw: dblarr(nwave), qso_comp: dblarr(nwave), $
           nebular: dblarr(nwave, nage, ntau), $
           sf_peak: sf_peak, sf_start: sf_start, $
           gri_01: fltarr(3, nage, ntau)} 

; Add powerlaw
m_burst.plaw = (m_burst.wave/5500.0)^(-1.6)

; Store qso composite spectrum
readcol, hizea_path(/sdss) + 'qso_composite.txt', qwl, qfl, skip = 23
vactoair, qwl
linterp, qwl, qfl, m_burst.wave, mcqso
m_burst.qso_comp = mcqso / median(mcqso[where(m_burst.wave gt 5450 and $
        m_burst.wave lt 5550)])

; Create nebular spectrum
readcol, hizea_path(/models) + 'hummer_storey.txt', num, bwl, bflux
; nomarlize to H-alpha rather than H-beta
 blfux = bflux / 3.0
 
sig_a = 92.0/2.99792e5 * bwl ; line width in A
ampl = bflux / (sqrt(2*!PI)*sig_a)

nebular = 0.0
for ii = 0, n_elements(bwl) - 1 do $
  nebular = nebular + $
    ampl[ii] * exp(-1*(m_burst.wave - bwl[ii])^2/(2*sig_a[ii]^2))

color = [!white, !magenta, !pink, !dred, !red, !orange, !yellow, !green, $
         !cyan, !blue, !dblue, !purple, !brown, !lgray, !dgray]

color = [color, color]

l_sun = 3.827d33 ; erg/s
k98_ha2sfr = 7.9d-42  ; Kennicutt 1998 Salpeter IMF 0.1 - 100
salp_to_chab = 1.8    ; to convert to Chabrier divide SFRs by this
john_o2sfr = 10.D^(0.762) * 1d-41  ; From Moustakas et al. 2006 (highest L_B)

ha_wl = where(m_burst.wave gt 6523 and m_burst.wave lt 6603)
o2_wl = where(m_burst.wave gt 3687 and m_burst.wave lt 3767)

; Loop through different Tau models

for itau = 0, ntau - 1 do begin
  
  print, tau_val[itau]
  sfh = make_sfh(tau_val[itau], time=time, sf_start=sf_start, sf_peak=sf_peak) 

  linterp, time, sfh, out_ages, osfh  ; time & age in Gyr
  m_burst.sfr[*,itau] = osfh

  plot, time/1e9, sfh, /xs, xtitle = 'time', ytitle = 'SFR'
  oplot, out_ages/1e9, osfh, psym=4, color=!red, thick=2
 
  wait, 2

  ;Normalize spectra & plot
  plot, [0,10], [0,10], /nodat, xr=[2500,5000], /xs, yr=[0,3], $
        title = 'Tau = ' + tau_str[itau] + ' Myr'         
 
  for iage = 0, nage - 1 do begin

    spec = add_past_sf(time, sfh, m_ssp, out_ages[iage], $
                       wlrange=[minwl, maxwl], wave = wave, mass=mass, $
                       ok=ok, ssp_nlyc=10.D^m_tab.nly, total_nlyc=tnlyc)

    norm = median(spec[where(wave ge normwl[0] and wave le normwl[1])])

    m_burst.flux[*,iage,itau] = spec /norm
    m_burst.burst_norm[iage,itau] = norm * l_sun
    m_burst.mass[iage,itau] = mass
    m_burst.log_nlyc[iage,itau] = alog10(tnlyc)
    oplot, wave, spec/norm, color=color[iage]

    ; Store nebular models
    ha_lum_per_lyc = 1.46d-12 ; calculated for T_e=10^4 n_e=10^4
    f_esc = 0.10
    ha_lum_per_lyc = ha_lum_per_lyc * (1 - f_esc)
    m_burst.nebular[*,iage,itau] = nebular * tnlyc * ha_lum_per_lyc / $
      l_sun / norm

    ; Measure line lums, EWs   
    ;m_burst.ha_lum[iage,itau] = m_burst.sfr[iage,itau] / k98_ha2sfr * $
    ;  salp_to_chab / l_sun ; in units of L_sun
    m_burst.ha_lum[iage, itau] = tnlyc * ha_lum_per_lyc / l_sun

    ;m_burst.OII_lum[iage,itau] = m_burst.sfr[iage,itau] / john_o2sfr * $
    ;  salp_to_chab / l_sun ; in units of L_sun
    m_burst.OII_lum[iage,itau] = m_burst.ha_lum[iage,itau] * k98_ha2sfr $
        / john_o2sfr / l_sun

    m_burst.ha_ew[iage,itau] = m_burst.ha_lum[iage,itau] / $
      median(m_burst.flux[ha_wl,iage,itau] * norm)
    m_burst.OII_ew[iage,itau] = m_burst.OII_lum[iage,itau] / $
      median(m_burst.flux[o2_wl,iage,itau] * norm)

    ; Measure gri colors at z=0.1

    gmag = filtermag(wave*1.1, spec, gwl, gfl) + 4.72     
    rmag = filtermag(wave*1.1, spec, rwl, rfl) + 4.72    
    imag = filtermag(wave*1.1, spec, iwl, ifl) + 4.72    
      
    m_burst.gri_01[*,iage,itau] = [gmag,rmag,imag]
  endfor

  legend, string((out_ages - 5.5e9)/1e6, form='(I4)'), color=color, $
          linestyle=0, /right


  wait, 2
endfor

mwrfits, m_burst, hizea_path(/models) + modelname, /create

end

