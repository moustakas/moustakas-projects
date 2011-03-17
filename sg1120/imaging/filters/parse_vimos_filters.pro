pro parse_vimos_filters, debug=debug
; jm08may14nyu - write out the VIMOS filters in K-correct format

    path = '/Users/ioannis/home/research/projects/sg1120/filters/'

    quad = ['Q1','Q2','Q3','Q4']
    bandpass = ['U','B','V','R','I']
    oldfilter = bandpass+'_total_transmission.dat'
    newfilter_root = 'vimos_'+['U','B','V','R','I']

    for ii = 0L, n_elements(oldfilter)-1L do begin

       readcol, path+'ORIGINAL/'+oldfilter[ii], wave, q1, $
         q2, q3, q4, /silent, comment='#'
       qresp = [[q1],[q2],[q3],[q4]]
       nn = n_elements(wave)

       for iq = 0L, n_elements(quad)-1L do begin

;         openw, lun, path+newfilter_root[ii]+'_'+quad[iq]+'.par', /get_lun
          
          hdr0 = '# Units:'
          hdr0 = [hdr0, '#  "lambda" is in Angstroms']
          hdr0 = [hdr0, '#  "pass" is the contribution to the detector signal per photon']
          hdr0 = [hdr0, '#']
          hdr0 = [hdr0, '# Bandpass Name(s): '+bandpass[ii]]
          hdr0 = [hdr0, '#']
          hdr0 = [hdr0, '# Instrument: VIMOS/Quadrant '+quad[iq]]
          hdr0 = [hdr0, '#']
          hdr0 = [hdr0, '# Determined by: Mieske et al. 2007:']
          hdr0 = [hdr0, '#   http://www.eso.org/sci/facilities/paranal/instruments/vimos/inst/filters_2007-04-04/memo_filtrans_March26th.pdf']
          hdr0 = [hdr0, '#']
          hdr0 = [hdr0, '# Date of determination: 2007-Mar-26']
          hdr0 = [hdr0, '#']
          hdr0 = [hdr0, '# Notes: response function including filter+CCD+telescope optics']
          hdr0 = [hdr0, '#']
          hdr0 = [hdr0, '# Downloaded from http://www.eso.org/sci/facilities/paranal/instruments/vimos/inst/imaging.html']
          hdr0 = [hdr0, '#']

;         ngood = n_elements(wave) & good = lindgen(ngood)
;         good = where((qresp[*,iq] gt 1E-4),ngood)
          good = where((qresp[*,iq] gt 0.0),ngood)
          kfilter1 = {lambda:0.0D, pass:0.0D}
          hdr = hdr0

          kfilter = replicate(kfilter1,ngood)
          kfilter.lambda = wave[good]
          kfilter.pass = qresp[good,iq]

          if keyword_set(debug) then begin
             plot, kfilter.lambda, kfilter.pass, xsty=3, ysty=3, $
               title=bandpass[ii]+'/'+quad[iq], charsize=2.0
             cc = get_kbrd(1)
          endif
             
          yanny_write, path+newfilter_root[ii]+'_'+quad[iq]+'.par', $
            ptr_new(kfilter), hdr=hdr, stnames=bandpass[ii]+'FILTER'

;         for jj = 0L, n_elements(wave)-1L do printf, lun, $
;           bandpass[ii]+'FILTER', wave[jj], q1[jj], $
;           format='(A7,F15.5,3x,F15.12)'
;         free_lun, lun

       endfor

; also produce *average* bandpass filter response functions

       qresp_avg = total(qresp,2)/float(n_elements(quad))
       qresp_med = djs_median(qresp,2)
       
       hdr0 = '# Units:'
       hdr0 = [hdr0, '#  "lambda" is in Angstroms']
       hdr0 = [hdr0, '#  "pass" is the contribution to the detector signal per photon']
       hdr0 = [hdr0, '#']
       hdr0 = [hdr0, '# Bandpass Name(s): '+bandpass[ii]]
       hdr0 = [hdr0, '#']
       hdr0 = [hdr0, '# Instrument: VIMOS']
       hdr0 = [hdr0, '#']
       hdr0 = [hdr0, '# Determined by: Mieske et al. 2007:']
       hdr0 = [hdr0, '#   http://www.eso.org/sci/facilities/paranal/instruments/vimos/inst/filters_2007-04-04/memo_filtrans_March26th.pdf']
       hdr0 = [hdr0, '#']
       hdr0 = [hdr0, '# Date of determination: 2007-Mar-26']
       hdr0 = [hdr0, '#']
       hdr0 = [hdr0, '# Notes: response function including filter+CCD+telescope optics']
       hdr0 = [hdr0, '#']
       hdr0 = [hdr0, '# Downloaded from http://www.eso.org/sci/facilities/paranal/instruments/vimos/inst/imaging.html']
       hdr0 = [hdr0, '#']

       good = where((qresp_avg gt 0.0),ngood)
       kfilter1 = {lambda:0.0D, pass:0.0D}
       hdr = hdr0

       kfilter = replicate(kfilter1,ngood)
       kfilter.lambda = wave[good]
       kfilter.pass = qresp_avg[good]

       yanny_write, path+newfilter_root[ii]+'.par', $
         ptr_new(kfilter), hdr=hdr, stnames=bandpass[ii]+'FILTER'

       if keyword_set(debug) then begin
          plot, wave, qresp_avg, yr=[0,max(qresp)], charsize=2.0, $
            title=bandpass[ii]
          oplot, wave, qresp_med, line=2, thick=2
          for iq = 0L, n_elements(quad)-1L do oplot, wave, qresp[*,iq], line=1
          cc = get_kbrd(1)
       endif
          
    endfor

stop    

return
end
    
