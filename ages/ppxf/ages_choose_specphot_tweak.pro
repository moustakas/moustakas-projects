function smooth_tweak, rawwave, rawtweak, tweakwave, maxwave=maxwave
    inrange = where((rawwave gt -900.0) and (rawwave lt maxwave))
    tweak1 = interpol(rawtweak[inrange],rawwave[inrange],tweakwave) ; use interpol here!
    sset = bspline_iterfit(tweakwave,tweak1,$
      bkspace=150.0,yfit=tweak2)
    tweak = 10.0^(-0.4*tweak2) ; convert to linear
return, tweak ; note minus sign!
end

function ages_choose_specphot_tweak, pass, tweakwave=tweakwave
; jm09nov15ucsd - read the output from AGES_PPXF_SPECPHOT_TWEAK and
;   return the appropriate spectrophotometric tweak; the observed
;   spectrum must be *multiplied* by the output vector, after
;   interpolating onto the wavelength vector of the desired spectrum 

; note that five of the plates [106,110,209,310,311] were not fluxed
; at all, so the output from this routine for those plates is
; meaningless! 
    
    if (n_elements(pass) eq 0) then begin
       doc_library, 'ages_choose_specphot_tweak'
       return, -1
    endif

    version = ages_version(/ppxf_specfit)
    base_specfitpath = ages_path(/ppxf)
    infile = base_specfitpath+'ppxf_specphot_tweak_'+version+'.fits.gz'
    splog, 'Reading '+infile
    alltweak = mrdfits(infile,1)

    this = where(alltweak.pass eq pass,nthis)
    if (nthis eq 0) then message, 'No tweak found for pass '+$
      string(pass,format='(I0)')+'!'
    all = alltweak[this]

; the choice of which plates require tweaking was determined by
; examining the QAplot written out by AGES_PPXF_SPECPHOT_TWEAK
    good = where(all.wave gt -900.0)
    rawwave = all.wave
    rawtweak = all.tweak
        
    npix = ceil(((9200<max(all.wave[good]))-3700.0)/1.0)
    tweakwave = range(3700.0,(9200<max(all.wave[good])),npix)
;   tweakwave = im_array(3700.0,9200.0<max(all.wave[good]),1.0)
;   tweakwave = im_array(3700.0>min(all.wave[good]),$
;     9200.0<max(all.wave[good]),1.0)

; 104: the bump at 6800 is poor fluxing (telluric correction)
; 604: horrible!  never use!
    case strtrim(pass,2) of
       '104': tweak = smooth_tweak(rawwave,rawtweak,tweakwave,maxwave=6200.0) 
       '105': tweak = smooth_tweak(rawwave,rawtweak,tweakwave,maxwave=8000.0) 
       '113': tweak = smooth_tweak(rawwave,rawtweak,tweakwave,maxwave=8000.0) ; bad!
       '115': tweak = smooth_tweak(rawwave,rawtweak,tweakwave,maxwave=8000.0) ; bad!
       '201': tweak = smooth_tweak(rawwave,rawtweak,tweakwave,maxwave=8000.0) ; bad!
       '203': tweak = smooth_tweak(rawwave,rawtweak,tweakwave,maxwave=8000.0) ; bad!
       '204': tweak = smooth_tweak(rawwave,rawtweak,tweakwave,maxwave=8000.0) 
       '205': tweak = smooth_tweak(rawwave,rawtweak,tweakwave,maxwave=8000.0) ; bad!
       '206': tweak = smooth_tweak(rawwave,rawtweak,tweakwave,maxwave=8000.0) 
       '208': tweak = smooth_tweak(rawwave,rawtweak,tweakwave,maxwave=8000.0) ; bad!
       '212': tweak = smooth_tweak(rawwave,rawtweak,tweakwave,maxwave=8000.0) 
       '213': tweak = smooth_tweak(rawwave,rawtweak,tweakwave,maxwave=8000.0) ; bad!
       '214': tweak = smooth_tweak(rawwave,rawtweak,tweakwave,maxwave=8000.0)
       '301': tweak = smooth_tweak(rawwave,rawtweak,tweakwave,maxwave=8000.0) ; bad!
       '302': tweak = smooth_tweak(rawwave,rawtweak,tweakwave,maxwave=8000.0)
       '303': tweak = smooth_tweak(rawwave,rawtweak,tweakwave,maxwave=8000.0)
       '304': tweak = smooth_tweak(rawwave,rawtweak,tweakwave,maxwave=8000.0) ; bad!
       '305': tweak = smooth_tweak(rawwave,rawtweak,tweakwave,maxwave=8000.0)
       '306': tweak = smooth_tweak(rawwave,rawtweak,tweakwave,maxwave=8000.0)
       '307': tweak = smooth_tweak(rawwave,rawtweak,tweakwave,maxwave=8000.0) ; bad!
       '308': tweak = smooth_tweak(rawwave,rawtweak,tweakwave,maxwave=8000.0) ; bad!
       '309': tweak = smooth_tweak(rawwave,rawtweak,tweakwave,maxwave=8000.0) ; bad!
       '314': tweak = smooth_tweak(rawwave,rawtweak,tweakwave,maxwave=8000.0)
       '409': tweak = smooth_tweak(rawwave,rawtweak,tweakwave,maxwave=8000.0)
       '410': tweak = smooth_tweak(rawwave,rawtweak,tweakwave,maxwave=8000.0)
       '411': tweak = smooth_tweak(rawwave,rawtweak,tweakwave,maxwave=8000.0)
       '422': tweak = smooth_tweak(rawwave,rawtweak,tweakwave,maxwave=8000.0)
       else: tweak = tweakwave*0.0+1.0 ; no tweaking
    endcase

return, tweak
end
    
