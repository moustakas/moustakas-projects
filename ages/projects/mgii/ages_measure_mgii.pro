function ages_ismgii, specdata

    ismgii = [$
      '102/229',$ ; broad Mg II (z=0.426)
;     '107/019',$ ; *possible*, but too close to the edge
      '107/236',$ ; broad H-beta (Mg II on the edge) (z=0.354)
      '111/029',$ ; 4800 km/s (z=0.677), broad Mg II
      '114/004',$ ; gorgeous broad-line AGN!
      '202/119',$ ; narrow!
      '205/300',$ ; narrow!
      '207/254',$ ; narrow!
      '211/078',$ ; broad Mg II
      '302/141',$ ; broad Mg II
      '302/209',$ ; broad Mg II
      '302/283',$ ; broad Mg II
      '303/008',$ ; broad Mg II
      '303/044',$ ; broad Mg II
      '303/270',$ ; broad Mg II (z=0.662)
      '304/163',$ ; narrow; marginal detection
      '308/270',$ ; broad Mg II (but low S/N)
      '309/257',$ 
      '312/132',$ ; broad Mg II
      '312/227',$ ; broad Mg II
      '313/173',$ ; broad Mg II
      '314/008',$ ; broad Mg II
      '314/030',$ ; broad Mg II
      '314/284',$ ; broad Mg II, H-beta
      '315/102',$ ; broad Mg II
      '315/133',$ ; broad Mg II
      '315/149',$ ; broad Mg II
      '403/236',$ ; broad Mg II
      '404/045',$ ; broad Mg II
      '404/187',$ ; broad Mg II
      '408/093',$ ; broad Mg II
;     '409/191',$ ; broad Mg II(?), but low S/N; reject
      '416/275',$ ; broad Mg II, H-beta
      '417/105',$ ; broad Mg II (z=0.762)
      '418/060',$ ; broad Mg II
      '418/065',$ ; marginal
      '418/154',$ ; broad(?) Mg II (significant?) (z=0.785)
      '419/213',$ ; broad Mg II, but low S/N
      '422/083',$ ; broad Mg II
      '422/242',$ ; broad Mg II
      '502/242',$ ; broad Mg II, but low S/N
      '503/104',$ ; broad Mg II, but low S/N
      '506/295',$ ; broad Mg II, H-beta
;     '507/133',$ ; marginal; reject
      '508/077',$ ; broad H-beta and Mg II (z=0.426)
      '508/261',$ ; broad Mg II
      '510/139',$ ; narrow! low S/N
      '510/207',$ ; broad Mg II
      '510/281',$ ; broad Mg II
;     '525/278',$ ; possibly there, but *right* on the edge
      '526/290',$ ; broad Mg II
      '602/061',$ ; narrow!
      '602/238',$ ; broad Mg II
      '602/296',$ ; broad Mg II
;     '603/120',$ ; narrow; too low S/N
      '603/213',$ ; broad Mg II
      '606/245',$ ; broad Mg II
      '606/251',$ ; broad Mg II
      '607/220',$ ; broad Mg II
      '608/036',$ ; broad Mg II
      '608/253',$ ; broad Mg II
      '611/167',$ ; broad Mg II
      '709/281']  ; broad Mg II
;     '714/254',$ ; narrow; too low S/N

    all = string(specdata.pass,format='(I3.3)')+'/'+$
      string(specdata.aper,format='(I3.3)')
    indx = where_array(ismgii,all)
;   if (indx[0] ne -1) then out = specdata[indx] else out = -1
    
return, indx
end

pro ages_measure_mgii, qaplot_all=qaplot_all, debug=debug
; jm09dec17ucsd - measure Mg II in the AGES spectra

    common measure_mgii, specfit
    
    path = ages_path(/projects)+'mgii/'
    ss = read_ages_gandalf(/ppxf)
    aa = read_ages_gandalf(/ancillary)

; build a GANDALF/QAplot of all the candidate galaxies with Mg II
; emission; original list based on v1.0 PPXF/GANDALF fitting:
    if keyword_set(qaplot_all) then begin
       all = where(ss.mgii_2800[1] gt 0.0,nall)
       data = ss[all]
       ancillary = aa[all]
       allspecfit = read_ages_gandalf_specfit(data)
       qaplot_ages_gandalf_specfit, data, allspecfit, $
         psfile=path+'qaplot_mgii_all.ps'
       return
    endif

; having inspected the QAplot and removed crap measurements, go back
; through and measure the EWs and line-strengths more carefully using
; my splot code; write out the results
    bigindx = ages_ismgii(ss)
    data = ss[bigindx]
    ancillary = aa[bigindx]
    if (n_elements(specfit) eq 0) then specfit = $
      read_ages_gandalf_specfit(data)
    nobj = n_elements(data)

    linewave = 2799.495
    linename = 'MgII'
    nmonte = 100
    boxwidth = 50.0 ; [A]
    nterms = 5
    
    if keyword_set(debug) then qafile = path+'qaplot_mgii.ps'
    im_plotconfig, 0, pos, psfile=qafile
    for ii = 0, nobj-1 do begin
       splog, data[ii].pass, data[ii].aper
       good = where(specfit[ii].wave gt 0.0,npix)
       wave = exp(specfit[ii].wave[good])
       flux = specfit[ii].flux[good]
       ferr = specfit[ii].ferr[good]
       ivar = (1.0/ferr^2)*(ferr lt 1E5)
       
       res1 = im_splotew(wave,flux,ivar,linewave,nmonte=nmonte,$
         linename=linename,doplot=debug,/silent,boxwidth=25)
;      if keyword_set(debug) then cc = get_kbrd(1)
       if (ii eq 0) then res = res1 else res = [res,res1]
    endfor
    if keyword_set(debug) then im_plotconfig, $
      psfile=qafile, /psclose, /gzip
;   spawn, 'rsync -auv '+qafile+'.gz ~/', /sh

; now, having inspected ages_mgii.ps, restrict the list to the
; following spectra, as the others are too low S/N or spurious; INDX
; is the index list of good spectra (or pages of qaplot_mgii.ps.gz)
    indx = [0,2,5,7,8,9,10,12,13,17,$
      19,20,22,23,24,25,26,28,29,30,$
      36,37,40,44,45,46,48,50,51,52,$
      53,54,55,56,57]
    
; write out
    out = struct_trimtags(data[indx],select=['ages_id','pass',$
      'aper','z','zabs','vdisp','lumage','d4000_narrow','oii_3727_ew'])
    out = struct_addtags(out,struct_trimtags(ancillary[indx],select='ugriz_absmag'))
    out = struct_addtags(out,res[indx])
    
    outfile = path+'ages_mgii.fits'
    im_mwrfits, out, outfile, /clobber

    main = where((ss.oii_3727_ew[0]/ss.oii_3727_ew[1] gt 3.0) and $
      (ss.d4000_narrow[0]/ss.d4000_narrow[1] gt 3.0) and $
      (aa.main_flag eq 1))
    djs_plot, ss[main].d4000_narrow[0], ss[main].oii_3727_ew[0], $
      psym=3, xsty=3, ysty=3, /ylog
    djs_oplot, out.d4000_narrow[0], out.oii_3727_ew[0], $
      psym=6, color='red'

stop    
    
return
end
    
