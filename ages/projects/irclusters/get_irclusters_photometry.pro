pro get_irclusters_photometry, phot
; jm10nov08ucsd - retrieve the multiband photometry for
; Brodwin's IR-selected clusters project

    irpath = ages_path(/projects)+'irclusters/'
    datapath = getenv('DATA_DIR')+'/data/bootes/2010b/'
    ss = rsex(irpath+'FullSample.cat')
    iband = mrdfits(datapath+'bootes_I.fits.gz',1)
    good = where(iband.flag_duplicate eq 0)
    iband = iband[good]
    
    spherematch, iband.alpha_j2000, iband.delta_j2000, $
      ss.ra, ss.dec, 5.0/3600.0, m1, m2
;   m1 = im_spherematch(iband,ss,ratagname1='alpha_j2000',dectagname1='delta_j2000',$
;     match2=m2,radius=5.0,raoff=raoff,decoff=decoff)
    splog, raoff*3600, decoff*3600
    srt = sort(m2) & m1 = m1[srt] & m2 = m2[srt]
    help, ss, m1

    missing = lindgen(n_elements(ss))
    remove, m2, missing
    im_mwrfits, ss[missing], irpath+'missing.fits', /clobber

stop    
    
    
    bands = ['Bw','R','I','J','H','Ks','ch1','ch2','ch3','ch4']
    for jj = 0, n_elements(bands)-1 do begin
       phot1 = mrdfits(datapath+'bootes_'+bands[jj]+$
         '.fits.gz',rows=m1,1)
       tags = tag_names(phot1)
       phot1 = im_struct_trimtags(phot1,select=tags,$
         newtags=bands[jj]+'_'+tags)
       if (jj eq 0) then phot = phot1 else phot = $
         struct_addtags(temporary(phot),phot1)
    endfor

    im_mwrfits, phot, irpath+'irclusters_photometry.fits', /clobber
    
return
end
