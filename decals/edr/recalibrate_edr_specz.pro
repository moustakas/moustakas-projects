pro recalibrate_edr_specz
; jm14nov13siena - recalibrate the SDSS DR10 photometry structure
; written out by MATCH_EDR_SPECZ to the latest (DR13) photometric
; uber-calibration  

    topdir = '/project/projectdirs/cosmo/work/decam/release/edr/'
;   topdir = '/global/data/desi/decam/release/edr/'
    edrfile = '~/edr-specz-dr10-oldcalib.fits'
;   edrfile = topdir+'spAll-dr10-edr-oldcalib.fits'
    outfile = topdir+'edr-specz-dr10.fits'

; read the full catalog and then apply the recalibration in a loop
    splog, 'Reading '+edrfile
    tractor = mrdfits(edrfile,1)
    spec = mrdfits(edrfile,2)
    phot = mrdfits(edrfile,3)
    unwise = mrdfits(edrfile,4)
    nobj = n_elements(phot)

    outphot = phot
    for ii = 0L, nobj-1 do begin
       print, format='("Recalibrating ",I0,"/",I0, A10,$)', ii, nobj, string(13b)
       outphot1 = phot[ii]
       sdss_recalibrate, outphot1
       outphot[ii] = outphot1
    endfor

    splog, 'Writing '+outfile
    mwrfits, tractor, outfile, /create
    mwrfits, spec, outfile
    mwrfits, outphot, outfile
    mwrfits, unwise, outfile
    
stop    

return
end
    
