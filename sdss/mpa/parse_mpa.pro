;+
; NAME:
;   PARSE_MPA
;
; PURPOSE:
;   Parse the MPA files I care about.  Specifically, combine the O/H
;   and stellar mass files into one structure, and parse the
;   emission- and absorption-line measurements into iSPEC format. 
;
; COMMENTS:
;
; MODIFICATION HISTORY:
;   J. Moustakas, 2010 Mar 07, UCSD 
;-

pro parse_mpa, dr4=dr4

; I/O path and data file names; there are a few choices for the
; "indices" file: idxfix has been corrected for sky-line
; contamination; I don't know what "newidx" is, so I don't use it
; here; and "indx" is the table of raw measurements; for more details
; see http://www.mpa-garching.mpg.de/SDSS/DR7/raw_data.html#errors
    if keyword_set(dr4) then begin
       datapath = sdss_path(/mpa_dr4)
       suffix = 'dr4_v5_1b'
    endif else begin
       datapath = sdss_path(/mpa_dr7)
       suffix = 'dr7_v5_2'
    endelse
    
    infofile = datapath+'gal_info_'+suffix+'.fit.gz'
    linefile = datapath+'gal_line_'+suffix+'.fit.gz'
    indxfile = datapath+'gal_indx_'+suffix+'.fit.gz'
;   indxfile = datapath+'gal_idxfix_'+suffix+'.fit.gz'
    ohfile = datapath+'gal_fiboh_'+suffix+'.fits.gz'
    massfile = datapath+'totlgm_'+suffix+'.fit.gz'
    sfrfile = datapath+'gal_totsfr_'+suffix+'.fits.gz'

; output file names    
    ispeclinefile = datapath+'ispecline_'+suffix+'.fits'
    massohfile = datapath+'mpamassoh_'+suffix+'.fits'

; read the info file    
    splog, 'Reading '+infofile
    mpainfo = mrdfits(infofile,1)
    moretags = struct_trimtags(temporary(mpainfo),$
      select=['plateid','mjd','fiberid','ra','dec','z'])

; parse the stellar mass, metallicity, and SFR files
    splog, 'Reading '+massfile
    mpamass = mrdfits(massfile,1)
    splog, 'Reading '+ohfile
    mpaoh = mrdfits(ohfile,1)
    splog, 'Reading '+sfrfile
    mpasfr = mrdfits(sfrfile,1)

    ohtags = tag_names(mpaoh)
    sfrtags = tag_names(mpasfr)
    masstags = tag_names(mpamass)
    massoh = struct_addtags(im_struct_trimtags(mpamass,select=masstags,$
      newtags='mass_'+masstags),im_struct_trimtags(mpaoh,select=ohtags,$
      newtags='oh_'+ohtags))
    massoh = struct_addtags(temporary(massoh),im_struct_trimtags(mpasfr,$
      select=sfrtags,newtags='sfr_'+sfrtags))
    
    massoh = struct_addtags(moretags,temporary(massoh))
    im_mwrfits, temporary(massoh), massohfile, /clobber
    
; parse the emission- and absorption-line files into ispec format 
    splog, 'Reading '+linefile ; this file is big!
    line = mrdfits(linefile,1)
    splog, 'Reading '+indxfile ; this file is big!
    indx = mrdfits(indxfile,1)
    ispecline = parse_mpa_linefit(line,galindx=indx)
    ispecline = struct_addtags(moretags,temporary(ispecline))

    im_mwrfits, temporary(ispecline), ispeclinefile, /clobber

return
end
