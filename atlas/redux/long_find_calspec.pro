function long_find_calspec, allplan, datapath
; jm09dec17ucsd - given a PLAN, return the indices of CALSPEC
; standards 

    nall = n_elements(allplan)
    newplan = allplan

; grab the coordinates of the standards    
    std = where((strtrim(allplan.flavor,2) eq 'std'),nstd)
    if (nstd eq 0) then return, newplan
    plan = allplan[std]

    info = iforage(datapath+'Raw/'+strtrim(plan.filename,2))
;   struct_print, struct_trimtags(info,sel=['object','ra','dec'])
    
; read the calspec database info file
    allstd = rsex(getenv('XIDL_DIR')+'/Spec/Longslit/'+$
      'calib/standards/calspec/calspec_info.txt')

; spherematch generously    
    spherematch, 15.0*im_hms2dec(allstd.ra), im_hms2dec(allstd.dec), $
      15.0*im_hms2dec(info.ra), im_hms2dec(info.dec), 1.0, $
      m1, m2, max=0
;   indx[std[m2]] = 1

;   struct_print, allstd[m1]
;   struct_print, struct_trimtags(info[m2],$
;     sel=['object','ra','dec'])

; expand the plan with this new information
    starinfo = replicate({starfile: '...', starname: '...'},nall)
    starinfo[std[m2]].starfile = repstr(allstd[m1].file,'.fits.gz','')
    starinfo[std[m2]].starname = allstd[m1].star

    newplan = struct_addtags(newplan,starinfo)
;   struct_print, newplan

return, newplan
end
    
