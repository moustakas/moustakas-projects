function mips_path
; jm07jan31nyu
    
    path = '/Volumes/WDexternal/data/sg1120/mips/'
    if (file_test(path,/directory) eq 0L) then begin
       splog, 'External volume '+path+' not found.'
       return, -1L
    endif

return, path    
end
