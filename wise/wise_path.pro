function wise_path, ppxf=ppxf
; jm11apr13ucsd
    
    wisepath = getenv('IM_RESEARCH_DIR')+'/projects/wise/'
    if keyword_set(ppxf) then wisepath = wisepath+'ppxf/'

return, wisepath
end
