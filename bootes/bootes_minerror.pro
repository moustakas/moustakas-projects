function bootes_minerror
; jm10jan23ucsd - minimum photometric uncertainty on the BOOTES
;   bandpasses 
    minerr = [$
      0.05,$           ; U
      0.02,0.02,0.02,$ ; BwRI
      0.02,$           ; z
      0.02,$           ; y
      0.05,0.05,0.05,$ ; JHKs
      0.1,0.1,0.1,0.1] ; ch1-4
return, minerr
end
