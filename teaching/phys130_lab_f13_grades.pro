pro phys130_lab_f13_grades, alldata, test=test, sendit=sendit, final=final
; jm13aug25siena - parse the grades for this class 

    path = getenv('TEACHING_DIR')+'/Phys130/130-F13/grades/lab/'
    
    date = '13oct30' ; update this
    semester = 'Fall 2013'
    class = 'Physics 130 - General Physics I - Lab'
    
; read the grade spreadsheet downloaded from GoogleDocs
    gradefile = path+'phys130-lab-'+date+'.csv'
    data = read_gradefile(gradefile,unique_assignments=assign)

; specify the complete list of *possible* assignments and their
; relative weights
    allassign = 'Lab'
    weight = 1.0

;   data = data[where(strtrim(data.first_name,2) eq 'Victoria')]
    
    process_grades, data, assign=assign, allassign=allassign, $
      weight=weight, class=class, semester=semester, /lab, $
      test=test, sendit=sendit, alldata=alldata, final=final
    srt = sort(alldata.current_grade)
    struct_print, alldata[srt]

return
end
