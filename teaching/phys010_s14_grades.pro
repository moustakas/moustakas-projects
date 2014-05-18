pro phys010_s14_grades, alldata, test=test, sendit=sendit, final=final
; jm14feb02siena - parse the grades for this class 

    path = getenv('TEACHING_DIR')+'/010-S14/grades/'
    
    date = '14may13' ; update this
    semester = 'Spring 2014'
    class = 'Physics 010 - Introductory Astronomy'

; read the grade spreadsheet downloaded from GoogleDocs
    gradefile = path+'phys010_grades_'+date+'.csv'
    data = read_gradefile(gradefile,unique_assignments=assign)

; specify the complete list of *possible* assignments and their
; relative weights
    allassign = ['Homework','Activities','Position Paper','Midterm 1',$
      'Midterm 2','Midterm 3','Final Exam']
    weight = [0.15,0.15,0.2,0.10,0.10,0.10,0.2]
;   keep = where(strmatch(data.last_name,'*Ippo*'))
;   keep = where(data.final_exam gt 0.0)
;   data = data[keep]
    
    process_grades, data, assign=assign, allassign=allassign, $
      weight=weight, class=class, semester=semester, test=test, $
      sendit=sendit, alldata=alldata, final=final
    srt = sort(alldata.current_grade)
    struct_print, alldata[srt]

return
end
