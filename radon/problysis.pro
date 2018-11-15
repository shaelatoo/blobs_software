pro problysis,diffs_file,ang,angwidth,rmin,r_width,speeds, $
         t0s,nomspeeds,allgood=allgood

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;                                                            ;;
;; Purpose: Analyzes the diff images in diffs_file by forming ;;
;;          a height-time image, using prob_radon to calculate;;
;;          the probability of a mean statistically different ;;
;;          than 
;;                                                            ;;
;; Created: 01/23/08                                          ;;
;;                                                            ;;
;; Created by: Shaela Jones                                   ;;
;;                                                            ;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

; parameters
resz=4
outsz=1024
szstruc=5
mult=14
minprob=0.9
minmean=0.0
minlocs=0.15


; initializations
struc=[REPLICATE(1,szstruc)]



; restore diffs file
restore,diffs_file
firstfile=filelist[0]
diffs=0
filelist=0



; main
temp=' '
;r_strp,resz,rmin,r_width,ang,angwidth,diffs,dttms,ctr, $
;     firstfile,radpic=temp
;diffs=0
r_strp,resz,rmin,r_width,ang,angwidth,diffs_file=diffs_file, $
      radpic=temp
sztemp=SIZE(temp)
factor1=600./sztemp[1]
factor2=600./sztemp[2]
tempmean=MEAN(temp)
temp=temp-tempmean
axes,dttms,rmin,r_width,firstfile,ctr,ang,outsz,tms,alts
thetas=find_thetas(temp)
nthets=N_ELEMENTS(thetas)
rs=(FINDGEN(nthets)-nthets/2)/nthets
if sztemp[1] gt sztemp[2] then nrs=sztemp[1] else nrs=sztemp[2]
rs=rs*nrs
prob_radon,temp,tms,alts,thetas,rs,trans,rad
copytrans=trans
copyrad=rad
trans[where(trans lt minprob)]=0.0
rad[where(rad lt minmean)]=0.0
combined=trans*rad
cluster_mask,combined,mask
for i=0,MAX(mask)-1 do begin
  locs=WHERE(mask eq i)
  if N_ELEMENTS(locs) lt minlocs*nthets then combined[locs]=0
endfor
combined2=CONVOL(combined,struc,/edge_truncate)
cluster_mask,combined2,mask
nclusts=MAX(mask)
print,'nclusts= ',nclusts



; choose representative points for each cluster
edit1=rep2d(mask,combined)



;; backproject and display
;;edit2=REVERSE(edit1,2)
;;back=RADON(edit2,/double,/backproject,theta=thetas,rho=rs, $
;;      nx=sztemp[1],ny=sztemp[2])
;back=s_back(edit1,thetas,rs,tms,alts)
;window,retain=2,xsize=sztemp[1]*factor1,ysize=sztemp[2]* $
;         factor2,/free
;tvscl,enlarge_image2(temp+back*mult,factor1,factor2)
;;stop



; calculate speeds,t0s
ntms=N_ELEMENTS(tms)
nalts=N_ELEMENTS(alts)
deltr=AVERAGE(alts[1:*]-alts)
deltt=AVERAGE(tms[1:*]-tms)
inds=ARRAY_INDICES(edit1,WHERE(edit1 ne 0))
ct=COS(thetas[inds[0,*]])
st=SIN(thetas[inds[0,*]])
speeds=REFORM(-deltr/deltt*ct/st)
sps=(nalts/2.-st*rs[inds[1,*]])/ct
xsps=-rs[inds[1,*]]*ct+sps*st
timeinds=xsps+ntms/2.-1
;stop


; calculate nominal speeds
cthet=COS(thetas)
sthet=SIN(thetas)
nomspeeds=REFORM(-deltr/deltt*cthet/sthet)



; display and query
szedit=SIZE(edit1)
WINDOW,retain=2,xsize=sztemp[1]*factor1,ysize=sztemp[2]* $
         factor2,/free
winnum=!D.WINDOW
back=s_back(edit1,thetas,rs,tms,alts)
 tvscl,enlarge_image2(temp+back*1000,factor1,factor2)
;stop
nin=0
for i=0,N_ELEMENTS(timeinds)-1 do begin
  if timeinds[i] ge 0 then begin
    if NOT(KEYWORD_SET(allgood)) then begin
      new=FLTARR(szedit[1],szedit[2])
      new[inds[0,i],inds[1,i]]=2.
      backp=s_back(new,thetas,rs,tms,alts)
      WINDOW,retain=2,xsize=sztemp[1]*factor1,ysize=sztemp[2]* $
             factor2,/free
      winnum2=!D.WINDOW
      tvscl,enlarge_image2(temp+backp*1000,factor1,factor2)
      stop
      read,'real?',truth
      WDELETE,winnum2
    endif else truth=1
    if truth eq '1' then begin
      if N_ELEMENTS(newspeeds) eq 0 then begin
        newspeeds=speeds[i]
        newtimeinds=timeinds[i]
        newinds=inds[*,i]
        nin=1
      endif else begin
        newspeeds=[newspeeds,speeds[i]]
        newtimeinds=[newtimeinds,timeinds[i]]
        temp2=newinds
        newinds=INTARR(2,nin+1)
        newinds[*,0:nin-1]=temp2
        newinds[*,nin]=inds[*,i]
        nin++
      endelse
    endif else edit1[inds[0,i],inds[1,i]]=0.
  endif else edit1[inds[0,i],inds[1,i]]=0.
endfor
if N_ELEMENTS(newtimeinds) ne 0 then begin
  t0s=ANYTIM(INTERPOLATE(tms,newtimeinds),/ccsds)
;  t0s=dttms[newtimeinds]
  speeds=TEMPORARY(newspeeds)
  inds=TEMPORARY(newinds)
endif else t0s=-1



; backproject and display
;edit2=REVERSE(edit1,2)
;back=RADON(edit2,/double,/backproject,theta=thetas,rho=rs, $
;      nx=sztemp[1],ny=sztemp[2])
back=s_back(edit1,thetas,rs,tms,alts)
;window,retain=2,xsize=sztemp[1]*factor1,ysize=sztemp[2]* $
;         factor2,/free
tvscl,enlarge_image2(temp+back*mult,factor1,factor2)



;stop
WDELETE,winnum
end


; obsolete

;vals=DIST(szstruc)
;struc=vals
;struc[where(struc ge 1.)]=1.
;struc[where(struc lt 1.)]=0.
