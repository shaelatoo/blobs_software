pro r_strp,resz,rmin,r_width,ang,angwidth,diffs,dttms,ctr, $
        first,radpic=radpic,img=img,diffs_file=diffs_file, $
        interp_r=interp_r,take_data=take_data,hts=hts, $
        no_adjcntrst=no_adjcntrst,nodiv=nodiv,de_trend= $
        de_trend,match_sigs=match_sigs

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;                                                            ;;
;; Purpose: Makes height-time diagrams from a series of diff- ;;
;;          erence images.  In each image a region specified  ;;
;;          by the user is integrated over angle and the res- ;;
;;          ultant pixel strips are combined sequentially     ;;
;;          to form one image.  The image can be used to make ;;
;;          visual height-time measurements of a feature and  ;;
;;          the height-time image is optionally returned.     ;;
;;                                                            ;;
;; Inputs: resz - factor by which to increase size            ;;
;;         rmin,r_width - minimum radius, radial width of ROI ;;
;;         ang,angwidth - center,width angle of bin (degrees) ;;
;;                                                            ;;
;; Optional inputs: diffs - difference images                 ;;
;;                  dttms - times of images provided in diffs ;;
;;                  ctr - sun center location in first diffs  ;;
;;                        image                               ;;
;;                  first - filename of the source file used  ;;
;;                          to create the first image         ;;
;;                                                            ;;
;; Keywords: radpic - variable in which to store real-size    ;;
;;                    height-time image                       ;;
;;           img - if set a height-time image is displayed    ;;
;;           diffs_file - file in which optional inputs are   ;;
;;                        stored; ignored if all inputs are   ;;
;;                        provided                            ;;
;;           interp_r - number of radial points to extract    ;;
;;                      data for, regardless of r_width       ;;
;;           take_data - if set the user is prompted to manu- ;;
;;                       ally select height-time points from  ;;
;;                       image to be fit                      ;;
;;           hts - 3xn array containing (times,heights,herr-  ;;
;;                 ors) of user-selected points; ignored if   ;;
;;                 take_data not set                          ;;
;;           no_adjcntrst - if set program does not adjust    ;;
;;                 contrast of height-time image              ;;
;;           nodiv - if set program does not divide each comp-;;
;;                 osite picture by the number of pixels summ-;;
;;                 ed within it (useful for mass calculation) ;;
;;           de_trend - attempt to remove a radial trend from ;;
;;                 each time step                             ;;
;;           match_sigs - adjusts the range of values in each ;;
;;                 column to give similar standard devia-     ;;
;;                 tions; for use with cor2 data              ;;
;;                                                            ;;
;; Dependencies: adj_cntrst,enlarge_image2,axes,adj_pix       ;;
;;                                                            ;;
;; Created: 07/11/07                                          ;;
;;                                                            ;;
;; Modifications: 11/19/07 - modified to automatically call   ;;
;;                           adj_cntrst if radpic is set      ;;
;;                05/05/08 - removed fit calculations,disp-   ;;
;;                           lay; added hts keyword           ;;
;;                11/20/08 - added de_trend keyword           ;;
;;                                                            ;;
;; Created by: Shaela Jones                                   ;;
;;                                                            ;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


; parameters
errsz=7
nangmin=2
degree=5



; restore difference images
if keyword_set(diffs_file) then begin
  restore,diffs_file
  first=filelist[0]
endif



; input checks
szdiffs=SIZE(diffs)
if szdiffs[0] lt 3 then begin
  print,'Image array must contain multiple two-dimensional images.'
  hts=-1
  return
endif
if szdiffs[1] ne szdiffs[2] then begin
  print,'Images are not symmetric.'
  hts=-1
  return
endif
;if szdiffs[1] ne 1024 then begin
;  print,'Input images have wrong resolution'
;  hts=-1
;  return
;endif



;initializations
npixang=FIX(rmin*!DTOR*angwidth)
n=szdiffs[3]
outsz=szdiffs[1]
if keyword_set(interp_r) then proxy=interp_r else proxy=r_width
radial_picture=FLTARR(n,proxy)
radius=rmin+FINDGEN(proxy)*FLOAT(r_width)/proxy
rmax=MAX(radius)
npts=INTARR(proxy)
summ=FLTARR(n,proxy)



; integrate image over angle
angmin=FIXANG(ang-angwidth/2.)
angmax=FIXANG(ang+angwidth/2.)
polar=FLTARR(2,4)
polar[0,*]=[angmin,angmin,angmax,angmax]
polar[1,*]=[rmin,rmax,rmin,rmax]
rect=cv_coord(from_polar=polar,/degrees,/to_rect)
ymax=1.
xmax=1.
ymin=FLOOR(MIN(rect[1,*],max=ymax))
ymax=CEIL(ymax)
xmin=FLOOR(MIN(rect[0,*],max=xmax))
xmax=CEIL(xmax) 
for i=ymin,ymax do begin
  for j=xmin,xmax do begin
    pol=CV_COORD(/degrees,from_rect=[j,i],/to_polar)
    if (pol[0] ge angmin) and (pol[0] le angmax) then begin
      if (pol[1] ge rmin) and (pol[1] le rmax) then begin
        x=j+ctr[0]
        y=i+ctr[1]
        if (y ge 0) and (y lt outsz) then begin
          if (x ge 0) and (x lt outsz) then begin
            temp=MIN(ABS(pol[1]-radius),rind)
            if rind ge rmax then print,rind
            summ[*,rind]=summ[*,rind]+diffs[x,y,*]
            npts[rind]++
          endif
        endif
      endif
    endif else if (angmin gt angmax) and ((pol[0] le angmax) $
            or (pol[0] ge angmin)) then begin   ; assumes cv_coord angle is between -180,180 always
      if (pol[1] ge rmin) and (pol[1] lt rmax) then begin
        rind=ROUND(pol[1])-rmin
        npts[rind]++
        x=j+ctr[0]
        y=i+ctr[1]
        summ[*,rind]=summ[*,rind]+diffs[x,y,*]
      endif
    endif
  endfor
endfor
if NOT(KEYWORD_SET(nodiv)) then begin
  halt=WHERE(npts eq 0)
  if halt[0] ne -1 then begin
    if ang-180. gt 0.25 then begin
      print,'Trouble in r_strp.'
      radpic=0.
      stop
    endif else begin
      if N_ELEMENTS(halt) eq 1 then npts[halt]=1 else begin
        print,'Trouble in r_strp.'
        radpic=0.
        stop
      endelse
    endelse
  endif
  for t=0,n-1 do radial_picture[t,*]=summ[t,*]/npts
endif else begin
  radial_picture=summ
endelse
;stop



; adjust radial contrast
if NOT(KEYWORD_SET(no_adjcntrst)) then begin
  radpic=adj_cntrst(radial_picture)
endif else begin
  radpic=radial_picture
endelse
szradpic=SIZE(radpic)
if KEYWORD_SET(de_trend) then begin
  array=FINDGEN(szradpic[2])
  for i=0,szradpic[1]-1 do begin
    adj_pix,radpic[i,*],array,degree,adjust,fit
;    radpic[i,*]=radpic[i,*]-adjust+fit[0]
    radpic[i,*]=radpic[i,*]-adjust
  endfor
;stop
endif
if KEYWORD_SET(match_sigs) then begin
  sigs=FLTARR(szradpic[1])
  for i=0,szradpic[1]-1 do sigs[i]=STDDEV(radpic[i,*],/nan)
  mnsig=MEAN(sigs)
  for i=0,szradpic[1]-1 do radpic[i,*]=(radpic[i,*]- $
         MEAN(radpic[i,*],/nan))/sigs[i]*mnsig
endif
  



; resize result
if KEYWORD_SET(img) OR KEYWORD_SET(take_data) then begin
  factor=FIX(1500./n)
  if factor le 1.0 then factor=1.0
  bigger=enlarge_image2(radial_picture,factor,resz)
endif



; make height-time image
if n_elements(img) ne 0 then begin
  if keyword_set(interp_r) then begin
    axes,dttms,rmin,r_width,first,ctr,ang,outsz,tms,rs,/interp_r
  endif else axes,dttms,rmin,r_width,first,ctr,ang,outsz,tms,rs
  date=strmid(dttms[0],0,4)+strmid(dttms[0],5,2) $
        +strmid(dttms[0],8,2)
  time=strmid(dttms[0],11,2)+strmid(dttms[0],14,2) $
        +strmid(dttms[0],17,2)
  str='/home/shaela/stereo/images/httm'+date+$
        '_'+time+STRCOMPRESS(fix(ang),/remove_all)+'.gif'
  ttl='Height-Time Image: '+STRCOMPRESS(fix(ang),/remove_all)+' degrees'
  utplot_image,radial_picture,tms,rs,title=ttl, $
        ytitle='POS Distance from Sun Center (Solar Radii)'
endif



; let user select points off of image
if KEYWORD_SET(take_data) then begin
  marked=FLTARR(2,n)
  marks=FLTARR(2,N)
  count=0
  done=0
  windw=SIZE(bigger)
  window,1,xsize=windw[1],ysize=windw[2],xpos=500,retain=2
  factor2=windw[1]/n
  while done ne 1 do begin
    tvscl,bigger
    xyouts,0.1,0.9,'done',/normal
    if count ne 0 then for j=0,count-1 do xyouts,marks[0,j],$
                                       marks[1,j],'X',/normal
    print,'Click on location of feature.'
    cursor,x,y,/up,/normal
    if (x lt 0.15) and (y gt 0.85) then done=1
    if done ne 1 then begin
      if count ne 0 then begin
        already=WHERE(marked[0,*] eq FIX(x*n))
        if already[0] ne -1 then begin
          print,'A point has already been chosen for that time.  Please carefully',$
                 ' choose a point at a different time.'
          cursor,x,y,/up,/normal
          already=WHERE(marked[0,*] eq FIX(x*n))
          if already[0] ne -1 then begin
            print,'User error.'
            hts=-1
            return
          endif
        endif
      endif
      marks[*,count]=[x,y]
      marked[*,count]=[FIX(x*n),FIX(y*proxy)]    ;mult. b/c of normal coords.
      count++
    endif
  endwhile
  wdelete,1
  if count le 1 then begin
    print,'Insufficient number of points selected.'
    hts=-1
    return
  endif



;   transform into usable height-time data
  tms=utc2tai(dttms)
  marked=marked[*,0:count-1]
  times=tms[marked[0,*]]
  ind=SORT(times)
  temp=times
  times=temp[ind]
  temp=REFORM(marked[1,*])
  marked[1,*]=temp[ind]
  temp=marked[1,1:*]-marked[1,*]
  locs=where(temp le 0)
  if locs[0] ne -1 then begin          ; if there's a problem
    print,'User Error.  Feature was not monotonically rising.'
    hts=-1
    return
  endif
  x=FLTARR(2,count)
  x[0,*]=ctr[0]+radius[marked[1,*]]*COS(ang*!DTOR)
  x[1,*]=ctr[1]+radius[marked[1,*]]*SIN(ang*!DTOR)
  if N_ELEMENTS(wcshead) eq 0 then begin
    secchi_prep,first,hdr,im,/calfac_off,/calimg_off,/silent
    struct=FITSHEAD2WCS(hdr)
    convert=hdr.dsun_obs/206265.0/1000.0
  endif else begin
    struct=wcshead
    convert=struct.position.dsun_obs/206265.0/1000.0
  endelse 
  pos=WCS_GET_COORD(struct,x)*convert
  rs=REFORM(SQRT((pos[0,*])^2+(pos[1,*])^2))





;   assemble error bars
  x[0,*]=x[0,*]+errsz*COS(ang*!DTOR)
  x[1,*]=x[1,*]+errsz*SIN(ang*!DTOR)
  pos=WCS_GET_COORD(struct,x)*convert
  errsup=REFORM(SQRT((pos[0,*])^2+(pos[1,*])^2))
  x[0,*]=x[0,*]-2*errsz*COS(ang*!DTOR)
  x[1,*]=x[1,*]-2*errsz*SIN(ang*!DTOR)
  pos=WCS_GET_COORD(struct,x)*convert
  errsdn=REFORM(SQRT((pos[0,*])^2+(pos[1,*])^2))
  errors=(errsup-errsdn)/2.



; store results
  count=N_ELEMENTS(rs)
  hts=DBLARR(3,count)
  hts[0,*]=times
  hts[1,*]=rs
  hts[2,*]=errors


endif      ; if taking data was indicated



end



; obsolete

;    t0t=dttms[0]
;    date=strmid(t0t,0,4)+strmid(t0t,5,2)+strmid(t0t,8,2)
;    tm=strmid(t0t,11,2)+strmid(t0t,14,2)+strmid(t0t,17,2)
;    fname='/home/shaela/stereo/walts/'+ $
;            date+'_'+tm+strcompress(img,/remove_all)+spcrft+ $
;            strcompress(ang,/remove_all)+'.sav'
;    save,filename=fname,bigger,spcrft,rmin,r_width,ang,errsz, $
;            angwidth,nangmin,rmin


;angle=FINDGEN(npixang)
;angle=ang+angle*angwidth/FLOAT(npixang)-angwidth/2.
;angle=angle*!DTOR


;get information for image axes
;  tms=anytim(dttms)
;  secchi_prep,first,hdr,im,/calfac_off,/calimg_off,/silent
;  struct=FITSHEAD2WCS(hdr)
;  convert=hdr.dsun_obs/206265.0/1000.0/696000.  

;  ;make image axes
;  x=FLTARR(2,r_width)
;  x[0,*]=ctr[0]+radius*COS(ang*!DTOR)
;  x[1,*]=ctr[1]+radius*SIN(ang*!DTOR)
;  pos=WCS_GET_COORD(struct,x)*convert
;  rs=REFORM(SQRT((pos[0,*])^2+(pos[1,*])^2))


;interpolate image, sum over angles
;for t=0,n-1 do begin
;  for r=0,r_width-1 do begin
;    npixang=FIX(!dpi*radius[r]*(angwidth/360.0))
;;    angle=FINDGEN(npixang)
;;    angle=ang+angle*angwidth/FLOAT(npixang)-angwidth/2.
;;    angle=angle*!DTOR
;    x=ctr[0]+radius[r]*COS(angle)
;    y=ctr[1]+radius[r]*SIN(angle)
;    overx=WHERE(x gt sz[1])
;    overy=WHERE(y gt sz[2])
;    if (overx[0] ne -1) or (overy[0] ne -1) then begin
;      print,'Extracted region outside of image'
;      return
;    endif
;    radial_picture[t,r]=TOTAL(INTERPOLATE(diffs[*,*,t],x,y,cubic=-0.5))
;  endfor
;endfor

;  sav=' '
;  read,'Would you like to save this data?',sav
;  if sav eq 'y' then begin
;    date=strmid(t0t,0,4)+strmid(t0t,5,2)+strmid(t0t,8,2)
;    time=strmid(t0t,11,2)+strmid(t0t,14,2)+strmid(t0t,17,2)
;    sz=SIZE(diffs)
;  ;  fname='C:\Users\Shaela\Documents\Home\stereo\walts\radial\Analyzed2\'+$
;    fname='/home/shaela/stereo/walts/radial/faintfeatures/'+ $
;           date+'_'+time+spcrft+STRCOMPRESS(ang,/remove_all)+ $
;           STRCOMPRESS(sz[1],/remove_all)+'.sav'
;    save,filename=fname,t0t,times,rs,poly,lin,bigger,marked, $
;         errors,marks,siglin
;    print,'Data saved to: ',fname
;  endif

;  print,poly[1],sig[1],FORMAT= $
;         '("Initial Radial Speed: ",E11.2," +/- ",E11.2," solar radii/s")'

;if npixang lt nangmin then begin
;  print,'Angular width too small.'
;  return
;endif

;            rind=ROUND(pol[1])-rmin

;  factor=proxy/n

;rect=FLTARR(2,4)
;rect[0,*]=polar[1,*]*COS(polar[0,*])
;rect[1,*]=polar[1,*]*SIN(polar[0,*])

;  t0=times[0]
;  times=times-t0
