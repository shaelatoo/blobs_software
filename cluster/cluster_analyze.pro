pro cluster_analyze,file,angle,angwidth,nmean,rmin,rwidth, $
    speeds,t0s,edit1,display=display,savefile=savefile, $
    srad=srad

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;                                                            ;;
;; Purpose: Gets height-time image based on inputs, analyzes  ;;
;;          with radon transform, thresholds, closes, clus-   ;;
;;          ters, converts remaining tracks to physically     ;;
;;          meaningful parameters, and saves data to file.    ;;
;;          (Also sometimes performs circus tricks.)          ;;
;;                                                            ;;
;; Inputs: file - difference movie file from desired date/    ;;
;;                time                                        ;;
;;         angle - desired central angle of angular bin       ;;
;;         angwidth - desired width of angular bin            ;;
;;         nmean - multiple of mean value to threshold trans- ;;
;;                 form to                                    ;;
;;         rmin - minimum radius to include in bin (pixels)   ;;
;;         rwidth - width of radial range (in pixels)         ;;
;;                                                            ;;
;; Outputs: speeds,t0s - calculated speed, onset time for     ;;
;;            each detected event                             ;;
;;          edit1 - clustered transform backprojected for     ;;
;;            display                                         ;;
;;                                                            ;;
;; Keyword: display - displays backprojected of clustered     ;;
;;            result                                          ;;
;;          savefile - set to save analysis results to file   ;;
;;          srad - sets srad keyword when calling trans; uses ;;
;;            s_radon to perform radon transform              ;; 
;;                                                            ;;
;; Notes: Note able to deal with diff sizes other than 1024   ;;
;;                                                            ;;
;; Dependencies: trans,cluster_mask,axes,s_back,centroid      ;;
;;                                                            ;;
;; Created: 10/23/07                                          ;;
;;                                                            ;;
;; Created by: Shaela Jones                                   ;;
;;                                                            ;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

; parameters
outsz=1024
ratio_factor=50


; get transformed image
if keyword_set(srad) then begin
  trans,file,rmin,rwidth,angle,angwidth,nmean,orig_im,thetas,rs, $
         fun3=httmim,/srad
endif else begin
  trans,file,rmin,rwidth,angle,angwidth,nmean,orig_im,thetas,rs, $
         fun3=httmim
endelse
cluster_mask,httmim,mask


; initializations
szhttmim=SIZE(httmim)
nclusts=MAX(mask)
edit1=BYTARR(szhttmim[1],szhttmim[2])



; choose representative points for each cluster
for i=1,nclusts do begin
  locs=WHERE(mask eq i)
  if locs[0] eq -1 then begin
    print,'Problem with cluster_mask result'
    return
  endif
  ave=AVERAGE(httmim[locs])
  ind=ARRAY_INDICES(httmim,locs)
  minthet=MIN(ind[0,*],max=maxthet)
  minr=MIN((ind[1,*]),max=maxr)
  subarray=httmim[minthet:maxthet,minr:maxr]
  submask=mask[minthet:maxthet,minr:maxr]
  nthet=maxthet-minthet+1
  locs=WHERE(submask ne i)
  if locs[0] ne -1 then subarray[locs]=0.
;; weight subarray by slope value
;  for j=0,nthet-1 do subarray[j,*]=subarray[j,*]
  szsub=SIZE(subarray)
  if szsub[0] gt 1 then cent=CENTROID(subarray) else cent=[0,0]
  if szsub[0] gt 1 then maxpos=ARRAY_INDICES(subarray, $
      WHERE(subarray eq MAX(subarray))) else maxpos=[0,0]
  edit1[maxpos[0]+minthet,maxpos[1]+minr]=ave
endfor
;stop



; display point choice results
if keyword_set(display) then begin
  restore,file
  diffs=0
  axes,dttms,rmin,rwidth,filelist[0],ctr,angle,outsz,tms,alts
  szorig=SIZE(orig_im)
  if keyword_set(srad) then begin
    crazy1=s_back(edit1,thetas,rs,tms,alts,/bilinear)
  endif else begin
    crazy1=radon(edit2,/double,/backproject, $
             theta=thetas,rho=rs,nx=szorig[1],ny=szorig[2])
  endelse
;  rat=MEAN(orig_im[where(orig_im ne 0)])/MEAN(crazy1[where(crazy1 ne 0)])
  utplot_image,SQRT(SQRT(crazy1))*ratio_factor+orig_im,dttms,alts, $
          title='Enhanced Features Laid Over Original Image', $
          ytitle='Altitude (Solar Radii)'
  stop
endif



; analyze
if NOT(KEYWORD_SET(display)) then begin
  restore,file
  diffs=0
  first=filelist[0]
  filelist=0
  axes,dttms,rmin,rwidth,first,ctr,angle,outsz,tms,alts
  dttms=0
endif
if KEYWORD_SET(srad) then rs=REVERSE(rs)
ntms=N_ELEMENTS(tms)
nalts=N_ELEMENTS(alts)
;deltrs=alts[1:*]-alts[0:nalts-2]
;deltr=AVERAGE(deltrs)
;deltts=tms[1:*]-tms[0:ntms-2]
;deltt=AVERAGE(deltts)
maxx=MAX(tms,min=minx)
delty=MAX(alts)-MIN(alts)
deltx=maxx-minx
rexs=(tms-minx)/deltx-0.5
inds=ARRAY_INDICES(edit1,WHERE(edit1 ne 0))
ct=COS(thetas[inds[0,*]])
st=SIN(thetas[inds[0,*]])
speeds=-delty/deltx*ct/st
spt=(0.5-st)/ct
xsp=-ct*rs[inds[1,*]]+st*spt
num=N_ELEMENTS(WHERE(edit1 ne 0))
t0s=FLTARR(num)
for i=0,num-1 do begin
;  if xsp gt -0.5 then begin
;    garbage=MIN(ABS(xsp[i]-rexs),before)
;    if xsp[i]-rexs[before] lt 0 then before--
;    xind=before+(xsp[i]-rexs[before])/(rexs[before+1]-rexs[before])
;    t0s[i]=INTERPOLATE(tms,xind)
;  endif else begin
    t0s[i]=(xsp[i]+0.5)*deltx+tms[0]
endfor
;sps=(nalts/2.+st*rs[inds[1,*]])/ct
;xsps=+rs[inds[1,*]]*ct+sps*st
;t0s=dttms[xsps+ntms/2.-1]



; file analysis results

if keyword_set(savefile) then begin
  date=STRMID(dttms[0],0,4)+STRMID(dttms[0],5,2)+STRMID(dttms[0],8,2)
  moving=0
  place=STRPOS(file,'moving',27)
  if place eq -1 then place=STRPOS(file,'running',27) else moving=1
  if place eq -1 then place=STRPOS(file,'base',27)
  dir='/home/shaela/stereo/radon/'
  flnm=dir+date+spcrft+STRING(angle)+'_'+ $
      STRMID(file,place-1)
  flnm=STRCOMPRESS(flnm,/remove_all)
  save,filename=flnm,file,angle,angwidth,nmean,rmin, $
    rwidth,speeds,t0s,edit1,srad,httmim,thetas,rs
;  flnm=dir+date+spcrft+STRCOMPRESS(angle)+'_'+ $
;      STRMID(file,place-1)
;  flnm=STRMID(flnm,0,STRLEN(flnm)-4)+'.csv'
;  openw,unit,flnm,/get_lun
;  printf,unit,'cluster_analyze',file,angle,angwidth,nmean,rmin,rwidth
;  printf,unit,' '
;  printf,unit,'Event Number',',','t0',',','speed'
;  for j=1,nclusts do printf,unit,j,t0s[j-1],speeds[j-1]*696000., $
;          format='(I14,",",F14.3,",",F14.1)' 
;  free_lun,unit
endif



;stop
end





; obsolete

;print,'Event Number','t0','speed',format='(3A14)'
;for i=1,nclusts do print,i,t0s[i-1],speeds[i-1]*696000., $
;          format='(I14,F14.3,F14.1)'


;  crazy2=radon(edit2,/double,/backproject, $
;             theta=thetas,rho=rs,nx=szorig[1],ny=szorig[2])
;  crazy3=radon(edit3,/double,/backproject, $
;             theta=thetas,rho=rs,nx=szorig[1],ny=szorig[2])
;  window,10,retain=2,xsize=9*szorig[1],ysize=3*szorig[2]
;  tvscl,enlarge_image2(orig_im,9,3)  


;wait,1.
;  tvscl,enlarge_image2(crazy2,3,3)
;  wait,1.
;  tvscl,enlarge_image2(crazy3,3,3)
;  wait,1.

;  edit1[cent[0]+minthet,cent[1]+minr]=ave
;  edit2[maxthet,rmean]=ave
;  edit3[thetmean,rmean]=ave

;edit2=BYTARR(szhttmim[1],szhttmim[2])
;edit3=BYTARR(szhttmim[1],szhttmim[2])
