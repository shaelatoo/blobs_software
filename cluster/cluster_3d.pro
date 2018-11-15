pro cluster_3d,trans,thetas,rhos,tms,alts,as,nmean,speeds,t0s, $
      accels,edit1,display=display,savefile=savefile,srad=srad, $
      useaves=useaves

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;                                                            ;;
;; Purpose: Thresholds 3D transforms, backprojects and (opt-  ;;
;;          ionally) displays the result, converts parameters ;;
;;          of surviving cluster points to physically mean-   ;;
;;          ingful quantities.                                ;;
;;                                                            ;;
;; Inputs: trans - 3D transform                               ;;
;;         thetas,rhos - theta and rho values corresponding   ;;
;;           to dimensions 1,2 of transform in normalized     ;;
;;           coordinate system                                ;;
;;         tms,alts - time, height values corresponding to    ;;
;;           dimensions 1,2 of original height-time image     ;;
;;         as - acceleration values correspoding to dimension ;;
;;           3 of the transform space in normalized coords.   ;;
;;         nmean - multiple of mean value at which to thresh- ;;
;;           old transform                                    ;;
;;                                                            ;;
;; Outputs: speeds,t0s,accels - calculated speed, onset time, ;;
;;            and acceleration for each detected event        ;;
;;                                                            ;;
;; Optional Outputs: edit1 - clustered transform backproject- ;;
;;            ed for display                                  ;;
;;                                                            ;;
;; Keyword: display - displays backprojection of clustered    ;;
;;            result                                          ;;
;;          savefile - set to save analysis results to file   ;;
;;          useaves - uses average value of transform within  ;;
;;            cluster to determine intensity of representat-  ;;
;;            ive point instead of weighting all equally      ;; 
;;                                                            ;;
;; Dependencies: cluster_mask,axes,s_back                     ;;
;;                                                            ;;
;; Created: 04/29/08                                          ;;
;;                                                            ;;
;; Created by: Shaela Jones                                   ;;
;;                                                            ;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

; parameters
ratio_factor=50
szstruc=3



; input checking
sztrans=SIZE(trans)
nthets=N_ELEMENTS(thetas)
nrho=N_ELEMENTS(rhos)
nas=N_ELEMENTS(as)
if (sztrans[1] ne nthets) or (sztrans[2] ne nrho) or (sztrans[3] $
         ne nas) then begin
  print,'Variables thetas, rhos, and as must have the same length ',$
         'as the first, second, and third dimensions of trans.'
  stop
endif



; initializations
maxx=MAX(tms,min=minx)
maxy=MAX(alts,min=miny)
rexs=(tms-minx)/(maxx-minx)-0.5
reys=(alts-miny)/(maxy-miny)-0.5



; threshold,mask transform
;threshold,trans,nmean,szstruc,fun7=thresh
stop
cluster_mask3D,trans,mask
nclusts=MAX(mask)
nx=N_ELEMENTS(tms)
ny=N_ELEMENTS(alts)
edit1=BYTARR(nthets,nrho,nas)



; choose representative points for each cluster
for i=1,nclusts do begin
  locs=WHERE(mask eq i)
  if locs[0] eq -1 then begin
    print,'Problem with cluster_mask result'
    return
  endif
  if KEYWORD_SET(useaves) then ave=AVERAGE(trans[locs]) $
           else ave=1.
  inds=ARRAY_INDICES(trans,locs)
  minthet=MIN(inds[0,*],max=maxthet)
  minrho=MIN(inds[1,*],max=maxrho)
  mina=MIN(inds[2,*],max=maxa)
  subarray=trans[minthet:maxthet,minrho:maxrho,mina:maxa]
  submask=mask[minthet:maxthet,minrho:maxrho,mina:maxa]
  nthet=maxthet-minthet+1
  locs=WHERE(submask ne i)
  if locs[0] ne -1 then subarray[locs]=0.
  szsub=SIZE(subarray)
  if szsub[0] gt 1 then cent=CENTROID(subarray) else cent=[0,0,0]
  if szsub[0] gt 1 then maxpos=ARRAY_INDICES(subarray, $
      WHERE(subarray eq MAX(subarray))) else maxpos=[0,0,0]
  edit1[maxpos[0]+minthet,maxpos[1]+minrho,maxpos[2]+mina]=ave
endfor



; display point choice results
if keyword_set(display) then begin
  szdisp=SIZE(display)
  dttms=TAI2UTC(tms,/ccsds)
  crazy1=s_back(edit1,thetas,rhos,tms,alts,as)
  utplot_image,crazy1*ratio_factor+display,dttms,alts, $
          title='Enhanced Features Laid Over Original Image', $
          ytitle='Altitude (Solar Radii)'
  stop
endif



; analyze
delty=maxy-miny
deltx=maxx-minx
paramset=ARRAY_INDICES(edit1,WHERE(edit1 ne 0))
ct=COS(thetas[paramset[0,*]])
st=SIN(thetas[paramset[0,*]])
speeds=-delty/deltx*ct/st*696000.
spt=(0.5-st)/ct
xsp=-ct*rhos[paramset[1,*]]+st*spt
num=N_ELEMENTS(paramset[0,*])
t0sfov=FLTARR(num)
for i=0,num-1 do begin
  exact_time=(xsp[i]+0.5)*deltx+tms[0]
  garbage=MIN(ABS(exact_time-tms),loc)
  if exact_time-tms[loc] lt 0. then loc++
  t0sfov[i]=tms[loc]
endfor
accels=2.*paramset[2,*]/(deltx^2/nx^2)*(delty/ny)*696000.



; file analysis results
if KEYWORD_SET(savefile) then begin
  save,filename=savefile,trans,thetas,rhos,as,tms,alts,nmean,speeds, $
    t0sfov,accels,paramset,edit1
endif



end


