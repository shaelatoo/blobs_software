pro cluster_mask,trans,mask,minval,thetas,rs,szorig, $
        display=display

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;                                                            ;;
;; Purpose: Clusters contiguous groups of non-zero pixels in  ;;
;;          trans, creates a mask that numbers contiguous     ;;
;;          clusters sequentially.                            ;;
;;                                                            ;;
;; Inputs: trans - 2D array containing a radon transform to be;;
;;                masked                                      ;;
;;                                                            ;;
;; Optional Inputs: minval - the starting point for numbering ;;
;;                    of clusters                             ;; 
;;                  thetas, rs - radon transform parameters   ;;
;;                    for use with display keyword            ;;
;;                  szorig - size of image from which trans   ;;
;;                    was created (for use with display kywd) ;;
;;                                                            ;;
;; Outputs: mask - array of same dimensions as trans; entries ;;
;;                 contain the number of the contiguous group ;;
;;                 the corresponding entry in trans belongs to;;
;;                                                            ;;
;; Dependencies: enlarge_image2                               ;;
;;                                                            ;;
;; Created: 10/18/07                                          ;;
;;                                                            ;;
;; Created by: Shaela Jones                                   ;;
;;                                                            ;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


; error checks
if keyword_set(display) then begin
  if (N_ELEMENTS(thetas) eq 0) OR (N_ELEMENTS(rs) eq 0) OR $
       (N_ELEMENTS(szorig) eq 0) then begin
    print,'Not all information needed for display provided.'
    return
  endif
endif



; initializations
if N_ELEMENTS(minval) eq 0 then minval=1
sztrans=SIZE(trans)
locs=WHERE(trans ne 0)
nclusts=N_ELEMENTS(locs)
list=ARRAY_INDICES(trans,locs)
mask=INTARR(sztrans[1],sztrans[2])
mask[list[0,*],list[1,*]]=LINDGEN(nclusts)+minval
oldnclusts=nclusts



; aggregate adjacent pixels
adjacent=[[0,1,0],[1,1,1],[0,1,0]]
repeat begin
  changed=0
  for i=minval+0L,MAX(mask) do begin
    locs=WHERE(mask eq i)
    if locs[0] ne -1 then begin
      szclust=N_ELEMENTS(locs)
      pos=ARRAY_INDICES(mask,locs)
      for j=0L,szclust-1 do begin
        x=pos[0,j]
        y=pos[1,j]
        mins=pos[*,j]-1
        maxes=pos[*,j]+1
        temp=adjacent
        if x eq 0 then begin
          temp=temp[1:2,*]
          mins[0]=0
        endif else if x eq sztrans[1]-1 then begin
          temp=temp[0:1,*]
          maxes[0]=sztrans[1]-1
        endif
        if y eq 0 then begin
          temp=temp[*,1:2]
          mins[1]=0
        endif else if y eq sztrans[2]-1 then begin
          temp=temp[*,0:1]
          maxes[1]=sztrans[2]-1
        endif
        subarray=temp*mask[mins[0]:maxes[0],mins[1]: $
                  maxes[1]]
        subnonzero=WHERE(subarray ne 0)
        different=WHERE(subarray[subnonzero] ne i,count)
        if count ne 0 then begin
          submin=MIN(subarray[subnonzero])
          inds=ARRAY_INDICES(subarray,subnonzero)
          mask[mins[0]+inds[0,*],mins[1]+inds[1,*]]=submin
          changed++
        endif
      endfor
    endif
  endfor
;  oldoldnclusts=oldnclusts
;  oldnclusts=nclusts
  uniques=UNIQ(mask,SORT(mask))
  nclusts=N_ELEMENTS(uniques)-1
;endrep until (((nclusts eq oldoldnclusts) and $
;           (nclusts eq oldnclusts)) or (nclusts eq 1))
endrep until (changed eq 0)



; re-number clusters
vals=mask[UNIQ(mask,SORT(mask))]
for i=1L,nclusts do begin
  locs=WHERE(mask eq vals[i])
  mask[locs]=i+minval-1
endfor
  


; display results
if keyword_set(display) then begin
  for i=1,nclusts do begin
    copy=FLTARR(sztrans[1],sztrans[2])
    locs=WHERE(mask eq i)
    copy[locs]=i
    TVSCL,copy
    wait,1.
    crazy=radon(copy,/double,/backproject,theta=thetas, $
             rho=rs,nx=szorig[1],ny=szorig[2])
    tvscl,enlarge_image2(crazy,3,3)
    wait,1.
  endfor
endif



end
           
