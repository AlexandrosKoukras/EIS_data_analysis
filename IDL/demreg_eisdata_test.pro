;; Compare modeled intensities from DEM 
;; Testing two versions of the demreg method
; 
; ------------------------------------------------------------------------------------------------------------------------
; This example is focused one a signle pixel and for a particular EIS data set that contains 7 iron lines.
; These lines are: 
; fe 10 184.536
; fe 8 185.213
; fe 11 188.216
; fe 12 195.119
; fe 9 197.862
; fe 13 202.044
; fe 15 284.160
; 
; To save computing time and make this program more easy to share I have precalculated some of the necessary quantities (i.e. line fitting, contribution functions etc.)
; so I just restore their values. The corresponding sav files are shared together with this program. 
; 
; The contribution functions were computed using gofnt.pro from Chianti and with the following parameters
; density = 1e^8.5
; wvl_min = ion_wvl - 0.04 (Angstrom)
; wvl_max = ion_wvl + 0.04 (Angstrom)
; abund_file=concat_dir(concat_dir(!xuvtop,'abundance'),'sun_photospheric_1998_grevesse.abund')
; 
; Also, the intensities and the contribution functions are in energy units
; ------------------------------------------------------------------------------------------------------------------------

;; Load the intensities for a signle pixel 
;int_name = 'int_pix_list_ix20_iy30.sav'
;restore,int_name
;
;; Intensities for different pixels for additional tests (comment/uncomment below)
; ix=20, iy=30 ;QS
 int_pix_array = [112.92151057095538,37.79246486115742,116.98769207172637,123.94628769767091,20.79971852024349,52.49903295442659,14.218655749436259]
; ix=36, iy=208 ;CH
; int_pix_array = [13.740565632165717,9.037133976743078,23.741768296980567,26.260882299271206,3.461395990637454,18.19320810204372,11.322520988755027

; Load the contribution function
goft_name = 'g_i_list.sav'
restore,goft_name
g_i_array = g_i_list.ToArray()


;------------------------------------------------------
; Prepare inputs for DEM
;------------------------------------------------------

;; dn_in
dn_in = int_pix_array

;; edn_in
error_perc = 10.  ;manual selection for testing
edn_in = (error_perc * dn_in)/100.

;; temperature response matrix (n_tresp * nlines)
sz = size(g_i_array)
nf = sz[1]
nt = sz[2]
TRmatrix = fltarr(nt,nf)
for i=0,nf-1 do begin
  TRmatrix[*,i] = g_i_array[i,*]
endfor

;; tresp_logt : Temperature (log) binning of the contribution function
tt = findgen(31,increment=0.02,start=5.7) ; log(T): 5.7-6.3
tr_logt = tt

;; temps : The DEM 'target' temperatures (temperature bin edges)
tt2 = findgen(16,increment=0.04,start=5.7)
temps = 10^tt2


;------------------------------------------------------
; Call DEM (Hannah and Kontar 2013)
;------------------------------------------------------

dn2dem_pos_nb, dn_in, edn_in,TRmatrix,tr_logt,temps,dem1,edem1,elogt1,chisq1,dn_reg1,max_iter=200,/gloci


;------------------------------------------------------
; Call DEM (Hannah and Kontar 2012)
;------------------------------------------------------

min_logt = 5.7
max_logt = 6.3 
nt_bins = 25
order = 0
gloci = 1
pos = 1
channels = ['fe 10 184.536','fe 8 185.213','fe 11 188.216','fe 12 195.119','fe 9 197.862','fe 13 202.044','fe 15 284.160']

;bare minimum run
reg = data2dem_reg(tr_logt ,TRmatrix ,dn_in ,edn_in,channels=channels,gloci=gloci,pos=pos)

;modified run
;reg = data2dem_reg(tr_logt ,TRmatrix ,dn_in ,edn_in,mint=min_logt,maxt=max_logt,nt=nt_bins,channels=channels,gloci=gloci,order=order,pos=pos)

;------------------------------------------------------
; Compare the modeled and the observed intensities from both methods
;------------------------------------------------------

ratio1 = dn_reg1/dn_in 
ratio2 = reg.data_reg/dn_in 

; plot the intensities per line
window,1,xsize=1200
xind = [1,2,3,4,5,6,7]
plot,xind,ratio1,linestyle=0,xtickv=xind,xtickname=channels,xminor=1,xticks=7,xrange=[0,8],$
  psym=-2,xtitle='EIS lines',ytitle='Int_modeled /Int_observed' ;H&K 2013
oplot,xind,ratio2,linestyle=2,psym=-4 ;H&K 2012

; boxplot of the intensities for each method
methods_list = [' ','dn2dem_pos_nb','data2dem_reg',' '] 
bpd1 = createboxplotdata(ratio1)
bpd2 = createboxplotdata(ratio2)
boxes = boxplot([bpd1,bpd2],xminor=0,xmajor=4,xtickname=methods_list,ytitle='Int_modeled/Int_observed')



end
;------------------------------------------------------
; Plot the DEM from both methods
;------------------------------------------------------

;window,2
;plot,reg.logt,reg.dem,linestyle=2,/ylog

