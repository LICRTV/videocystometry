#pragma rtGlobals=3		// Use modern global access method and strict wave access.
﻿#pragma TextEncoding = "UTF-8"
#pragma rtGlobals=3		// Use modern global access method and strict wave access.
#include <Image Common>
#include <Image Threshold Panel>
#include <ImageSlider>
#include <ControlBarManagerProcs>

//VCAnalyser 1.1
//By Thomas Voets 20190806


menu "Videocystometry"

"VCAnalyser"
"SaveVCAnalysis"
"menuVCA"
"PressureFlowCurve"

end



macro menuVCA()

variable/G numpizzaslices = 12, pixelsize=51, length=50, bordershown=0, fillrate=30, framerate=1, numcores=1
make/O/N=1 xcenter=50, ycenter=50 

string/G loadedmovie = "no file loaded"



dowindow/K VCA
NewPanel/W=(100,50, 600,300)/N=VCA as "VCA"

Button LoadMovie title="Load Movie",size={100,20},pos={20,20},fColor=(65535,0,0), proc=loadmovie
Button Findcenter title="Find center",size={100,20},pos={130,20}, proc=findcenter
Button AnalyzeMovie title="Analyze Movie",size={100,20},pos={20,80},fColor=(65535,0,0), proc=startanalysis
setvariable numslices title = "Borderpoints",size={120,20}, pos={350,80}, value = numpizzaslices
setvariable pixsize title = "Pixel size",size={100,20}, pos={20,140}, value = pixelsize
setvariable filrate title = "Filling rate (µl/min)",size={150,20}, pos={140,140}, value = fillrate
setvariable frameratev title = "Frame rate (Hz)",size={140,20}, pos={310,140}, value = framerate, limits={0.1,100,0.01 }
setvariable loadedfile title ="Loaded file: ", size={250,20}, pos={20,50}, value = loadedmovie
PopupMenu Method title = "Method: ",pos={130,80}, size={200, 20}, value="Slope (fast);Fit (slow); Fit background subtracted (slow)"
setvariable multicore title = "How many cores?", pos={280,20},size={150, 20}, value=numcores, limits={1,ThreadProcessorCount,1}

setvariable xcent title = "Center x",size={100,20}, pos={20,110}, value = xcenter[0], proc=drawcircle
setvariable ycent title = "Center y",size={100,20}, pos={130,110}, value = ycenter[0], proc=drawcircle
setvariable diameter title = "Diameter",size={120,20}, pos={260,110}, value = length, proc=drawcircle

//setvariable borderdisplay4 title = "Borderpoint",size={170,20}, pos={150,190}, value = bordershown, proc=updatepizza


if (exists("movie2")==1)
loadedmovie = "Movie in memory"
button Loadmovie  fcolor=(0,10000,0), win=VCA
endif

end



macro loadmovie(test):buttoncontrol
string test

killwindow/Z bordermovement

killwaves/Z movie2

imageLoad/T=tiff/O/S=0/C=-1/Q/N=movie2
Button/Z AnalyzeMovie fColor=(65535,0,0)
Button/Z showresults disable =1

if (V_flag==1)
loadedmovie = S_filename
button Loadmovie  fcolor=(0,10000,0), win=VCA

endif
end

macro startanalysis(test):buttoncontrol
string test

variable hoofdtimer=STARTMStimer
variable/G bgintensity, usedframerate
controlinfo Method

pizza(numpizzaslices,pixelsize,(V_value-1))
button AnalyzeMovie  fcolor=(0,10000,0), win=VCA
bgintensity=calculatevolume(movie2,framerate, fillrate)
usedframerate = framerate

variable/G usedps=pixelsize, usedfr=fillrate, usedbg = bgintensity
duplicate/O newvolume, usednv
duplicate/O volume, usedv
Button showresults disable =0, title="Show results",size={100,20},pos={20,170},fcolor=(0,10000,0), win=VCA, proc=showresults

variable totaltime = stopMStimer(hoofdtimer)
print "Analysis took "+num2str(totaltime/1000000) +" seconds using "+ num2str(numcores) +" out of "+num2str(threadProcessorCount)+ " cores."
end




macro showresults(test):buttoncontrol
string test

string borderwave, borderFF

makesplineborder(bordermatrix,500)

variable/G slidervalue, bordershown2, firsttodelete=0, lasttodelete=0, maxspectrum=20, startsa=0, endsa=numpnts(borderpoint0), peethreshold=-3, spectb=0.1, spectw=0.05
make/O/N=1 indexlinex, indexlineyv, indexlineys, indexlineym

indexlinex=slidervalue
indexlineyv=wavemax(volume)
indexlineys=0
indexlineym=0


dowindow/K bordermovement

slidervalue=0

borderwave = "borderpointdif"+num2str(bordershown)
borderFF = "FFplot"+num2str(bordershown)

duplicate/O $borderFF, spectrum
duplicate/O $borderwave, bordermove
Smooth 5, bordermove
//bordermove *=framerate*pixelsize/1000
matrixop/O spectrum = layer(spectrum3d,bordershown)
setscale/P y 0,dimdelta(spectrum3D,1), "Hz", spectrum

newimage/N=bordermovement/G=1 spectrum
ModifyGraph axisEnab(left)={0.05,0.7},axisEnab(top)={0.05,0.4}
modifyGraph fSize(top)=12
ModifyImage spectrum ctab= {*,maxspectrum,ColdWarm,0}
SetAxis left *,*
duplicate/O bordermove framewave
framewave=p
AppendToGraph/L=topleft/T bordermove vs framewave
ModifyGraph axisEnab(topleft)={0.75,0.95}
SetAxis topleft -10,2
ModifyGraph freePos(topleft)=0;DelayUpdate


WMAppendAxisSlider()
SetVariable setvar3,pos={100,40},size={122,17},limits={-1,numpizzaslices,1},value = bordershown,title= "slice:", proc = updatepizza



Matrixop/O framem = movie2[][][slidervalue]

appendimage/R=righttop/T=topright  framem
ModifyGraph axisEnab(righttop)={0.5,1},axisEnab(topright)={0.6,1}
Matrixop/O framep = points[][][slidervalue]

AppendImage/R=righttop/T=topright framep
ModifyImage framep ctab={-1,1,Rainbow256,0}
makexy(newborder,slidervalue)
AppendToGraph/R=righttop/T=topright yvalues vs xvalues

make/O/N=200 spectrumband=0


controlbar 100
  slider sld, limits={0,(dimsize(movie2,2)-1),1}, ticks =0, pos={1200,150}, size={300,20},vert=0,variable = slidervalue,proc=ActionProcName
    SetVariable setvar0,pos={100,60},size={122,17},limits={0,(dimsize(movie2,2)-1),1},title = "Frame", value= slidervalue, proc=updateframe2
 SetVariable pizsize2,pos={300,60},title = "Pixel size (µm)",size={200,20},  value = pixelsize, proc=recalculatevolsur
 setvariable filrate2 title = "Filling rate (µl/min)",size={200,20}, pos={300,40}, value = fillrate, proc=recalculatevolint
  setvariable bgcorrection title = "Background",size={200,20}, pos={300,20}, value = bgintensity, proc=recalculatevolint
  slider spectrumscale vert=1, variable= maxspectrum, size={50,50}, pos={20,30}, limits={0, wavemax(spectrum), 0.5}, proc=changespectrumscale
  button replaceborder title= "Replace point by neighbours", size={200,20},pos={960,30}, proc = replaceborderpoint
  button restoreborder title= "Restore original point", size={200,20},pos={1180,30}, proc = restoreborderpoint
 button replaceseveral title= "Replace multiple points", size={200,20},pos={960,60}, proc = replaceborderpoints
 setvariable framerateadapt title = "Frame rate (Hz)",size={200,20}, pos={300,80}, value = framerate, limits = {0.5,100, 0.1}, proc = recalculatefr
setvariable from title = "from",size={100,20}, pos={1180,60}, value = firsttodelete,limits={0,(numpizzaslices-1),1}
setvariable to title = "to",size={100,20}, pos={1280,60}, value = lasttodelete,limits={1,(numpizzaslices-1),1}

button findvoids title= "Find voids for spectrum", size={200,20},pos={560,30}, proc=voidspectrum
button voidmovie title= "Make void movie", size={200,20},pos={560,80}, disable=2, proc=makevoidmovie
button bordermovie title= "Make bordermovie", size={200,20},pos={1020,80},  proc=makebordermovie
setvariable thresholdplas title="Threshold", size={150,20},pos={780,30}, value=peethreshold, limits={-10, 0, 0.1}
button showprofile title="Profile", size={60,15},pos={850,60}, proc=profileline
setvariable sb title = "Freq",size={140,20}, pos={810,110}, value = spectb,limits={0.01,5,0.01},proc = updatepizza
setvariable sw title = "Width",size={140,20}, pos={980,110}, value = spectw,limits={0,2.5,0.01},proc = updatepizza
button allprofiles title="Create all profiles", size={120,15},pos={850,80}, proc=allspectrumbands

setvariable fromv title = "from frame",size={140,20}, pos={560,60}, value = startsa,limits={0,numpnts(volume)-1,1}
setvariable tov title = "till",size={100,20}, pos={710,60}, value = endsa,limits={1,numpnts(volume),1}


ModifyGraph axThick(righttop)=0,freePos(righttop)=0
ModifyGraph noLabel(righttop)=2,noLabel(topright)=2,axThick(topright)=0

AppendToGraph/R=righttop/T=topright ywave vs xwave
xwave[1] = xwave[0]+length*cos(bordershown*2*pi/numpizzaslices)
ywave[1] = ywave[0]+length*sin(bordershown*2*pi/numpizzaslices)




MoveWindow/W=bordermovement 2, 2, 2, 2
ModifyGraph nticks(left)=4,minor(left)=0,sep(left)=1,fSize(left)=12
ModifyGraph fSize(topleft)=12,btLen(topleft)=3,freePos(topleft)=-10
ModifyGraph axisEnab(top)={0.02,0.4}
ModifyGraph axisEnab(topright)={0.7,0.95}
ModifyGraph axisEnab(righttop)={0.4,0.95}

Appendtograph/R=rightbottom/B=bottomright volume
ModifyGraph axisEnab(bottomright)={0.6,0.95}
ModifyGraph axisEnab(rightbottom)={0.05,0.35}
Appendtograph/R=rightbottom/B=bottomright indexlineyv vs indexlinex
AppendToGraph/R=rightbottom/B=bottomright newvolume
ModifyGraph rgb(newvolume)=(16386,65535,16385)



ModifyGraph mode(indexlineyv)=1
SetAxis rightbottom 0,*
ModifyGraph lsize(indexlineyv)=2,rgb(indexlineyv)=(0,0,0),lOptions(indexlineyv)=1
ModifyGraph freePos(bottomright)={0,left}

AppendToGraph/L=topleft/T indexlineym vs indexlinex
ModifyGraph mode(indexlineym)=3
ModifyGraph usePlusRGB(indexlineym)=1,plusRGB(indexlineym)=(0,0,0)

AppendToGraph/T indexlineys vs indexlinex
ModifyGraph mode(indexlineys)=1
ModifyGraph fSize(bottomright)=12
ModifyGraph highTrip(bottomright)=100000,lowTrip(bottomright)=1,btLen(bottomright)=3
modifyGraph mode(indexlineys)=1

ModifyGraph marker(indexlineym)=23,rgb(indexlineym)=(16386,65535,16385)
ModifyGraph mode(indexlineym)=8
indexlineym=1
indexlineym=3
ModifyGraph plusRGB(indexlineym)=(16386,65535,16385)
ModifyGraph mode(indexlineys)=3,marker(indexlineys)=17,msize(indexlineys)=5,rgb(indexlineys)=(16386,65535,16385)
ModifyGraph mode(indexlineys)=8,usePlusRGB(indexlineys)=1,plusRGB(indexlineys)=(16386,65535,16385)
ModifyGraph mode(indexlineyv)=8,marker(indexlineyv)=23,usePlusRGB(indexlineyv)=1,rgb(indexlineyv)=(16386,65535,16385),plusRGB(indexlineyv)=(16386,65535,16385)

//Legend/C/N=text1/J/F=0/M=1/A=RC/Y=-10/X=+10 "\\Z12\\s(volume) volume (based on border)\r\\s(newvolume) volume (based on intensity)\r\\Z12"
ModifyGraph margin(left)=40, margin(right)=40
Label left "Frequency"
ModifyGraph lblPos(left)=30,tlOffset(left)=-2
Label rightbottom "\\Z12 Volume (microl)";DelayUpdate
ModifyGraph fSize(rightbottom)=12,lblPos(rightbottom)=40,btLen(rightbottom)=2
ModifyGraph axisEnab(righttop)={0.45,0.95},axRGB(righttop)=(65535,65535,65535)

ModifyGraph margin(top)=39
Label top "Frame";DelayUpdate
ModifyGraph lblPos(top)=100
AppendToGraph/L=topy/T=topc spectrumband
ModifyGraph axisEnab(topy)={0.8,1}
ModifyGraph axisEnab(topc)={0.45,0.65}
SetAxis topy 0,*
ModifyGraph noLabel(topc)=2,axThick(topc)=0,freePos(topy)={0.45,kwFraction},freePos(topc)=0
end

Function changespectrumscale(S_Struct) : SliderControl
	STRUCT WMSliderAction &S_Struct
	NVAR maxspectrum
	
	ModifyImage spectrum ctab= {*,maxspectrum,ColdWarm,0}
	return 0
End


macro SaveVCAnalysis()


wavestats/Q/Z pressurepeak
if (exists("pressurepeak") <1)
doalert/T="VCA alert" 0,"Perform analysis first" 
else
Save/J/W/U={0,0,1,0}/I pressurepeak,pressureamplitude, maxconductanceevent, maxconductancetimepoint, nonvoiding, voidNo, voidnumber,voidstart,voidmaxtime,voidend,startpressure,maxpressure,endpressure,  maxdebiet,maxconductance,voidvolume,voidefficiency,voidduration, infusionspeed, debiet, conductance as filename

endif

end


threadsafe Function el3(w,x) : FitFunc
	Wave w
	Variable x

	//CurveFitDialog/ These comments were created by the Curve Fitting dialog. Altering them will
	//CurveFitDialog/ make the function less convenient to work with in the Curve Fitting dialog.
	//CurveFitDialog/ Equation:
	//CurveFitDialog/ variable result
	//CurveFitDialog/ if(x<w1)
	//CurveFitDialog/ result = w0*sqrt(w1^2-x^2) + w2
	//CurveFitDialog/ else
	//CurveFitDialog/ result = w2
	//CurveFitDialog/ endif
	//CurveFitDialog/ 
	//CurveFitDialog/ f(x) = result
	//CurveFitDialog/ End of Equation
	//CurveFitDialog/ Independent Variables 1
	//CurveFitDialog/ x
	//CurveFitDialog/ Coefficients 3
	//CurveFitDialog/ w[0] = w0
	//CurveFitDialog/ w[1] = w1
	//CurveFitDialog/ w[2] = w2

	variable result
	if(x<w[1])
	result = w[0]*sqrt(w[1]^2-x^2) + w[2]
	else
	result = w[2]
	endif
	
	return result
End


function calculatediameters(w)
wave w

make/O/N=(dimsize(w,0), (dimsize(w,1)/2)) diameters

diameters = w[p][q] + w[p][q+(dimsize(w,1)/2)]
end


function slopes (w, wt,n,ni)
wave w,wt
variable n,ni

wave W_coef = root:W_coef
CurveFit/Q/M=2/W=0 line, w[ni-n,ni+n]/X=wt[ni-n,ni+n]/D
return W_coef[1]
end

function fitslope(w,wt,n)
wave w,wt
variable n

duplicate/O w, fitdif
fitdif = -slopes(w,wt,n,p)

end

function calculatesurfacefromspline(border)
wave border

make/O/N=(dimsize(border,0)) surfacefromspline

surfacefromspline = calculatesurface(border,p)
end



function calculatesurface(border,frame)
wave border
variable frame

variable ww=dimsize(border,0), hh=dimsize(border,1)

wave xcenter, ycenter, xvalues, yvalues, M_ROIMask
makexy(border, frame)
imageboundaryToMask width = ww, height=hh, xwave=xvalues, ywave=yvalues, seedx=xcenter[0], seedy=ycenter[0]
return mean(M_ROIMask)
end

makexy(border, count)
imageboundaryToMask width = ww, height=hh, xwave=xvalues, ywave=yvalues, seedx=xcenter[0], seedy=ycenter[0]

function calculateinout(movie, border)
wave movie,border

wave xcenter, ycenter, xvalues, yvalues, M_ROIMask
make/O/N=(dimsize(movie,2)) intime, outtime

make/O/N=(dimsize(movie,0), dimsize(movie,1)) framein, frameout, volintensity
Make/O/D/C/N=(dimsize(movie,2)) gravity
variable count=0, ww=dimsize(movie,0), hh=dimsize(movie,1), tt=dimsize(movie,2)

do
makexy(border, count)
imageboundaryToMask width = ww, height=hh, xwave=xvalues, ywave=yvalues, seedx=xcenter[0], seedy=ycenter[0]
multithread framein = movie[p][q][count]*M_ROIMask[p][q]
intime[count]=mean(framein)

gravity[count]=centreofmass2D(framein)


multithread frameout = movie[p][q][count]*(1-M_ROIMask[p][q])

outtime[count]=mean(frameout)
//print count
count +=1

while (count<tt)

volintensity = framein+frameout
end

function calculatevolume(movie,framerate, fillrate)
wave movie
variable framerate, fillrate


variable pixelpf = dimsize(movie,0)*dimsize(movie,1)

make/O/N=(dimsize(movie,2)) intensity 
matrixOP/O columns = sumcols(movie)
matrixOP/O rows = sumrows(columns)
intensity = rows[0][0][p]
intensity /=pixelpf



duplicate/O intensity, newvolume
wave intensity_L, intensity_L_DIF, intensity_L_DIF_Hist, w_coef
Interpolate2/T=1/N=200/Y=intensity_L intensity
Differentiate intensity_L/D=intensity_L_DIF
Histogram/B={0,wavemax(intensity_L_DIF)/100,100}/DEST=intensity_L_DIF_Hist intensity_L_DIF,intensity_L_DIF_Hist
wavestats/Q intensity_L_DIF_Hist
variable/G intensityslope = V_maxloc
variable bgintensity = wavemin(intensity)
newvolume=(intensity-bgintensity)*fillrate/(intensityslope*framerate*60)

return bgintensity

end

function makexy(borderdata, frame)
wave borderdata
variable frame

wave xcenter, ycenter
variable numslices = dimsize(borderdata,1)

Differentiate/DIM=0 borderdata/D=differ;DelayUpdate
differ = abs(differ)

Smooth/B 10, differ;DelayUpdate
make/O/N=(numslices) xvalues, yvalues, zvalues

xvalues= borderdata[frame][p]*cos(p*2*pi/numslices)+xcenter[0]
yvalues= borderdata[frame][p]*sin(p*2*pi/numslices)+ycenter[0]
zvalues= differ[frame][p]
InsertPoints numslices,1, yvalues,xvalues
yvalues[numslices]=yvalues[0]

xvalues[numslices]=xvalues[0]
end



function makesplineborder(borderdata, points)
wave borderdata
variable points

variable numframes = dimsize(borderdata,0)
variable numslices = dimsize(borderdata,1)
duplicate/O borderdata bmplus
Insertpoints/M=1 numslices, 3, bmplus
bmplus[][numslices]=bmplus[p][0]
bmplus[][numslices+1]=bmplus[p][1]
bmplus[][numslices+2]=bmplus[p][2]

make/O/N=(numslices+3) border
make/O/N=(numframes, points) newborder
make/O/N=(points) xdataborder
xdataborder = p/points*numslices

variable count=0

do
border=bmplus[count][p]

Interpolate2/T=2/E=1/J=1/I=3/N=(points)/Y=border_SS/X=xdataborder border
newborder[count][] = border_SS[q]
count +=1
while (count<numframes)


end







function findborder(profiles,number)
wave profiles
variable number

make/O/N=(dimsize(profiles,0)) profile
profile = profiles[p][number]

Smooth/S=2 25, profile


Differentiate profile/D=profile_DIF
wavestats/Q profile_DIF



return V_minloc
end



function findend(sdprofiles,number)
wave sdprofiles
variable number

make/O/N=(dimsize(sdprofiles,0)) sdprofile
sdprofile = sdprofiles[p][number]
wavestats/Q sdprofile

//findlevel/Q/R=(V_maxloc) sdprofile, mean(sdprofile,0,0.75*V_maxloc)
//if (V_flag==0)
//return V_levelX
//else
//return 0
//endif
return V_maxloc
end

threadsafe function findendfit(profiles, number)
wave profiles
variable number

make/O/N=(2,3) Cmatrix=0
Cmatrix[0][2]=+1
Cmatrix[1][2]=-1
Make/O/N=2 Dwave


matrixOP/O profile = col(profiles,number)
Dwave[0]=wavemin(profile)+(wavemax(profile)-wavemin(profile))/10
DWave[1]=wavemin(profile)-(wavemax(profile)-wavemin(profile))/10
Smooth/S=2 5, profile
Make/D/N=3/O W_coef
W_coef = {(wavemax(profile)-wavemin(profile))/40,40,wavemin(profile)}
//print W_coef
//Make/O/T/N=2 T_Constraints
//T_Constraints[0] = {"K2 > "+num2str(wavemin(profile)-500),"K2 < "+num2str(wavemin(profile)+500)}
FuncFit/Q/N/W=2 el3 W_coef profile  /D/C={Cmatrix,Dwave}
//duplicate/O W_coef fitresults
//wave fit_profile = root:fit_profile
//duplicate/O fit_profile, fullfit
//variable crossing =abs(W_coef[1]-W_coef[3])

//CurveFit/Q/M=2/W=2 line, profile[round(crossing-8),round(crossing) ]/D

return W_coef[1]
end

threadsafe function findendfitb(profiles, number)
wave profiles
variable number

make/O/N=(2,3) Cmatrix=0
Cmatrix[0][2]=+1
Cmatrix[1][2]=-1
Make/O/N=2 Dwave

matrixOP/O profile = col(profiles,number)
Dwave[0]=wavemax(profile)/10
DWave[1]=wavemax(profile)/10


Smooth/S=2 5, profile
Make/D/N=3/O W_coef
W_coef = {wavemax(profile)/40,40,0}
//print W_coef
Make/O/T/N=2 T_Constraints
T_Constraints[0] = {"K2 > -5","K2 < 5"}
FuncFit/Q/N/W=2 el3 W_coef profile  /D/C={Cmatrix,Dwave}
//duplicate/O W_coef fitresults
//wave fit_profile = root:fit_profile
//duplicate/O fit_profile, fullfit
//variable crossing =abs(W_coef[1]-W_coef[3])

//CurveFit/Q/M=2/W=2 line, profile[round(crossing-8),round(crossing) ]/D

return W_coef[1]
end




function findslope(profiles,number)
wave profiles
variable number

make/O/N=(dimsize(profiles,0)) profile
profile = profiles[p][number]

Smooth/S=2 25, profile

Differentiate profile/D=profile_DIF
wavestats/Q profile_DIF


return V_min
end




function necky(movie,x1,x2,y1,y2,how)
wave movie
variable x1,x2,y1,y2, how


variable timer=STARTMStimer
Make/O/N=2 xwave,ywave
xwave = {x1,x2}
ywave = {y1,y2}

make/O/N=(dimsize(movie,2)) border, slope, ends, fitslope1,fitslope2

ImageLineprofile/S/P=-2 srcwave = movie, xwave = xwave, ywave=ywave, width = 10
wave M_imagelineprofile = root:M_imagelineprofile
wave M_LineProfileStdv = root:M_LineProfileStdv
Differentiate/DIM=0 M_ImageLineProfile/D=DIFprofile


wave fitresults = root:fitresults
variable count=0
make/O/N=(dimsize(M_ImageLineProfile,0)) DIFline
do

if(how==0)
DIFline=DIFprofile[p][count]
Findpeak/Q/R=[numpnts(DIFline)-1,0]/N/M=(0.6*(wavemin(DIFline))) DIFline
	if(V_flag==0)
	ends[count]=max(V_peakloc,0)
	else
	wavestats/Q DIFline
	ends[count]=V_minloc
	endif
endif
if(how==1)
ends[count]= max(min(findendfit(M_ImageLineProfile,count)-5,150),0)
endif

if(how==2)
ends[count]= max(min(findendfitb(M_ImageLineProfile,count)-5,150),0)
endif


count+=1
while (count<numpnts(fitslope1))


print STOPMSTimer(timer)/1000000
end

function necky2(movie,x1,x2,y1,y2,how,cores)
wave movie
variable x1,x2,y1,y2, how, cores

NVAR length
variable timer=STARTMStimer
Make/O/N=2 xwave,ywave
xwave = {x1,x2}
ywave = {y1,y2}
controlinfo multicore
variable mc=V_value
make/O/N=(dimsize(movie,2)) border, slope, ends, fitslope1,fitslope2

ImageLineprofile/S/P=-2 srcwave = movie, xwave = xwave, ywave=ywave, width = 10
wave M_imagelineprofile = root:M_imagelineprofile
wave M_LineProfileStdv = root:M_LineProfileStdv
Differentiate/DIM=0 M_ImageLineProfile/D=DIFprofile


//wave fitresults = root:fitresults
variable count=0
make/O/N=(dimsize(M_ImageLineProfile,0)) DIFline

switch(how)
case 0:
do

DIFline=DIFprofile[p][count]
Findpeak/Q/R=[numpnts(DIFline)-1,0]/N/M=(0.6*(wavemin(DIFline))) DIFline
	if(V_flag==0)
	ends[count]=max(V_peakloc,0)
	else
	wavestats/Q DIFline
	ends[count]=V_minloc
	endif
count+=1
while (count<numpnts(fitslope1))
break

case 1:

multithread/NT=(cores) ends= max(min(findendfit(M_ImageLineProfile,p)-5,length),0)

break

case 2:

multithread/NT=(cores) ends= max(min(findendfitb(M_ImageLineProfile,p)-5,length),0)

break

endswitch
NVAR framerate


print STOPMSTimer(timer)/1000000 
end






function calculatelength (DIFprofile, borderpoints)
wave DIFprofile, borderpoints


duplicate/O DIFprofile, sqrtm

sqrtm=sqrt(DIFprofile[p][q]^2+1)
sqrtm = p<borderpoints[q] ? sqrtm[p][q] : 0


matrixOP/O lengthse = sumcols(sqrtm)
make/O/N=(dimsize(lengthse,1)) lengthe = lengthse[0][p]

lengthe=lengthse[0][p]
end



function deletemovie(points,bordertc,angle, x0,y0)
wave points, bordertc
variable angle, x0,y0


variable count=0, xr, yr

do
xr = x0+round(bordertc[count]*cos(angle))
yr = y0+round(bordertc[count]*sin(angle))
points[xr-1,xr+1][yr-1,yr+1][count] = NaN

count+=1

while(count<dimsize(bordertc,0))

end



function makemovie(points,bordertc,slopetc,angle, x0,y0)
wave points, bordertc, slopetc
variable angle, x0,y0


variable count=0, xr, yr

do
xr = x0+round(bordertc[count]*cos(angle))
yr = y0+round(bordertc[count]*sin(angle))
points[xr-1,xr+1][yr-1,yr+1][count] = abs(slopetc[count])<0.5 ? 0 : slopetc[count]

count+=1

while(count<dimsize(bordertc,0))

end


macro replacemborderpoint(startp,endp)
variable startp, endp



string bordername, borderdifname, FFplotname
replacebymneighbours(bordermatrix, startp, endp)
makesplineborder(bordermatrix,200)
surfformula[]=surfaces2(bordermatrix,p)

variable count=0

do
bordername = "borderpoint"+num2str(count)
borderdifname = "borderpointdif"+num2str(count)
FFplotname = "FFplot"+num2str(count)
deletemovie(points, $bordername,count*2*pi/numpizzaslices, xcenter[0],ycenter[0])
$bordername = bordermatrix[p][count]
Differentiate $bordername/D=$borderdifname
//Smooth 5, $borderdifname
STFT/SEGS=64/OUT=6/HOPS=1/WINF=Cos1/DEST=$FFplotname $borderdifname

SetScale/P y 0,1, "Frame", $FFplotname
Spectrum3D[][][count]=$FFplotname[p][q]
setscale/P y 0,dimdelta($FFplotname,1), "Hz", spectrum3D
spectrum = $FFPlotname
bordermove = $borderdifname
Smooth 5, bordermove
//bordermove *=(framerate*pixelsize/1000)

makemovie(points, $bordername,$borderdifname,count*2*pi/numpizzaslices, xcenter[0],ycenter[0])
count+=1

while (count<dimsize(bordermatrix,1))

makexy(newborder,slidervalue)

Matrixop/O framep = points[][][slidervalue]
volume = (4/3)*pi*(surfformula/pi)^1.5
volume *=(pixelsize/1000)^3
duplicate/O volume, bettervolume

bettervolume = (10^(-9))*(pixelsize^3)*(4/3)*surfformula*(borderpoint3+borderpoint15)/2
volume = bettervolume
end


macro replaceborderpoints(test): buttoncontrol
string test

replacemborderpoint(firsttodelete, lasttodelete)
end



macro replaceborderpoint(test): buttoncontrol
string test



string bordername, borderdifname, FFplotname
replacebyneighbours(bordermatrix, bordershown)
makesplineborder(bordermatrix,200)
surfformula[]=surfaces2(bordermatrix,p)

bordername = "borderpoint"+num2str(bordershown)
borderdifname = "borderpointdif"+num2str(bordershown)
FFplotname = "FFplot"+num2str(bordershown)
deletemovie(points, $bordername,bordershown*2*pi/numpizzaslices, xcenter[0],ycenter[0])
$bordername = bordermatrix[p][bordershown]
Differentiate $bordername/D=$borderdifname
//Smooth 5, $borderdifname
STFT/SEGS=64/OUT=6/HOPS=1/WINF=Cos1/DEST=$FFplotname $borderdifname
SetScale/P x 0,1,"Frame", $FFplotname
Spectrum3D[][][count]=$FFplotname[p][q]
setscale/P y 0,dimdelta($FFplotname,1), "Hz", spectrum3D

matrixop/O spectrum = layer(spectrum3d,count)
setscale/P y 0,dimdelta(spectrum3D,1), "Hz", spectrum
bordermove = $borderdifname
Smooth 5, bordermove
//bordermove *=framerate*pixelsize/1000
makexy(newborder,slidervalue)

makemovie(points, $bordername,$borderdifname,bordershown*2*pi/numpizzaslices, xcenter[0],ycenter[0])
Matrixop/O framep = points[][][slidervalue]
volume = (4/3)*pi*(surfformula/pi)^1.5
volume *=(pixelsize/1000)^3

duplicate/O volume, bettervolume

bettervolume = (10^(-9))*(pixelsize^3)*(4/3)*surfformula*(borderpoint2+borderpoint14)/2
volume = bettervolume

end

macro restoreborderpoint(test): buttoncontrol
string test

string bordername, borderdifname, FFplotname
bordermatrix[][bordershown]=bmoriginal[p][bordershown]
makesplineborder(bordermatrix,200)
surfformula[]=surfaces2(bordermatrix,p)

bordername = "borderpoint"+num2str(bordershown)
borderdifname = "borderpointdif"+num2str(bordershown)
FFplotname = "FFplot"+num2str(bordershown)
deletemovie(points, $bordername,bordershown*2*pi/numpizzaslices, xcenter[0],ycenter[0])
$bordername = bordermatrix[p][bordershown]
Differentiate $bordername/D=$borderdifname
//Smooth 5, $borderdifname
STFT/SEGS=64/OUT=6/HOPS=1/WINF=Cos1/DEST=$FFplotname $borderdifname
SetScale/P x 0,1,"Frame", $FFplotname
Spectrum3D[][][count]=$FFplotname[p][q]
setscale/P y 0,dimdelta($FFplotname,1), "Hz", spectrum3D
matrixop/O spectrum = layer(spectrum3d,count)
setscale/P y 0,dimdelta(spectrum3D,1), "Hz", spectrum
bordermove = $borderdifname
Smooth 5, bordermove
//bordermove *=framerate*pixelsize/1000
makexy(newborder,slidervalue)

makemovie(points, $bordername,$borderdifname,bordershown*2*pi/numpizzaslices, xcenter[0],ycenter[0])
Matrixop/O framep = points[][][slidervalue]
volume = (4/3)*pi*(surfformula/pi)^1.5
volume *=(pixelsize/1000)^3

duplicate/O volume, bettervolume

bettervolume = (10^(-9))*(pixelsize^3)*(4/3)*surfformula*(borderpoint3+borderpoint15)/2
volume = bettervolume
end



function replacebyneighbours(bw,sl)
wave bw
variable sl

variable nc=dimsize(bw,1)



if (sl==0)
bw[][sl]=(bw[p][nc-1] + bw[p][1])/2

elseif (sl==(nc-1))

bw[][sl]=(bw[p][nc-2] + bw[p][0])/2

else
bw[][sl]=(bw[p][sl-1] + bw[p][sl+1])/2
endif
end


function replacebymneighbours(bw,sl1,sl2)
wave bw
variable sl1,sl2

matrixop/O bw2 = bw

variable np = dimsize(bw,1)

insertpoints/M=1 np,(2*np), bw2

bw2[][np,(2*np-1)] = bw[p][q-np]
bw2[][(2*np),(3*np-1)] = bw[p][q-2*np]


variable nc=dimsize(bw,1), distance=abs(sl1-sl2), count=0

if (sl1==sl2)
replacebyneighbours(bw,sl1)
endif


if (sl1<sl2)

do

bw[][sl1+count]=(((distance+1-count)/(distance+2))*bw2[p][sl1+np-1] + ((Count+1)/(distance+2))*bw2[p][sl2+np+1])


count +=1

while (count<(distance+1))
endif


if (sl1>sl2)

distance = np-distance
do


bw[][mod((sl1+count),np)]=(((distance+1-count)/(distance+2))*bw2[p][sl1+np-1] + ((Count+1)/(distance+2))*bw2[p][sl2+2*np+1])


count +=1

while (count<(distance+1))
endif





end




macro makepolygons(npoints, nframes)
variable npoints, nframes

string title
make/O/N=(nframes, npoints) xvalues, yvalues

variable counter=0

do
title = "borderpoint"+num2str(counter)
xvalues[][counter] = $title[p]*cos(counter*2*pi/npoints)
yvalues[][counter] = $title[p]*sin(counter*2*pi/npoints)

counter +=1

while (counter < npoints)

end

macro surfacefromfit(npoints)
variable npoints

pauseupdate
variable counter=0
make/O/N=24 ellipseX, ellipseY, ellipseXFit,ellipseYFit
Make/D/O ellipseCoefs={10,10,0,0, 110*pi/180}	
make/O/N=(npoints) fittedsurface
do
ellipseCoefs={10,10,0,0, 110*pi/180}	
ellipseX = xvalues[counter][p]
ellipsey = yvalues[counter][p]

FuncFit/Q/ODR=3 FitEllipse2, ellipseCoefs /X={ellipseX, ellipseY}/XD={ellipseXFit,ellipseYFit}
fittedsurface[counter]=pi*ellipseCoefs[0]*ellipsecoefs[1]
//print w_coef[0], w_coef[1]
counter +=1
while (counter<npoints)
end

function curvature(l0,l1,l2,angle)
variable l0,l1,l2,angle

variable x0,x1,x2,y0,y1,y2

x0 = l1
y0=0
x1=l1*cos(angle)
y1=l1*sin(angle)
x2=l2*cos(2*angle)
y2=l2*sin(2*angle)

variable s0,s1, d
s0=(y1-y0)/(x1-x0)
s1=(y2-y1)/(x2-x1)

d=sqrt( ((y2-y0)/2)^2 + ((x2-x0)/2)^2)

return ((s1-s0)/d)

end

macro VCAtofront(test):buttoncontrol
string test

dowindow/F  VCA

movewindow/W=VCA 1,1,1,1
movewindow/W=bordermovement 2,2,2,2
dowindow/B=VCA bordermovement
end

macro findcenter(test): buttoncontrol
string test



matrixOP/O average=sumbeams(movie2)
maxframe(movie2,average)
Dowindow/K centerwindow
Display/W=(50,300,450,600)/N=centerwindow
AppendImage average


xcenter = dimsize(average,0)/2
ycenter = dimsize(average,1)/2

//print xcenter, ycenter


ModifyGraph mode=3,marker=1,msize=12
appendtograph/C=(65535,0,0) ycenter vs xcenter

ModifyGraph mode=3,marker=1,msize=14

length = min(xcenter[0],ycenter[0])
length -=10

drawcircle("",0,"","")

appendtoGraph uppercircle vs xcircle
appendtoGraph lowercircle vs xcircle




Dowindow/K centerwindow2
Display/W=(450,300,850,600)/N=centerwindow2
AppendImage maximum
appendtoGraph uppercircle vs xcircle
appendtoGraph lowercircle vs xcircle
//WMAppend3DImageSlider()



end

function maxframe(w,vb)
wave w,vb

matrixOP/O maximum=vb 
multithread maximum = maxbeam(w,p,q)

end

threadsafe function maxbeam(w,x1,y1)
wave w
variable x1,y1

matrixOP/O be=beam(w,x1,y1)
return wavemax(be)

end

Function Show() //displays the movie with slider and slice
    Wave movie2, points
    NVAR bordershown, numpizzaslices
   
    Matrixop/O framem = movie2[][][slidervalue]

    appendimage/W=bordermovement/L/T framem
    variable/G slidervalue, bordershown2
        Matrixop/O framep = points[][][slidervalue]
    appendimage/W=bordermovement/T/L framep
    ModifyImage framep ctab={-1,1,Rainbow256,0}
    ModifyGraph axisEnab(top)={0,0.5}
       ModifyGraph axisEnab(left)={0.5,1}
    movewindow/W=win 50, 50, 1400, 800
    controlbar 100
  
    slider sld, limits={0,(dimsize(movie2,2)-1),1}, ticks =0,  size={300,20},vert=0,variable = slidervalue,proc=ActionProcName
    SetVariable setvar0,pos={520,18},size={122,17},limits={0,(dimsize(movie2,2)-1),1},value= slidervalue
     
End


Function updateframe()
    Wave movie2, points, indexlinex, newborder
    NVAR slidervalue
    Controlinfo sld
    Variable value = v_value
    Matrixop/O framem = movie2[][][value]
    Matrixop/O framep = points[][][value]
    slidervalue = value
    makexy(newborder,value)
   
    
    doupdate
    indexlinex=slidervalue
    
end

Function updateframebis(x0)
	variable x0
    Wave movie2, points, indexlinex
   NVAR slidervalue
   if((x0>-1)&(x0<(dimsize(points,2)-1)))
   slidervalue = x0
    
    Matrixop/O framem = movie2[][][x0]
    Matrixop/O framep = points[][][x0]
    doupdate
   else
   x0= dimsize(points,2)-1
 	slidervalue = x0
   
    Matrixop/O framem = movie2[][][x0]
    Matrixop/O framep = points[][][x0]
    doupdate
   
   
    endif
     indexlinex=slidervalue
end




Function updatepizza(SV_Struct):setvariablecontrol
		STRUCT WMSetVariableAction &SV_Struct
	
	
	
	
	Wave xwave,ywave, spectrum, bordermove,  spectrum3D, framerate, pixelsize
	
	
	
	Nvar numpizzaslices, length, bordershown
	
 
	
	
	controlinfo setvar3
	Variable value = V_value
	if(value==-1)
	value = numpizzaslices-1
	bordershown = value
	endif
	
	if(value==numpizzaslices)
	value=0
	bordershown = value
	endif
	

	xwave[1] = xwave[0]+length*cos(value*2*pi/numpizzaslices)
	ywave[1] = ywave[0]+length*sin(value*2*pi/numpizzaslices)

		
	
string borderwave, borderFF, bordermeanFF
borderwave = "borderpointdif"+num2str(value)
borderFF = "FFplot"+num2str(value)
bordermeanFF = "FFmean"+num2str(value)
WAVE w= $borderFF, pw=$borderwave

matrixop/O spectrum = layer(spectrum3d,bordershown)
setscale/P y 0,dimdelta(spectrum3D,1), "Hz", spectrum
bordermove = pw
Smooth 5, bordermove
//bordermove *=framerate*pixelsize/1000
if (exists("normspec")==1)
WAVE qw=$bordermeanFF, normspec, profileyt, profileyb
NVAR spectb, spectw
normspec =qw
makeprofile(normspec, spectb,spectw)
wave W_imagelineprofile, spectrumband
spectrumband = W_imagelineprofile 

	profileyt=spectb+spectw/2
	profileyb=spectb-spectw/2


endif

		
end	

macro allspectrumbands(test): buttoncontol
string test

string bordermeanFF, spectrumbandname
variable count=0

do
bordermeanFF = "FFmean"+num2str(count)
spectrumbandname="Spectrumband"+num2str(count)
makeprofile($bordermeanFF, spectb,spectw)
duplicate/O W_imagelineprofile, $spectrumbandname

count+=1
while(count<numpizzaslices)
end

Function updateframe2(SV_Struct):setvariablecontrol
		STRUCT WMSetVariableAction &SV_Struct
 NVAR slidervalue

updateframe()


getaxis/Q top	
Variable dx= (V_max-V_min)/2

SetAxis top slidervalue-dx, slidervalue+dx
		
end	




Function ActionProcName(sa) : SliderControl
    STRUCT WMSliderAction &sa

    sa.blockReentry=1
    if(sa.eventcode == 9)
        updateFrame()
    endif
End
Function ActionProcName2(sa) : SliderControl
    STRUCT WMSliderAction &sa

    sa.blockReentry=1
    if(sa.eventcode == 9)
     //  updatepizza()
    endif
End

function drawcircle(ctrlName,varNum,varStr,varName) : SetVariableControl
	String ctrlName
	Variable varNum	// value of variable as number
	String varStr		// value of variable as string
	String varName	// name of variable
NVAR length
WAVE xcenter = root:xcenter, ycenter = root:ycenter


make/O/N=(length*2+1) uppercircle, lowercircle, xcircle
xcircle = p+xcenter[0]-length
uppercircle = sqrt(length^2-(p-length)^2 )+ycenter[0]
lowercircle = -sqrt(length^2-(p-length)^2 )+ycenter[0]

end


function recalculatevolint(ctrlName,varNum,varStr,varName) : SetVariableControl
	String ctrlName
	Variable varNum	// value of variable as number
	String varStr		// value of variable as string
	String varName	// name of variable
NVAR fillrate, usedfr, bgintensity, usedbg, framerate, usedframerate, intensityslope
WAVE newvolume, usednv, intensity

//newvolume =(usednv*fillrate/usedfr - (bgintensity-usedbg)*fillrate/usedfr)*(usedframerate/framerate)
newvolume=(intensity-bgintensity)*fillrate/(intensityslope*framerate*60)
end

function recalculatevolsur(ctrlName,varNum,varStr,varName) : SetVariableControl
	String ctrlName
	Variable varNum	// value of variable as number
	String varStr		// value of variable as string
	String varName	// name of variable
NVAR pixelsize, usedps
WAVE volume, usedv

volume =usedv*((pixelsize/usedps)^3)

end

function recalculatefr(ctrlName,varNum,varStr,varName) : SetVariableControl
	String ctrlName
	Variable varNum	// value of variable as number
	String varStr		// value of variable as string
	String varName	// name of variable
NVAR framerate, usedframerate, bordershown, fillrate, usedfr, bgintensity, usedbg, intensityslope, numpizzaslices
wave spectrum3D, newvolume, usednv, intensity, normspec
variable spacing = dimdelta(spectrum3D, 1)
setscale/P y 0,(spacing*framerate/usedframerate), "Hz", spectrum3D


matrixOP/O spectrum=layer(spectrum3D, bordershown)
setscale/P y 0,dimdelta(spectrum3D,1), "Hz", spectrum

string meanspectrumname
variable count=0

do
meanspectrumname="FFmean"+num2str(count)

SetScale/P y 0,dimdelta(spectrum3D,1),"Hz", $meanspectrumname
SetScale/P y 0,dimdelta(spectrum3D,1),"Hz", normspec
 count+=1
 
while (count<numpizzaslices) 



//newvolume =(usednv*fillrate/usedfr - (bgintensity-usedbg)*fillrate/usedfr)*(usedframerate/framerate)
newvolume=(intensity-bgintensity)*fillrate/(intensityslope*framerate*60)
usedframerate = framerate
end


macro pizza (numslices,  pixelsize, method)
variable numslices=8,  pixelsize = 51.423, method=0


 






variable x0,y0

y0 = ycenter[0]
x0 = xcenter[0]







string bordername, borderdifname, slopename, slopedifname, FFplotname, arcname

variable counter = 0
variable x1,y1


//duplicate/O movie2, points
//Redimension/S points

matrixOP/O points = movie2*NaN

pauseupdate
do
print counter
bordername = "borderpoint"+num2str(counter)
borderdifname = "borderpointdif"+num2str(counter)
arcname = "arc"+num2str(counter)
slopename = "slope"+num2str(counter)
slopedifname="slopedif"+num2str(counter)
x1 = x0+length*cos(counter*2*pi/numslices)
y1 = y0+length*sin(counter*2*pi/numslices)
necky2(movie2,x0,x1,y0,y1, method, numcores)
//print x0,x1,y0,y1
duplicate/O ends, $bordername
//calculatelength(DIFprofile, $bordername)
//duplicate/O lengthe, $arcname
duplicate/O fitslope1, $slopename
Differentiate $slopename/D=$slopedifname
Differentiate $bordername/D=$borderdifname
//Smooth 5, $borderdifname
makemovie(points, $bordername,$borderdifname,counter*2*pi/numslices, x0,y0)



counter +=1
while (counter<numslices)


duplicate/O ends, surface

surface =0
counter=0

make/O/N=(numpnts(ends), numslices) bordermatrix
make/O/N=(dimsize(movie2, 2), 65, numslices) spectrum3D
do
bordername = "borderpoint"+num2str(counter)
borderdifname = "borderpointdif"+num2str(counter)
FFplotname = "FFplot"+num2str(counter)
surface +=$bordername^2*pi/(numslices)

bordermatrix[][counter]=$bordername[p]

duplicate/O $borderdifname movementtrace
setscale/P x 0,1/framerate, movementtrace

STFT/SEGS=64/OUT=6/HOPS=1/WINF=Cos1/DEST=$FFplotname movementtrace
SetScale/P x 0,1,"Frame", $FFplotname
spectrum3D[][][counter] = $FFplotname[p][q]
setscale/P y 0,dimdelta($FFplotname,1), "Hz", spectrum3D
SetScale/P x 0,1,"Frame", $FFplotname
counter+=1
while (counter<numslices)
duplicate/O surface volume


make/O/N=(numpnts(ends)) surfformula
make/O/N=(numslices) coordinates

counter=0


surfformula[]=surfaces2(bordermatrix,p)




volume = (4/3)*pi*(surfformula/pi)^1.5
volume *=(pixelsize/1000)^3
duplicate/O volume, bettervolume

bettervolume = (10^(-9))*(pixelsize^3)*(4/3)*surfformula*(borderpoint3+borderpoint8)/2
volume = bettervolume


duplicate/O bordermatrix, bmoriginal

//XLLoadWave/S="Sheet1"/R=(A1,B10002)/O/D/T/N=CTtime 

//display volume vs CTtime0
end



macro calculatecurves(numslices)
variable numslices

variable volgend, vorig
variable counter =0

string vorigborder, huidigborder, volgendborder, curvaturename, contractname
numslices -=1
do
curvaturename = "curvature"+num2str(counter)
contractname = "contract"+num2str(counter)
volgend = counter +1
if (volgend>numslices)
volgend = 0
endif 
vorig =counter-1
if (vorig<0)
vorig = numslices
endif

vorigborder = "borderpoint"+num2str(vorig)
huidigborder = "borderpoint"+num2str(counter)
volgendborder = "borderpoint"+num2str(volgend)

duplicate/O $huidigborder, $curvaturename
$curvaturename = ($huidigborder-($vorigborder[p]+$volgendborder[p])/2)/$huidigborder
duplicate/O $curvaturename, $contractname
Differentiate $curvaturename/D=$contractname
Smooth/B 5, contractname
counter +=1
while (counter<numslices)


end


function surfaces2(dw,row)
wave dw
variable row

make/O/N=(dimsize(dw,1)) d

d[]=dw[row][p]

variable angle = 2*pi/numpnts(d)

variable counter=0

duplicate/O d dx,dy, shiftup, shiftdown
dx=d[p]*cos(angle*p)
dy=d[p]*sin(angle*p)
InsertPoints (numpnts(dy)),1, dy
dy[numpnts(dy)-1]=dy[0]

shiftup=dx[p]*dy[p+1]


deletePoints (numpnts(dy)-1),1, dy
InsertPoints 0,1, dx
dx[0]=dx[numpnts(dx)-1]
shiftdown=dx[p+1]*dy[p]
deletePoints 0,1, dx
variable areas

areas = sum(shiftup)-sum(shiftdown)

return areas
end




function surfaces(d)
wave d

variable angle = 2*pi/numpnts(d)

variable counter=0

duplicate/O d dx,dy, shiftup, shiftdown
dx=d[p]*cos(angle*p)
dy=d[p]*sin(angle*p)
InsertPoints (numpnts(dy)),1, dy
dy[numpnts(dy)-1]=dy[0]

shiftup=dx[p]*dy[p+1]


deletePoints (numpnts(dy)-1),1, dy
InsertPoints 0,1, dx
dx[0]=dx[numpnts(dx)-1]
shiftdown=dx[p+1]*dy[p]
deletePoints 0,1, dx
variable areas

areas = sum(shiftup)-sum(shiftdown)

return areas
end




macro VCAnalyser()

doWindow/K firstgraph
dowindow/K resultstable
dowindow/K scrollfigure
dowindow/K overviewfigure
XLLoadWave/S="Sheet1"/R=(A1,E200000)/O/D/N=wave
deletepoints 0,3, wave0, wave1, wave2, wave3, wave4

string/G filename = "VCA_"+S_filename
Display/W=(10,0,460,390)/N=firstgraph wave1 vs wave0
AppendToGraph wave2 vs wave0
ModifyGraph rgb(wave2)=(16385,16388,65535)
TextBox/C/N=text0/A=LB/X=2.73/Y=57.24 "\\K(0,65535,0)Put cursors for reference period!"

SetAxis left -5,*
showinfo
cursor A wave1 10
cursor B wave1 1000
variable/G minpamp = 10, debietthresholdmultiplier = 3, filterfactor=1

button startanalysis, proc=voidanalyzer, title="Start analysis", size = {150,30}, pos={100,10}, fsize=16, fcolor=(65535,0,0)
setvariable minamp size={180,20}, limits={0,100,5}, pos={300,10}, title = "Min pressure amplitude", value = minpamp
setvariable filterfact size={180,20}, limits={0.1,5,0.1}, pos={300,30}, title = "Filter factor", value = filterfactor
setvariable dtmp size={180,20}, limits={1,100,.5}, pos={300,50}, title = "Flow rate sensitivity", value = debietthresholdmultiplier
end




macro voidanalyzer(test) : buttoncontrol
string test
dowindow/K resultstable
dowindow/K scrollfigure
dowindow/K overviewfigure


silent 1; pauseupdate
variable/G refstart, refend, samplingrate=wave0[2]-wave0[1]
refstart=min(xcsr(A),xcsr(B))
refend=max(xcsr(A),xcsr(B))


make/O/N=200 gausswidth=0, pressurestart=0, maxconductanceevent=0, maxconductancetimepoint=0, pressurestop=0,pressurepeak=0, pressurepeaklevel=0, pressureamplitude=0, voidNo=NaN, voidmaxx=0, maxdebiet=0, voidmaxtime=0, voidstart=0, voidend=0, residualvolume=0,  maxconductance=0, startpressure=0, maxpressure=0, endpressure=0, voidvolume=0, voidefficiency=0, capacity=0

variable count=0

variable lastpeak, timespace, slopenumber

//Differentiate wave1/X=wave0/D=debiet

duplicate/O wave1, debiet
timespace = wave0[201]-wave0[200]

slopenumber = round(0.5/timespace)

debiet = slopes(wave1,wave0,slopenumber*filterfactor,p)

debiet = -debiet
make/O debiet_Hist
Histogram/B={-20,0.2,200} debiet,debiet_Hist;DelayUpdate
CurveFit/Q/M=2/W=0 gauss, debiet_Hist[0,200]

debiet -=W_coef[2]
Make/O/N=1 infusionspeed
infusionspeed[0] = -W_coef[2]


duplicate/O debiet, conductance

variable threshold = debietthresholdmultiplier*W_coef[3]
print threshold
conductance = debiet/wave2
duplicate/O conductance, conductancec
conductancec = debiet>threshold ? conductance : 0
conductance = conductancec



variable startx, endx


count=0
lastpeak = refstart
do

findpeak/Q/M=(threshold)/R=[min(lastpeak+round(2/timespace), numpnts(debiet)-1),numpnts(debiet)] debiet
if(V_flag==0)
voidmaxx[count] = V_Peakloc
voidmaxtime [count] = wave0(V_Peakloc)

findlevel/Q/R=(voidmaxx[count],0) debiet,(V_peakval/5)
voidstart[count] = wave0(V_LevelX)
startx=V_levelx
startpressure[count]=wave2[V_levelX]
findlevel/Q/R=(voidmaxx[count],) debiet,0
voidend[count] = wave0(V_LevelX) 
endx=V_levelx




endpressure[count]=wave2[V_levelX]
maxpressure[count]=wavemax(wave2,startx, endx)
maxdebiet[count] = wavemax(debiet, startx, endx)

CurveFit/Q/M=2/W=0 gauss, debiet[max(0,round(voidmaxx[count]-20/samplingrate)),round(voidmaxx[count]+20/samplingrate)]/X=wave0[max(0,round(voidmaxx[count]-20/samplingrate)),round(voidmaxx[count]+20/samplingrate)]/D
gausswidth[count]=W_coef[3]*2.355
maxconductance[count]=wavemax(conductance,startx,endx)
voidvolume[count] = wave1[startx]-wave1[endx]
voidefficiency[count] = voidvolume[count]/wave1[startx]
residualvolume[count]=wave1[endx]
capacity[count]=wave1[startx]
//voidedvolume 
else
break




endif


//lastpeak = V_peakloc + 50
lastpeak = endx

count+=1
while((count<200)&&(lastpeak+round(2/timespace)<numpnts(debiet)))

deletepoints count, 200-count, voidmaxx,voidmaxtime, voidstart, voidend,gausswidth,  maxconductance, maxdebiet, maxpressure, startpressure, endpressure, voidvolume, voidefficiency, residualvolume, capacity
duplicate/O voidend, voidduration
voidduration -= voidstart

duplicate/O voidmaxtime, ici
ici[1,]=voidmaxtime[p]-voidmaxtime[p-1]
ici[0]=Nan
duplicate/O voidend, voidnumber
voidnumber[] = p

differentiate wave2/X=wave0/D=pressurechange
duplicate/O pressurechange, sq

pressurechange =slopes(wave2,wave0,slopenumber,p)
sq=pressurechange^2

smooth/B 10,sq
wavestats/Q/R=[refstart, refend] sq
variable/G pressurechangethreshold = V_sdev*3



lastpeak = refstart
//print firstpoint
variable pp, pressurem, pressuren, tempx1, tempx2
count=0
do
findpeak/Q/B=(1/timespace)/M=3/R=(lastpeak,) sq
//print timespace
pp = V_peakloc

if(V_flag==0)


findlevel/Q/R=(pp-1,0) sq, pressurechangethreshold
pressurestart[count]=wave0(V_levelX)
tempx1=V_levelX
pressurem=wave2[V_levelx]
findlevel/Q/R=(pp+1,) sq, pressurechangethreshold
tempx2 = V_levelX
pressuren=wave2[V_levelx]
pressurestop[count]=wave0(V_levelX)
wavestats/Q/R=(tempx1,tempx2) wave2
pressurepeaklevel[count] = V_max
pressureamplitude[count]=V_max-min(pressurem, pressuren)
pressurepeak[count]=wave0(V_maxloc)
wavestats/Q/R=(tempx1, tempx2) conductance
maxconductanceevent[count]=V_max
maxconductancetimepoint[count]=wave0(V_maxloc)
lastpeak=V_levelX+1
count +=1

else
break
endif

if(pressureamplitude[count-1]<minpamp)
pressurestart[count-1]=0
pressurestop[count-1]=0
pressureamplitude[count-1]=0
pressurepeak[count-1]=0
count-=1
endif



while(count<200)

deletepoints count, 200-count, pressurepeak,pressurepeaklevel, pressurestart, maxconductanceevent, maxconductancetimepoint, pressurestop, pressureamplitude, voidNo
duplicate/O pressurestart, onesp, nonvoiding

wavestats/Q pressurepeak
if (V_npnts <1)
doalert/T="VCA alert" 0,"No pressure peaks detected! Adjust the reference period!"
abort
endif
button startanalysis fcolor=(0,0,0)
nonvoiding=1
duplicate/O voidNo,eventnumber
eventnumber[]=p
onesp=100

//delete too small contractions



//find nonvoiding contractions

variable countp=0, countv=0

do
if((voidmaxtime[countv]>pressurestart[countp])&&(voidmaxtime[countv]<pressurestop[countp]))
nonvoiding[countp]=0
voidNo[countp]=countv
countp +=1
countv=0
else
	if(countv==(numpnts(voidmaxtime)-1))
	countp +=1
	countv=0
	else
	countv +=1
	endif
endif

while(countp <numpnts(pressurestart))





Edit/N=resultstable/W=(20,460,1400,650) voidnumber,maxconductance,maxdebiet, maxpressure, startpressure,voidduration,gausswidth, voidefficiency,voidmaxtime,voidstart,voidend,voidvolume, residualvolume, capacity, ici, endpressure, eventnumber,pressurepeak, nonvoiding, pressureamplitude, voidNo, infusionspeed

Duplicate/O voidstart ones
ones=100

Display/N=overviewfigure/W=(470,0,920,390) wave2 vs wave0
ModifyGraph rgb(wave2)=(0,0,0)
ModifyGraph axisEnab(left)={0,0.3333}
AppendToGraph pressurepeaklevel vs pressurepeak
ModifyGraph mode(pressurepeaklevel)=3,marker(pressurepeaklevel)=23,rgb(pressurepeaklevel)=(1,4,52428)
ModifyGraph zColor(pressurepeaklevel)={nonvoiding,0,1,RedWhiteBlue256,0}
AppendToGraph/R debiet vs wave0
setaxis right -wavemax(debiet)/5, wavemax(debiet)
ModifyGraph rgb(debiet)=(1,4,52428)
ModifyGraph axisEnab(right)={0.36,0.6667}
AppendToGraph/R maxdebiet vs voidmaxtime
ModifyGraph mode(maxdebiet)=3,marker(maxdebiet)=1,mrkThick(maxdebiet)=1,rgb(maxdebiet)=(0,65535,0)
//AppendToGraph/L=L2 maxconductance vs voidmaxtime
Appendtograph/L=L2 conductancec vs wave0
//ModifyGraph mode(maxconductance)=3,marker(maxconductance)=1,rgb(maxconductance)=(1,4,52428)
//ModifyGraph mrkThick(maxconductance)=2
ModifyGraph axisEnab(L2)={0.7,1}
ModifyGraph freePos(L2)=0
Label left "Pressure (cmH20)"
Label L2 "Conductance\r(microl/s/cmH2O)"
SetAxis L2 0,*
Label right "Flow rate (microl/s)"
ModifyGraph highTrip(bottom)=100000
SetAxis bottom wave0(refstart),*
Label bottom "Time (s)"
ModifyGraph lblPos(L2)=60
ModifyGraph lblMargin(right)=20
ModifyGraph textMarker(pressurepeaklevel)={eventnumber,"default",0,0,5,0.00,8.00}
ModifyGraph textMarker(maxdebiet)={voidnumber,"default",0,0,5,0.00,5.00}
AppendToGraph/R=R2 capacity vs voidstart
ModifyGraph axisEnab(R2)={0.7,1},freePos(R2)=0
ModifyGraph mode(capacity)=3,marker(capacity)=19,rgb(capacity)=(0,0,65280)
SetAxis R2 0,*
ModifyGraph lblPos(R2)=45;DelayUpdate
Label R2 "Capacity (microl)"



Display/N=scrollfigure/W=(930,0,1380,390) wave2 vs wave0
AppendToGraph/R debiet vs wave0
ModifyGraph rgb(debiet)=(1,4,52428)
AppendToGraph/R maxdebiet vs voidmaxtime
ModifyGraph mode(maxdebiet)=3,marker(maxdebiet)=1,mrkThick(maxdebiet)=1,rgb(maxdebiet)=(0,65535,0)
AppendToGraph ones vs voidstart
ModifyGraph mode(ones)=1,rgb(ones)=(16385,16388,65535)
AppendToGraph ones vs voidend
ModifyGraph mode(ones#1)=1, rgb(ones#1)=(16385,16388,65535)
AppendToGraph onesp vs pressurestart
ModifyGraph mode(onesp)=1
ModifyGraph lsize(onesp)=3
AppendToGraph onesp vs pressurestop
ModifyGraph mode(onesp#1)=1
Label bottom "Time (s)"
Label left "Pressure (cmH2O)"
Label right "Flow rate (microl/s)"
ModifyGraph highTrip(bottom)=100000
ModifyGraph msize(maxdebiet)=7;DelayUpdate
ModifyGraph textMarker(maxdebiet)={voidnumber,"default",1,0,5,0.00,2.00}

variable/G scrollnumber = 0
//SetAxis/W=scrollfigure bottom pressurestart[0]-20, pressurestop[0]+20
adjustxaxis("",0,"","")

SetAxis left 0,wavemax(wave2)
setaxis right -wavemax(debiet)/5, wavemax(debiet)

setvariable scroller size={120,20}, limits={0,(numpnts(pressurestart)-1),1}, pos={100,10}, title = "Event #", fsize=16
setvariable scroller value=scrollnumber, proc = adjustxaxis

button deletebadevent size={120,20}, pos={360,10}, title = "Delete event", fsize=16, proc = deleteevent






end

macro pressureflowcurve(startvoid, endvoid, prepoints, postpoints)
variable startvoid=0, endvoid=1, prepoints=100, postpoints=100

make/O/N=((prepoints+postpoints+1),(endvoid-startvoid+1)) meanpressurem=0, meanflowm=0

variable count=startvoid

string meanpressurename = "meanpressure" + num2str(startvoid)+"_"+num2str(endvoid), meanflowname="meanflow"+num2str(startvoid)+"_"+num2str(endvoid)
string SEMpressurename = "SEMpressure" + num2str(startvoid)+"_"+num2str(endvoid), SEMflowname="SEMflow"+num2str(startvoid)+"_"+num2str(endvoid)
do
duplicate/O/R=[(voidmaxx[count]-prepoints),(voidmaxx[count]+postpoints)] debiet, cutoutflow
duplicate/O/R=[(voidmaxx[count]-prepoints),(voidmaxx[count]+postpoints)] wave2, cutoutpressure
//print count
meanpressurem[][count-startvoid]=cutoutpressure[p]
meanflowm[][count-startvoid] =cutoutflow[p]
count +=1
while(count<(endvoid+1))


matrixop/O meanpressure = (averagecols(meanpressurem^t))^t
matrixop/O meanflow = (averagecols(meanflowm^t))^t

matrixop/O varpressure = sqrt((varcols(meanpressurem^t))^t)/sqrt(endvoid-startvoid+1)
matrixop/O varflow = sqrt((varcols(meanflowm^t))^t)/sqrt(endvoid-startvoid+1)


//matrixop/O meanflow = sumrows(meanflowm)/(endvoid-startvoid+1)licate/O meanpressure, $meanpressurename
duplicate/O meanflow, $meanflowname
duplicate/O meanpressure, $meanpressurename
SetScale/P x 0,(wave0[2]-wave0[1]),"", $meanflowname, $meanpressurename


duplicate/O varpressure, $SEMpressurename
duplicate/O varflow, $SEMflowname
dowindow/K pressureflowwindow
display/N=pressureflowwindow $meanflowname vs $meanpressurename
//dowindow/K pressureflowwindow2
//display/N=pressureflowwindow2 $meanflowname, $meanpressurename
ModifyGraph lsize=3,rgb=(1,4,52428);DelayUpdate
ErrorBars/RGB=(56797,56797,56797) $meanflowname XY,wave=($SEMpressurename,$SEMpressurename),wave=($SEMflowname,$SEMflowname)
ErrorBars/T=0/X=1/Y=1 $meanflowname XY,wave=($SEMpressurename,$SEMpressurename),wave=($SEMflowname,$SEMflowname)
Label bottom "Pressure (cmH2O)"
Label left "Flow rate (μl/s)"
movewindow/W=pressureflowwindow 100,100,400,300


dowindow/K meanpressureflowwindow
display/N=meanpressureflowwindow $meanflowname
AppendToGraph/R  $meanpressurename
ModifyGraph rgb($meanpressurename)=(1,4,52428)
Label bottom "Time (s)"
Label left "Flow rate (μl/s)"
Label right "Pressure (cmH2O)"
movewindow/W=meanpressureflowwindow 100,330,400,630


end


macro adjustxaxis(ctrlName, varnum, varstr, varname) : setvariablecontrol
string ctrlname
variable varnum
string varstr
string varname
variable lx, rx

SetAxis/W=scrollfigure bottom pressurestart[scrollnumber]-20, pressurestop[scrollnumber]+20



findlevel/Q wave0,pressurestart[scrollnumber]-20
lx=V_levelX
findlevel/Q wave0,pressurestop[scrollnumber]+20
rx = V_levelx
duplicate/O/R=[lx,rx] wave0,waart
duplicate/O/R=[lx,rx] wave2,waar

AppendToGraph/W=overviewfigure/L waar vs waart
ModifyGraph/W=overviewfigure lsize(waar)=1,rgb(waar)=(65280,0,0)


end




macro deleteevent(test) : buttoncontrol
string test

if (nonvoiding[scrollnumber]==0)
deletepoints voidNo[scrollnumber],1,voidmaxx,voidmaxtime, voidstart, voidend, maxconductance, maxdebiet, maxpressure, startpressure, endpressure, voidvolume, voidefficiency, voidduration, voidnumber, ici, residualvolume, capacity
voidNo[scrollnumber,] -=1
voidnumber[] = p
endif

deletepoints scrollnumber,1,  pressurepeak,pressurepeaklevel, pressurestart, pressurestop, pressureamplitude, voidNo, nonvoiding, eventnumber
eventnumber[]=p
setvariable scroller limits={0,(numpnts(pressurestart)-1),1}

adjustxaxis("",0,"","")

end


Function WMAxisSliderInGraph(grfName)
	String grfName	// or ""
	
	if( strlen(grfName) == 0 )
		grfName= WinName(0, 1)
		if( strlen(grfName) == 0 )
			return 0
		endif
	endif
	
	DFREF dfr= root:Packages:WMAxisSlider:$(grfName)
	if( DataFolderRefStatus(dfr) == 0 )
		return 0
	endif
	ControlInfo/W=$grfName WMAxSlSl
	return V_Flag == 7 // returns 1 if Slider is in the graph.
End

Function/S WMAxisSliderExpandMenuItem()
	
	String item= "" // disappears
	if( WMAxisSliderInGraph("") )
		item= "-;Zoom to selection;"
	endif
	return item
End

Function WMAxisSliderZoomToMarquee()

	String dfSav= GetDataFolder(1)
	String grfName= WinName(0, 1)
	SetDataFolder root:Packages:WMAxisSlider:$(grfName)

	NVAR gLeftLim,gRightLim
	SVAR gAxisName
	//GetAxis/Q $gAxisName
	GetMarquee/W=$grfName/K/Z $gAxisName
	SetAxis $gAxisName,V_left,V_right
	DoUpdate/W=$gAxisName

	// set zoom and resync slider
	WMAxisSliderSetAxis(gAxisName,0,gLeftLim,gRightLim)
	SetDataFolder dfSav
End


Function WMAxisSliderProc(name, value, event)
	String name	// name of this slider control
	Variable value	// value of slider
	Variable event	// bit field: bit 0: value set; 1: mouse down, //   2: mouse up, 3: mouse moved

	String dfSav= GetDataFolder(1)
	String grfName= WinName(0, 1)
	SetDataFolder root:Packages:WMAxisSlider:$(grfName)

	NVAR gLeftLim,gRightLim
	SVAR gAxisName
	GetAxis/Q $gAxisName
	//print gAxisname
	Variable dx= (V_max-V_min)/2
	Variable x0= value*(gRightLim-gLeftLim)+gLeftLim
	//print gLeftlim, grightlim, dx, x0, value, trunc(x0)
	SetAxis $gAxisName,x0-dx,x0+dx
	 

    
	
	
	SetDataFolder dfSav
	updateframebis(trunc(x0))
			
	return 0	// other return values reserved
End

static constant kSliderLMargin= 70
static constant ControlBarDelta = 46

#if Exists("PanelResolution") != 3
Static Function PanelResolution(wName)			// For compatibility with Igor 7
	String wName
	return 72
End
#endif

Function WMAppendAxisSlider()
	String grfName= WinName(0, 1)
	DoWindow/F $grfName
	if( V_Flag==0 )
		return 0			// no top graph, exit
	endif
	ControlInfo WMAxSlSl
	if( V_Flag != 0 )
		return 0			// already installed, do nothing
	endif
	String dfSav= GetDataFolder(1)
	NewDataFolder/S/O root:Packages
	NewDataFolder/S/O WMAxisSlider
	NewDataFolder/S/O $grfName
	
	Variable/G gLeftLim,gRightLim
	String/G gAxisName=""

	String CBIdent=""
	Variable/G gOriginalHeight = ExtendControlBar(grfName, ControlBarDelta, CBIdent) // we append below original controls (if any)
	String/G CBIdentifier = CBIdent
	
	// Find first axis of type bottom or top
	String axList= AxisList("" )
	Variable i,nax= ItemsInList(axList )
	for(i=0;i<nax;i+=1)
		gAxisName= StringFromList(i,axList)
		String axtype= StringByKey("AXTYPE",AxisInfo("", gAxisName))
		if( CmpStr(axtype,"bottom")==0 || CmpStr(axtype,"top")==0 )
			break
		endif
	endfor

	GetWindow kwTopWin,gsize	// returns points, and controls are positioned in pixels
	V_left *= ScreenResolution/PanelResolution(grfName)	// v1.02 - convert points to pixels
	V_right *= ScreenResolution/PanelResolution(grfName)
	
	Slider WMAxSlSl,pos={V_left+80,gOriginalHeight+9},size={V_right-V_left-kSliderLMargin,16},proc=WMAxisSliderProc
	Slider WMAxSlSl,limits={0,1,0},value= .5,vert= 0,ticks= 0,side=0
	PopupMenu WMAxSlPop,pos={V_left+10,gOriginalHeight+5},size={30,5},proc=WMAxSlPopProc
	PopupMenu WMAxSlPop,mode=0,title = "Zoom", value= #"\"Zoom Factor...;\""

	Variable/G gLastAuto=1
	WMAxisSliderSetAxis(gAxisName,gLastAuto,0,0)
	
	SetDataFolder dfSav
End

Function WMAxisSliderSetAxis(axName,doAutoscale,vmin,vmax)
	String axName
	Variable doAutoscale,vmin,vmax

	NVAR gLeftLim,gRightLim
	SVAR gAxisName

	gAxisName= axName
	GetAxis/Q $axName
	Variable origV_min= V_min,origV_max=V_max
	if( doAutoscale )
		SetAxis/A $axName
		DoUpdate
		GetAxis/Q $axName
	else
		V_min= vmin
		V_max= vmax
	endif
	gLeftLim= V_min
	gRightLim= V_max
	// JW 180530 Added doAutoscale check so that this function doesn't cause infinite updates if called from
	// a window hook modified event. The gLeftLim and gRightLim globals keep the limits of the axis if it were
	// autoscaled so that the slider has the right limit values. Needed if the data changes its X range.
	if (doAutoscale)
		SetAxis $axName,origV_min,origV_max
	endif

	Variable value= (((origV_max+origV_min)/2)-gLeftLim)/(gRightLim-gLeftLim)
	Slider WMAxSlSl,value=value
End

Function WMAxSlPopProc(ctrlName,popNum,popStr) : PopupMenuControl
	String ctrlName
	Variable popNum
	String popStr

	String dfSav,grfName
	
	StrSwitch(popStr)
		case "Set Axis...":
			dfSav= GetDataFolder(1)
			grfName= WinName(0, 1)
			SetDataFolder root:Packages:WMAxisSlider:$(grfName)
			
			String axList= AxisList("" )
			NVAR gLeftLim,gRightLim,gLastAuto
			SVAR gAxisName
			
			String newAx= gAxisName
			Prompt newAx,"axis",popup,axList
			Variable doAuto= gLastAuto ? 1 : 2
			Prompt doAuto,"autoscale",popup,"Yes;No"
			Variable axMin= gLeftLim
			Prompt axMin,"axis minimum"
			Variable axMax= gRightLim
			Prompt axMax,"axis maximum"
			
			DoPrompt "Set Axis for Slider",newAx,doAuto,axMin,axMax
			if( V_Flag==0 )
				gLastAuto= doAuto==1
				WMAxisSliderSetAxis(newAx,gLastAuto,axMin,axMax)
			endif

			SetDataFolder dfSav
			break
		case "Zoom Factor...":
			dfSav= GetDataFolder(1)
			grfName= WinName(0, 1)
			SetDataFolder root:Packages:WMAxisSlider:$(grfName)
			
			NVAR gLeftLim,gRightLim
			SVAR gAxisName
			GetAxis/Q $gAxisName
			// JP 7.09: added presetting of popup value to current zoom.
			Variable displayedRange= V_Max-V_Min
			Variable maxRange= gRightLim-gLeftLim
			Variable oldZoom= round(maxRange/displayedRange)
			String zoomFactor= num2istr(oldZoom) 				// JP 7.09: was = "10"
			String choices= "1;4;10;40;100;400;1000;"		// AL 1.03: Added 1 and 4 as options.
			Variable inList = WhichListItem(zoomFactor, choices) >= 0
			if( !inList )
				zoomFactor= "Other"
			endif
			Prompt zoomFactor,"Zoom Factor",popup,choices+"Other;"
			Variable other = oldZoom
			Prompt other, "Other Zoom Factor"
			DoPrompt "Set Axis Zoom",zoomFactor, other
			if( V_Flag==0 )
				Variable x0= (V_Max+V_Min)/2
				Variable zoom= str2num(zoomFactor)
				if( numtype(zoom) != 0 )
					zoom = other
				endif
				Variable dx= (gRightLim-gLeftLim)/(2*zoom)
				SetAxis $gAxisName,x0-dx,x0+dx
			endif

			SetDataFolder dfSav
			break
		case "Resync position":
			dfSav= GetDataFolder(1)
			grfName= WinName(0, 1)
			SetDataFolder root:Packages:WMAxisSlider:$(grfName)
			
			NVAR gLeftLim,gRightLim
			SVAR gAxisName
			WMAxisSliderSetAxis(gAxisName,0,gLeftLim,gRightLim)

			SetDataFolder dfSav
			break
		case "Resize":
			GetWindow kwTopWin,gsize	// returns points, and controls are positioned in pixels
			grfName= WinName(0, 1)
			V_left *= ScreenResolution/PanelResolution(grfName)	// v1.02 - convert points to pixels
			V_right *= ScreenResolution/PanelResolution(grfName)
			Slider WMAxSlSl,size={V_right-V_left-kSliderLMargin,16}
			break
		case "Instructions...":
			Execute/Q/Z "WMAxisSliderInstructions()"
			break
		case "Remove":
			dfSav= GetDataFolder(1)
			grfName= WinName(0, 1)
			SetDataFolder root:Packages:WMAxisSlider:$(grfName)
//			NVAR gOriginalHeight
//			Variable searchTop = gOriginalHeight+ControlBarDelta
//			String moveCList = ListControlsInControlBar(grfName, searchTop)
			KillControl WMAxSlSl
			KillControl WMAxSlPop
//			ControlBar gOriginalHeight
//			MoveControls(moveCList, 0, -ControlBarDelta)	
			SVAR CBIdentifier
			ContractControlBar(grfName, CBIdentifier, ControlBarDelta)
			KillDataFolder :
			if( CountObjects(":",4) == 0 )
				KillDataFolder :
				Execute/P "DELETEINCLUDE  <AxisSlider>"
				Execute/P "COMPILEPROCEDURES "
			endif
			SetDataFolder dfSav
			break
	endswitch
End


Function/c centreofMass2D(w)
wave w
//w is a 2D image wave.
//xx and yy are abscissa and ordinates of image wave
//xx and yy should be 1 point larger than w, as we are considering an image wave
//return C-O-M via a complex value

make/O/N=(dimsize(w,0)+1) xx
make/O/N=(dimsize(w,1)+1) yy
xx=p
yy=p
variable ii,jj,totalmass=0,sumrmx=0,sumrmy=0

for(ii=0;ii<dimsize(w,0);ii+=1)
    for(jj=0;jj<dimsize(w,1);jj+=1)
        totalmass+=w[ii][jj]
        sumrmx += w[ii][jj]*(xx[ii]+xx[ii+1])/2
        sumrmy+=w[ii][jj]*(yy[jj]+yy[jj+1])/2
    endfor
endfor
return cmplx(sumrmx/totalmass,sumrmy/totalmass)

End





function meanspectrum(spectrum, voidtimes, pts)
wave spectrum, voidtimes 
variable pts

variable numvoids = numpnts(voidtimes)-1

make/O/N=(pts, dimsize(spectrum,1), numvoids) tnormalizedspectra
variable count=0, count2=0


do

duplicate/O/R=[voidtimes(count), voidtimes(count+1)] spectrum, spectrumpiecespecttp
count2=0

	
count+=1
while (count<numvoids)

matrixOP/O normalizedspectrum=sumbeams(tnormalizedspectra)/numvoids

end


function xinterpolate (image,yp, pts)
wave image
variable yp, pts



matrixOP/O fr=col(image,yp)
Interpolate2/T=2/N=(pts)/E=2/Y=fri fr
end


function findvoidsfromborder(w,threshold, startx, endx)
wave w
variable threshold, startx, endx

findlevels/DEST=plasjes/M=20/R=(startx, endx) w, threshold

return V_levelsfound

end

macro voidspectrum (test): buttoncontrol
string test

string borderdifname, borderdifnamez, FFplotnamez


make/O/N=2 profilex,profileyt, profileyb
profilex={0,100}
profileyt=spectb+spectw/2
profileyb=spectb-spectw/2

findvoidsfromborder(bordermove, peethreshold, startsa, endsa)


make/O/N=(dimsize(movie2, 2), 65, numpizzaslices) spectrum3Dz
variable count=0, countvoids
string spectrumname, meanspectrumname
do
borderdifname = "borderpointdif"+num2str(count)
borderdifnamez= "borderpointdifz"+num2str(count)
duplicate/O $borderdifname, $borderdifnamez
	countvoids=0
	do
	//$borderdifnamez[round(plasjes[countvoids])-15*framerate, round(plasjes[countvoids])+15*framerate]=$borderdifnamez[p+30*framerate]
	$borderdifnamez[round(plasjes[countvoids])-15*framerate, round(plasjes[countvoids])+15*framerate]=0
	
	countvoids+=1
	while(countvoids<numpnts(plasjes))

FFplotnamez = "FFplotz"+num2str(count)

duplicate/O $borderdifnamez movementtrace
setscale/P x 0,1/framerate, movementtrace


STFT/SEGS=64/OUT=6/HOPS=1/WINF=Cos1/DEST=$FFplotnamez movementtrace
SetScale/P x 0,1,"Frame", $FFplotnamez
spectrum3Dz[][][count] = $FFplotnamez[p][q]
setscale/P y 0,dimdelta($FFplotnamez,1), "Hz", spectrum3Dz
meanspectrumname="FFmean"+num2str(count)
spectrumname = "FFplotz"+num2str(count)
//meanspectrum($spectrumname, plasjes, 100)
//duplicate/O normalizedspectrum, $meanspectrumname
//SetScale/P y 0,dimdelta($FFplotnamez,1),"Hz", $meanspectrumname

makebetteraveragespectrum($borderdifname, plasjes, 20, 20, 100)
duplicate/O meanspectrum2, $meanspectrumname
SetScale/P y 0,dimdelta($FFplotnamez,1),"Hz", $meanspectrumname
count+=1
while (count<numpizzaslices)

meanspectrumname="FFmean"+num2str(bordershown)
duplicate/O $meanspectrumname, normspec




RemoveImage/Z normspec
AppendImage/R=middley/B=middlex normspec

ModifyGraph axisEnab(middley)={0.45,0.75},axisEnab(middlex)={0.45,0.65}
ModifyImage normspec ctab= {*,maxspectrum,ColdWarm,0}
ModifyGraph noLabel(middlex)=2,axThick(middlex)=0

appendtograph/R=middley/B=middlex profileyt vs profilex
appendtograph/R=middley/B=middlex profileyb vs profilex
Label middley "\\Z10"
ModifyGraph btLen(middley)=2,freePos(middley)={0.35,kwFraction}
ModifyGraph fSize(middley)=12
button voidmovie disable=0

duplicate/O plasjes, plasjes0
plasjes0=0
removefromgraph/Z plasjes0
AppendToGraph/L=topleft/T plasjes0 vs plasjes
ModifyGraph mode(plasjes0)=3,marker(plasjes0)=25,rgb(plasjes0)=(3,52428,1)

end


function makebetterspectrum (movementtrace, voidtiming, voidnumber, postcut, precut, numpoints)
wave movementtrace, voidtiming
variable voidnumber, postcut, precut, numpoints


duplicate/R=[voidtiming[voidnumber]+postcut, voidtiming[voidnumber+1]-precut]/O movementtrace, cutoutmt

Interpolate2/T=1/N=200/Y=cutoutmt_L cutoutmt
duplicate/O cutoutmt_L, cutoutmt

STFT/SEGS=64/OUT=6/HOPS=1/WINF=Cos1/Dest=specttemp cutoutmt

duplicate/O specttemp, spectrumcutout


end

function makebetteraveragespectrum(movementtrace, voidtiming, postcut, precut, numpoints)
wave movementtrace, voidtiming
variable postcut, precut, numpoints

makebetterspectrum(movementtrace, voidtiming, 0, postcut, precut, numpoints)
wave spectrumcutout
spectrumcutout=0
duplicate/O spectrumcutout, meanspectrum2


variable count=0
do

makebetterspectrum(movementtrace, voidtiming, count, postcut, precut, numpoints)
print count, dimsize(spectrumcutout,0), dimsize(spectrumcutout,1)

meanspectrum2 += spectrumcutout

print "done"
count+=1

while (count<dimsize(voidtiming,0)-2)

meanspectrum2 /=(dimsize(voidtiming,0)-1)
SetScale/P x 1,1,"", meanspectrum2

end



macro profileline(test): buttoncontrol
string test

makeprofile(normspec, spectb,spectw)

spectrumband =W_imagelineprofile

end


function makeprofile(spectrum, position, width)
wave spectrum
variable position, width

make/O/N=2 profilex, profiley
profiley=position
profilex[0]=0
profilex[1]=dimsize(spectrum,0)
print position, width, scaletoindex(spectrum,width,1)
imagelineProfile xwave = profilex, ywave=profiley, srcwave=spectrum, width=scaletoindex(spectrum,width,1)
end


function averagevoid(voidtiming, movie, width)
wave voidtiming, movie
variable width
variable count=0
variable numberofvoids = numpnts(voidtiming)-1

mAKE/O/D/N=(dimsize(movie,0), dimsize(movie,1), 2*width+1) voidmovie=0

do

duplicate/O/R=[][][(round(voidtiming[count]))-width,(round(voidtiming[count]))+width] movie tee

voidmovie +=tee[p][q][r]/numberofvoids
count +=1
//print count
while(count<numberofvoids)

end

macro makevoidmovie(test): buttoncontrol
string test

averagevoid(plasjes, movie2,50)

killwindow/Z thevoidmovie
Display/N=thevoidmovie;DelayUpdate
AppendImage voidmovie
ModifyImage voidmovie ctab= {*,*,ColdWarm,0}
WMAppend3DImageSlider()
end


function xytoimage(xwave,ywave, zwave, xdim, ydim)
wave xwave, ywave, zwave
variable xdim, ydim

make/O/N=(xdim,ydim) imagexy=0

variable count=0
do
imagexy[xwave[count]][ywave[count]] = zwave[count]*1000
count+=1
while (count<numpnts(xwave))
end


macro makebordermovie(test): buttoncontrol
string test

silent 1; pauseupdate
make/O/N=(dimsize(movie2,0), dimsize(movie2,1), dimsize(movie2,2)) bordermovie
bordermovie = NaN
variable count=0

do

makexy(newborder,count)
xytoimage(xvalues,yvalues,zvalues,  dimsize(movie2,0), dimsize(movie2,1))
bordermovie[][][count] = imagexy[p][q]
count+=1
while(count< dimsize(movie2,2))
end
