#pragma rtGlobals=3		// Use modern global access method and strict wave access.
#pragma rtGlobals=3		// Use modern global access method and strict wave access.

	


macro offtimes2()

	variable j,k
	string	waver, waverh, waver2,
	
	
	make /o /n=0	taunums
	make /o /n=0	wavenam
	make /o /n=0	taus
	make /o /n=0	means
edit
	
	make /n=50/o waverhist
	k=0
	j=0
	do
		waver= "wave"+ num2istr(j)
		waver2= "wave"+ num2istr(j)+"b"
		waverh="wave"+ num2istr(j)+"h"
		
		duplicate /o $waver $waver2
		
		sort $waver2, $waver2
		findlevel /q $waveR2, 0.3  // no values <2x GCaMP6f timeconstant
		if (V_Flag==0)
			deletepoints 0, V_LevelX+1, $waver2
		endif
		wavestats /q $waver2
		if(V_npnts>=61)
			make /o /n=50 $waverh
			Histogram/B={0,3,50} $waver2,$waverh
			insertpoints k,1, taunums
			taunums[k]=V_npnts
			insertpoints k,1, wavenam
			wavenam[k]=j
			wavestats /q $waver2
			insertpoints k,1, means
			means[k]=V_avg
			CurveFit/q/M=2/W=0 exp, $waverh
			insertpoints k,1, taus
			taus[k]=1/W_coef[2]
			append $waver2
			
			k+=1
		endif
		
		if (j==0)
			duplicate /o $waver wg
		else
		 	concatenate /np {$waver}, wg
		endif
		
	
		j+=1
			
	while (j<491)
	
	wavestats taus
	wavestats means

end
