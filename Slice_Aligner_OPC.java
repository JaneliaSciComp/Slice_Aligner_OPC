//slice aligner based on a template slice
//written by Hideo Otsuna (HHMI.janelia, otsunah@janelia.hhmi.org) 2017 July 10

//Copyright (c) 2017, Howard Hughes Medical Institute, All rights reserved.

//Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

//Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
//Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
//Neither the name of the Howard Hughes Medical Institute nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.
//THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, ANY IMPLIED WARRANTIES OF MERCHANTABILITY, NON-INFRINGEMENT, OR FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; REASONABLE ROYALTIES; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

import ij.*;
import ij.plugin.filter.*;
import ij.plugin.PlugIn;
import ij.process.*;
import ij.gui.*;
import java.math.*;
import java.io.*;
import java.util.*;
import java.net.*;
import ij.Macro.*;
import java.awt.*;
import java.util.concurrent.atomic.AtomicInteger;
import ij.gui.GenericDialog.*;
import java.math.BigDecimal;
import java.math.BigInteger;



public class Slice_Aligner_OPC implements PlugInFilter
{
	ImagePlus isamp, itemp,dupsampleimp3,newimp1,newimp2,newimp3,isampOri;
	ImageProcessor ip2, ipTemp, ipnew,ipSampDup3,ipSignal1,ipSignal2,ipSignal1dup,ipSignal1rot,ipSignal2dup,ipSignal2rot;
	ImageStack Signal1, Signal2;
	
	double maxvalue=0,MaxrotateP=0;
	
	int  Maxshiftx, Maxshifty,Maxshiftx2, Maxshifty2,MaxshiftxZ, MaxshiftyZ,tempSliceNo; 
	int secondTime=0,BaseChannel,thread_num_, CalM, sampleDominant,nslice;
	
	int[] values = new int[10];
	double[] Dvalues = new double[2];
	BigDecimal[] bigvalues = new BigDecimal[2];
	
	//	long SumGap;
	public int setup(String arg, ImagePlus isamp)
	{
		IJ.register (Slice_Aligner_OPC.class);
		if (IJ.versionLessThan("1.32c")){
			IJ.showMessage("Error", "Please Update ImageJ.");
			return 0;
		}
		
		//	IJ.log(" wList;"+String.valueOf(wList));
		
		this.isamp = isamp;
		if(isamp.getType()!=isamp.GRAY8 && isamp.getType()!=isamp.GRAY16){
			IJ.showMessage("Error", "Plugin requires 8- or 16-bit image");
			return 0;
		}
		return DOES_8G+DOES_16;
		
	}
	
	public void run(ImageProcessor ip){
		
		int wList [] = WindowManager.getIDList();
		if (wList==null) {
			IJ.showMessage("There must be at least one window open");
			return;
		}
		int imageno = 0;
		String titles [] = new String[wList.length];
		for (int i=0; i<wList.length; i++) {
			ImagePlus isamp = WindowManager.getImage(wList[i]);
			if (isamp!=null){
				titles[i] = isamp.getTitle();//Samp.tif and Data.tif
				imageno = imageno +1;
			}else
			titles[i] = "";
		}
		
		isampOri = WindowManager.getCurrentImage();
		isamp = WindowManager.getCurrentImage();
		String TitleName=isamp.getTitle();
		
		int nChannels = isamp.getNChannels();
		String Color [] = new String[3];
		Color[0]="C1";
		Color[1]="C2";
		//Color[2]="C3";
		Color[2]="All Channels";
		/////Dialog//////////////////////////////////////////////		
		
		String tempimgSelect=(String)Prefs.get("tempimgSelect.String","First slice");
		double poRot=(double)Prefs.get("poRot.double",3);
		double miRot=(double)Prefs.get("miRot.double",3);
		double Overlap=(int)Prefs.get("Overlap.double",90);
		thread_num_=(int)Prefs.get("thread_num.int",6);
		BaseChannel=(int)Prefs.get("BaseChannel.int",0);
		//		String BrightSelect=(String)Prefs.get("BrightSelect.String","none");
		if(BaseChannel==3)
		BaseChannel=0;
		
		//CalM = (int)Prefs.get("CalM.int",0);
		//	sampleDominant = (int)Prefs.get("sampleDominant.int",0);
		double RotIncri=(double)Prefs.get("RotIncri.double",0.2);
		boolean thresholding=(boolean)Prefs.get("thresholding.boolean",false);
		
		GenericDialog gd = new GenericDialog("Slice Aligner OPC");
		
		String tempimg [] = new String[3];
		tempimg[0]="Current slice"; tempimg[1]="First slice"; tempimg[2]="input slice number";
		
		//		String BrightST [] = new String[3];
		//		BrightST[0]="average"; BrightST[1]="median"; BrightST[2]="none";
		
		gd.addRadioButtonGroup("Template slice: ", tempimg, 1, 3, tempimgSelect);
		
		gd.addNumericField("+ rotaion angle", poRot,0);
		gd.addNumericField("- rotaion angle", miRot,0);
		
		gd.addNumericField("Overlap percentage (higher is faster, max 100)", Overlap,0);
		gd.addNumericField("Parallel Threads", thread_num_, 0);
		
		gd.addNumericField("Rotation increment", RotIncri, 1);
		
		if(nChannels!=1){
			gd.addChoice("Reference channel", Color, Color[BaseChannel]);
		}
		
		gd.addCheckbox("70% background subtraction", thresholding);
		//	gd.addRadioButtonGroup("Brightness equalization", BrightST, 1,3, BrightSelect);
		
		//	String []	CalMST = {"Subtraction", "ABS gap", "OBJ peasonCoeff"};
		//	gd.addRadioButtonGroup("Calculation method", CalMST, 0, 3, CalMST[CalM]);
		
		//	String []	CalMST2 = {"Equal weight (temp and sample)", "Sample dominant"};
		//	gd.addRadioButtonGroup("Weight method", CalMST2, 0, 2, CalMST2[sampleDominant]);
		
		gd.showDialog();
		if(gd.wasCanceled()){
			return;
		}
		
		tempimgSelect=gd.getNextRadioButton();
		
		poRot=(double)gd.getNextNumber();
		miRot=(double)gd.getNextNumber();
		Overlap=(double)gd.getNextNumber();
		thread_num_ = (int)gd.getNextNumber();
		RotIncri=(double)gd.getNextNumber();
		
		if(nChannels!=1){
			BaseChannel= gd.getNextChoiceIndex();
		}
		thresholding=(boolean)gd.getNextBoolean();
		
		//	BrightSelect=gd.getNextRadioButton();
		//	String STdominant = (String)gd.getNextRadioButton();
		
		IJ.log("tempimgSelect; "+tempimgSelect);
		if(tempimgSelect=="input slice number"){
			tempSliceNo=(int)Prefs.get("tempSliceNo.int",1);
			GenericDialog gd2 = new GenericDialog("Template slice number");
			gd2.addNumericField("Template slice number", tempSliceNo, 0);
			
			gd2.showDialog();
			if(gd2.wasCanceled()){
				return;
			}
			tempSliceNo = (int)gd2.getNextNumber();
			Prefs.set("tempSliceNo.int",tempSliceNo);
		}else if(tempimgSelect=="First slice")
		tempSliceNo=1;
		else if(tempimgSelect=="Current slice"){
			int nSlices = isamp.getNSlices();
			if(nSlices!=1)
			tempSliceNo=isamp.getSlice();
			else 
			tempSliceNo=isamp.getFrame();
		}
		
		
		
		
		if(RotIncri==0)
		RotIncri=1;
		final double RotIncrifin=RotIncri;
		Prefs.set("RotIncri.double",RotIncri);
		Prefs.set("thresholding.boolean",thresholding);
		Prefs.set("tempimgSelect.String",tempimgSelect);
		Prefs.set("BaseChannel.int",BaseChannel);
		//	Prefs.set("BrightSelect.String",BrightSelect);
		Prefs.set("poRot.double",poRot);
		Prefs.set("miRot.double",miRot);
		Prefs.set("Overlap.double",Overlap);
		
		//	if(STdominant == "Equal weight (temp and sample)")
		sampleDominant = 0;
		//	else if(STdominant == "Sample dominant")
		//	sampleDominant = 1;
		
		//	Prefs.set("sampleDominant.int", sampleDominant);
		
		
		if(Overlap>100){
			IJ.log(String.valueOf(Overlap)+"overlap %, must be less than 100%");
			return;
		}
		
		IJ.log("Template slice; "+String.valueOf(tempSliceNo));
		IJ.log("Positive Rot; "+String.valueOf(poRot));
		IJ.log("Negative Rot; "+String.valueOf(miRot));
		IJ.log("RotIncri; "+String.valueOf(RotIncri));
		IJ.log("Overlap; "+String.valueOf(Overlap)+" %");
		IJ.log("thread_num_; "+String.valueOf(thread_num_)+" CPU");
		IJ.log("70% thresholding; "+thresholding);
		IJ.log("nChannels; "+nChannels+" Ch");
		
		int neuronChannel1=0, neuronChannel2=0;
		String AllCh="No";
		if(nChannels!=1){
			
			IJ.log("BaseChannel; C"+String.valueOf(BaseChannel+1));
			if(nChannels==3){
				if(BaseChannel==0){
					neuronChannel1=1;
					neuronChannel2=2;
				}else if(BaseChannel==1){
					neuronChannel1=0;
					neuronChannel2=2;
				}else if(BaseChannel==2){
					neuronChannel1=0;
					neuronChannel2=1;
				}
			}
			
			if(nChannels==2){
				if(BaseChannel==0)//C1
				neuronChannel1=1;
				else if(BaseChannel==1)
				neuronChannel1=0;
				else if(BaseChannel==2){//"All Channels"
					nChannels=3;
					neuronChannel1=0;
					neuronChannel2=1;
					BaseChannel=1;
					AllCh="Yes";
					IJ.log("All channel mode "+nChannels+" channels");
				}
			}
		}//if(nChannels!=1){
		
			ImagePlus[] channels = new ImagePlus[nChannels];
			
		channels = splitChannels(isamp, nChannels);
		
		if(nChannels!=1)
		isamp=channels[BaseChannel];

		
		if(thread_num_ <= 0) thread_num_ = 1;
		Prefs.set("thread_num.int", thread_num_);
		
		int widthTemp = isamp.getWidth();
		int heightTemp = isamp.getHeight();
		
		ImageStack dcStackfinal = new ImageStack (widthTemp,heightTemp);
		ImageStack	dcStackfinal2 = new ImageStack (widthTemp,heightTemp);
		ImageStack	dcStackfinal3 = new ImageStack (widthTemp,heightTemp);
		
		// shifting value decision ////////////////////////
		double bb=(double)Overlap/100;
		double aa=bb*(double)widthTemp;
		int OverlapW=Math.round((int)aa);//50~100% seach 42 at 60
		
		double cc=bb*(double)heightTemp;
		int OverlapH=Math.round((int)cc);//50~100% seach
		
		final int MaxXshift=Math.round(widthTemp-OverlapW);// 18 FOR 60 at 70%
		final int MaxYshift=Math.round(heightTemp-OverlapH);// 20 for 100 and 70%, from 0 is 50% overlap
		
		int maxshift=0;
		if(MaxXshift>MaxYshift)
		maxshift=MaxXshift;
		else
		maxshift=MaxYshift;
		
		final int MaxShift=maxshift;
		
		/////////////////////Start: signal detection///////////////////////////////
		double OBJscore = 0;
		//	IJ.log("MaxXshift; "+String.valueOf(MaxXshift));
		//	IJ.log("MaxYshift; "+String.valueOf(MaxYshift));
		
		double fullgap=(double)MaxXshift*2; double incriGap=0;
		nslice = isamp.getNSlices();
		if(nslice==1)
		nslice = isamp.getNFrames();
		
		int widthSamp = isamp.getWidth();
		int heightSamp = isamp.getHeight();
		
		if(AllCh=="Yes" && nChannels!=1){
			IJ.log("Creating max ch...");
			channels[2]=channels[0];
			
			ImagePlus ch1i=channels[0];
			ImagePlus ch2i=channels[1];
			ImagePlus ch3i=channels[2];
			
			ImageStack ch1istack=ch1i.getStack();
			ImageStack ch2istack=ch2i.getStack();
			ImageStack ch3istack=ch3i.getStack();
			
			for(int iplus=1; iplus<=nslice; iplus++){//slice
				ImageProcessor ch1ip=ch1istack.getProcessor(iplus);
				ImageProcessor ch2ip=ch2istack.getProcessor(iplus);
				ImageProcessor ch3ip=ch3istack.getProcessor(iplus);
				
				for(int ipixAdd=0; ipixAdd<widthSamp*heightSamp; ipixAdd++){//image, pix
					int ch1pix=ch1ip.get(ipixAdd);
					int ch2pix=ch2ip.get(ipixAdd);
					
					int MaxPixAdd=0;
					if(ch1pix>ch2pix)
					MaxPixAdd=ch1pix;
					else
					MaxPixAdd=ch2pix;
					
					ch3ip.set(ipixAdd,MaxPixAdd);
				}
			}
			isamp=channels[2];
			IJ.log("Max ch creation done");
		}//	if(AllCh=="Yes" && nChannels!=1){
		
		
		
		ImageProcessor tempo = isamp.getProcessor(); 
		int sumpx = tempo.getPixelCount();
		double sumMaxLowBri=0;
		
		ImagePlus DUP =isamp.duplicate();
		
		ImageStack stackOri = DUP.getStack();
		int[] Detected_Thre = new int[nslice+1];
		int [] pureGap = new int[nslice+1];
		int [] maxPIXarray = new int[nslice+1];
		
		if(thresholding){
			int maxPIX=0;
			double maxi=0;
			for(int islicen=1; islicen<=nslice; islicen++){
				
				maxPIX=0;
				maxi=0;
				
				ImageProcessor ipOri = stackOri.getProcessor(islicen); 
				
				int[] HistoG = new int[65536];
				int ioripx=0, maxcounts=0;
				
				for(int pxscan=0; pxscan<sumpx; pxscan++){// histogram creation
					ioripx= ipOri.get(pxscan);
					
					int Grayval=HistoG[ioripx];
					HistoG[ioripx]=Grayval+1;
					
					if(ioripx>maxPIX)
					maxPIX=ioripx;
				}
				
				for(int i3=5; i3<maxPIX; i3++){// shifting histogram, max amount value decision
					
					int sumVal20=0; 
					if(i3<maxPIX-30){
						for(int aveval=i3; aveval<i3+30; aveval++){
							int Val20=HistoG[aveval];
							
							sumVal20=sumVal20+Val20;
						}
						int AveVal20=sumVal20/30;
						
						if(AveVal20>maxcounts){
							maxcounts=AveVal20;
							maxi= (double) (i3+15);
						}
					}//if(i3<280){
				}
				
				Detected_Thre[islicen] = (int) (maxi* 0.7);//clip 70% of dimmer signals
				
				sumMaxLowBri=sumMaxLowBri+(maxi* 0.7);//
				pureGap[islicen] =maxPIX- (int) maxi;
				maxPIXarray[islicen]=maxPIX;
				//	IJ.log("pureGap; "+String.valueOf(pureGap[islicen])+"   maxPIX; "+String.valueOf(maxPIX)+"   maxi; "+String.valueOf(maxi));
				
			}//for(int islicen=1; islicen<=nslice; islicen++){
			int aveLowbri=(int) (sumMaxLowBri/nslice);
			
			for(int islicen2=1; islicen2<=nslice; islicen2++){
				for(int pxscan2=0; pxscan2<sumpx; pxscan2++){
					ImageProcessor ipOri = stackOri.getProcessor(islicen2);
					int ioripx= ipOri.get(pxscan2);
					
					if(ioripx<aveLowbri)
					ipOri.set(pxscan2,0);
					else{
						
						double percentage=0;
						if(Detected_Thre[islicen2]<aveLowbri*2)
						percentage= (double) (ioripx - Detected_Thre[islicen2]) / (double) pureGap[islicen2];
						else
						percentage= (double) (ioripx - aveLowbri) / (double) pureGap[islicen2];
						
						double finalval= percentage* (double) ioripx;
						int finalval2 = (int) Double.parseDouble(String.format("%.0f", finalval));
						
						if(finalval2>maxPIXarray[islicen2])
						finalval2=maxPIXarray[islicen2];
						else if(finalval2<0)
						finalval2=0;
						
						if(IJ.escapePressed())
						return;
						
						
						ipOri.set(pxscan2,finalval2);
						//	IJ.log("finalval2; "+String.valueOf(finalval2));
					}
				}//for(int pxscan2=0; pxscan<sumpx; pxscan++){
				//	IJ.log("islicen2; "+String.valueOf(islicen2)+"   Detected_Thre[islicen2]; "+String.valueOf(Detected_Thre[islicen2]));
			}//for(int islicen2=1; islicen2<=nslice; islicen2++){
			IJ.log("Nslices; "+String.valueOf(nslice)+"   aveLowbri; "+String.valueOf(aveLowbri));
		}//if(thresholding!=0){
		final ImageStack stack = stackOri;
		
		
		//	newimp1 = new ImagePlus(TitleName+"_70% subtracted.tif", stackOri);
		//	newimp1.setCalibration(isampOri.getCalibration());
		//	newimp1.show();
		
		ipTemp = stack.getProcessor(tempSliceNo); //Samp
		
		
		double AVEtempB=0;
		
		long sumVXtemp = 0;
		long sumpxT = ipTemp.getPixelCount();
		
		if(sampleDominant!=1){
			
			for(int Gn0=0; Gn0<sumpx; Gn0++){
				int pixTemp = ipTemp.get(Gn0);//template
				
				if(pixTemp!=0)
				sumVXtemp = sumVXtemp+pixTemp;//total sum temp
			}//for(int n=0; n<sumpx; n++){
		}//if(sampleDominant!=1){
		
		AVEtempB=sumVXtemp/sumpxT;
		
		final double [] Detected_OBJscore = new double [nslice+1];//40000000
		final int[] Detected_BestX = new int[nslice+1];
		final int[] Detected_BestY = new int[nslice+1];
		final double[] Detected_BestRot = new double[nslice+1];
		
		IJ.log("Maxshift; "+String.valueOf(maxshift)+" pixel");
		
		final AtomicInteger ai = new AtomicInteger(1);
		final Thread[] threads = newThreadArray();
		final int Fsumpx=sumpx;
		final int FwidthSamp=widthSamp;
		final int FwidthTemp=widthTemp;
		final ImageProcessor FipTemp=ipTemp;
		final int FheightSamp=heightSamp;
		final int FsampleDominant=sampleDominant;
		final double FAVEtempB=AVEtempB;
		final double FmiRot=miRot;
		final double FpoRot=poRot;
		final int Fnslice=nslice;
		//	final ImageStack FdcStackfinal =dcStackfinal;
		
		for (int ithread = 0; ithread < threads.length; ithread++) {
			threads[ithread] = new Thread() {
				
				{ setPriority(Thread.NORM_PRIORITY); }
				
				public void run() {
					
					for(int SliceShift=ai.getAndIncrement(); SliceShift<=Fnslice; SliceShift=ai.getAndIncrement()){
						
						ImageProcessor ipSamp;
						
						if(Fnslice>1){
							ipSamp = stack.getProcessor(SliceShift); //Temp
							
							IJ.showProgress (SliceShift, Fnslice);
							IJ.showStatus("Slice; "+SliceShift);
						}else
						ipSamp=isamp.getProcessor();
						
						if(IJ.escapePressed()){
							//		newimp = new ImagePlus("XY_shift_&_rotation_fixed.tif", FdcStackfinal);
							//		newimp.show();
							return;
						}
						
						double [] ArrayScore = new double [1];//40000000
						int[] BestX = new int[1];
						int[] BestY = new int[1];
						double MaxOBJ=0;
						double[] BestRot = new double[1];
						
						for(int shiftx=-MaxXshift; shiftx<=MaxXshift; shiftx++){// x pixel shift
							for(int shifty=-MaxYshift; shifty<=MaxYshift; shifty++){// y pixel shift
								//	IJ.log("Y shift; "+String.valueOf(shifty));
								
								for(double rotateP=-FmiRot; rotateP<=FpoRot; rotateP+=RotIncrifin){
									int rotNum=0;
									if(rotateP<0)
									rotNum=(360*10)+(int) (rotateP*10);
									else
									rotNum=(int) (rotateP*10);
									
									//		IJ.log("rotateP; "+String.valueOf(rotateP));
									
									ImageProcessor ipSampDup = ipSamp.duplicate();
									ImageProcessor ipSampDup2 =  ipSamp.duplicate();//abs-subtracted image for total brightness
									
									double Rd = Double.parseDouble(String.format("%.2f", rotateP));
									ipSampDup.rotate(Rd);
									
									
									for(int nn=0; nn<Fsumpx; nn++){
										ipSampDup2.set(nn, 0);
									}
									int posipx=0;
									for(int iSx=0; iSx<FwidthSamp; iSx++){
										for(int iSy=0; iSy<FheightSamp; iSy++){
											int pixsa=ipSampDup.get(iSx,iSy);
											
											if(iSx+shiftx>=0 && iSx+shiftx<FwidthSamp){
												if(iSy+shifty>=0 && iSy+shifty<FheightSamp){
													if(pixsa>0){
														posipx=1;
														ipSampDup2.set(iSx+shiftx, iSy+shifty, pixsa);
														
													}
												}
											}
										}
									}
									
									if(posipx==1){
										OBJpearson(Rd,shiftx,shifty,ipSampDup2,FipTemp,FsampleDominant,FAVEtempB,0,ArrayScore,BestX,BestY,BestRot);
										
									}//if(posipx==1){
								} //for(rotateP=Rot; rotateP<=EndRot; rotateP++){
							}//	for(int shifty=-MaxYshift; shifty<=MaxYshift; shifty++){// y pixel shift
						}//for(int shiftx=-MaxXshift; shiftx<=MaxXshift; shiftx++){// x pixel shift
						
						Detected_OBJscore[SliceShift]=ArrayScore[0];
						Detected_BestX[SliceShift]=BestX[0];
						Detected_BestY[SliceShift]=BestY[0];
						Detected_BestRot[SliceShift]=BestRot[0];
						double MaxOBJdouble = Double.parseDouble(String.format("%.4f", Detected_OBJscore[SliceShift]));
						
					}//for(int SliceShift=2; SliceShift<=nslice; SliceShift++){
			}};
		}//	for (int ithread = 0; ithread < threads.length; ithread++) {
		startAndJoin(threads);
		
		//	ImageProcessor ipSignal1,ipSignal2,ipSignal1dup,ipSignal1rot,ipSignal2dup,ipSignal2rot;
		
		ImageStack stackOri2 = isamp.getStack();
		if(nChannels!=1){
			Signal1=channels[neuronChannel1].getStack();
			
			//	channels[neuronChannel1].show();
			
			if(nChannels==3){
				Signal2=channels[neuronChannel2].getStack();
				//	channels[neuronChannel2].show();
			}
		}//if(nChannels!=1){
		
		
		for(int addslice=1; addslice<=nslice; addslice++){
			if(addslice!=tempSliceNo){//if(addslice==tempSliceNo){
				
				ImageProcessor ipSamp2= stackOri2.getProcessor(addslice);
				ImageProcessor ipSampDup2 = ipSamp2.duplicate();// single slice for final image
				ImageProcessor ipSampDup = ipSamp2.duplicate();// for rotation
				for(int nn=0; nn<Fsumpx; nn++){
					ipSampDup2.set(nn, 0);
				}
				ipSampDup.rotate(Detected_BestRot[addslice]);
				
				if(nChannels>1){
					ipSignal1 = Signal1.getProcessor(addslice);
					ipSignal1dup = ipSignal1.duplicate();// single slice for final image
					ipSignal1rot = ipSignal1.duplicate();// for rotation
					for(int nn=0; nn<Fsumpx; nn++){
						ipSignal1dup.set(nn, 0);
					}
					ipSignal1rot.rotate(Detected_BestRot[addslice]);
					
					if(nChannels==3){
						ipSignal2 = Signal2.getProcessor(addslice);
						ipSignal2dup = ipSignal2.duplicate();// single slice for final image
						ipSignal2rot = ipSignal2.duplicate();// for rotation
						for(int nn=0; nn<Fsumpx; nn++){
							ipSignal2dup.set(nn, 0);
						}
						ipSignal2rot.rotate(Detected_BestRot[addslice]);
					}
				}//if(nChannels>1){
				
				
				for(int iSx=0; iSx<FwidthSamp; iSx++){
					for(int iSy=0; iSy<FheightSamp; iSy++){
						int pixsa=ipSampDup.get(iSx,iSy);//after rotation, get pixel
						int pixsa1=0, pixsa2=0;
						
						if(nChannels>1){
							pixsa1=ipSignal1rot.get(iSx,iSy);//after rotation, get pixel
							if(nChannels==3)
							pixsa2=ipSignal2rot.get(iSx,iSy);//after rotation, get pixel
						}
						if(iSx+Detected_BestX[addslice]>=0 && iSx+Detected_BestX[addslice]<FwidthSamp){
							if(iSy+Detected_BestY[addslice]>=0 && iSy+Detected_BestY[addslice]<FheightSamp){
								ipSampDup2.set(iSx+Detected_BestX[addslice], iSy+Detected_BestY[addslice], pixsa);
								if(nChannels>1){
									ipSignal1dup.set(iSx+Detected_BestX[addslice], iSy+Detected_BestY[addslice], pixsa1);
									if(nChannels==3)
									ipSignal2dup.set(iSx+Detected_BestX[addslice], iSy+Detected_BestY[addslice], pixsa2);
								}
							}
						}
					}
				}
				
				double MaxOBJdouble = Double.parseDouble(String.format("%.4f", Detected_OBJscore[addslice]));
				
				dcStackfinal.addSlice(addslice+" X; "+Detected_BestX[addslice]+"   Y; "+Detected_BestY[addslice]+"   Rotate; "+Detected_BestRot[addslice]+"  OBJ; "+MaxOBJdouble, ipSampDup2);
				
				if(nChannels>1){
					dcStackfinal2.addSlice(addslice+" X; "+Detected_BestX[addslice]+"   Y; "+Detected_BestY[addslice]+"   Rotate; "+Detected_BestRot[addslice]+"  OBJ; "+MaxOBJdouble, ipSignal1dup);
					if(nChannels==3)
					dcStackfinal3.addSlice(addslice+" X; "+Detected_BestX[addslice]+"   Y; "+Detected_BestY[addslice]+"   Rotate; "+Detected_BestRot[addslice]+"  OBJ; "+MaxOBJdouble, ipSignal2dup);
				}
				
				IJ.log(String.valueOf(addslice)+";  shiftx;"+String.valueOf(Detected_BestX[addslice])+"  shifty;"+String.valueOf(Detected_BestY[addslice])+"  rotation;"+String.format("%.2f", Detected_BestRot[addslice])+"  OBJ score;"+String.valueOf(MaxOBJdouble));//"   lower_bri_val; "+String.valueOf(Detected_Thre[addslice])
			}	else if(addslice==tempSliceNo){
				
				ImageProcessor ipTemp2= stackOri2.getProcessor(tempSliceNo);
				dcStackfinal.addSlice("Template", ipTemp2);
				
				
				ImageProcessor FirstSlice1, FirstSlice2;
				if(nChannels!=1){
					Signal1=channels[neuronChannel1].getStack();
					FirstSlice1= Signal1.getProcessor(tempSliceNo);
					dcStackfinal2.addSlice("Template", FirstSlice1);
					
					if(nChannels==3){
						Signal2=channels[neuronChannel2].getStack();
						FirstSlice2= Signal2.getProcessor(tempSliceNo);
						dcStackfinal3.addSlice("Template", FirstSlice2);
					}
				}//if(nChannels!=1){
				IJ.log(String.valueOf(addslice)+" ; Template slice ");
				
			}	//else adding temp
		}//for(int addslice=2; addslice<=nslice; addslice++){
		
		int endchannel=0;
		if(nChannels!=1){
			if(AllCh=="No"){
				newimp1 = new ImagePlus(TitleName+"_shift_rotation_fixed.tif", dcStackfinal);
				
				if(nChannels>1){
					endchannel=2;
					newimp2 = new ImagePlus(TitleName+"_shift_rotation_fixed.tif", dcStackfinal2);
					
					if(nChannels==3){
						newimp3 = new ImagePlus(TitleName+"_shift_rotation_fixed.tif", dcStackfinal3);
						endchannel=3;
						channels[2]=newimp3;
					}
				}
			}//	if(AllCh=="No"){
			
			if(AllCh=="Yes"){
				newimp1 = new ImagePlus(TitleName+"_shift_rotation_fixed.tif", dcStackfinal2);
				
				if(nChannels>1){
					endchannel=2;
					newimp2 = new ImagePlus(TitleName+"_shift_rotation_fixed.tif", dcStackfinal3);
					channels[2]=null;
					if(nChannels==3){
						newimp3 = new ImagePlus(TitleName+"_shift_rotation_fixed.tif", dcStackfinal);
						//	endchannel=2;
						nChannels=2;
					}
				}
			}
			
			channels[0]=newimp1;
			channels[1]=newimp2;
			ImagePlus impCombined = combineChannels(channels, 0, endchannel,Detected_BestX, Detected_BestY, Detected_BestRot, Detected_OBJscore,tempSliceNo);
			ImagePlus newimp = new CompositeImage(impCombined);
			newimp.setDisplayMode(IJ.COMPOSITE);
			newimp.setCalibration(isamp.getCalibration());
			newimp.setDimensions(nChannels,isampOri.getNSlices(),isampOri.getNFrames());
			newimp.show();
		}else{
			
			newimp1 = new ImagePlus(TitleName+"_shift_rotation_fixed.tif", dcStackfinal);
			newimp1.setCalibration(isampOri.getCalibration());
			newimp1.show();
			
		}
		//			isamp.setTitle(String.valueOf(totalmax));
	} //public void run(ImageProcessor ip){
	
	private Thread[] newThreadArray() {
		int n_cpus = Runtime.getRuntime().availableProcessors();
		if (n_cpus > thread_num_) n_cpus = thread_num_;
		if (n_cpus <= 0) n_cpus = 1;
		return new Thread[n_cpus];
	}
	
	public static void startAndJoin(Thread[] threads)
	{
		for (int ithread = 0; ithread < threads.length; ++ithread)
		{
			threads[ithread].setPriority(Thread.NORM_PRIORITY);
			threads[ithread].start();
		}
		
		try
		{   
			for (int ithread = 0; ithread < threads.length; ++ithread)
			threads[ithread].join();
		} catch (InterruptedException ie)
		{
			throw new RuntimeException(ie);
		}
	}
	
	ImagePlus[] splitChannels(ImagePlus imp, int nChannels) {
		//	int nChannels = imp.getNChannels();
		
		int nSlices = imp.getNSlices();
		int nFrames = imp.getNFrames();
		
		ImageStack stack = imp.getStack();
		ImagePlus[] channels = new ImagePlus[nChannels];
		for (int c = 0; c < nChannels; ++c) {
			ImageStack channelStack = new ImageStack(imp.getWidth(), imp.getHeight());
			
			if(nSlices>1){
				for (int s = 0; s < nSlices; ++s)
				channelStack.addSlice(stack.getProcessor(imp.getStackIndex(c+1, s + 1, 1)));
			}else{
				for (int s = 0; s < nFrames; ++s)
				channelStack.addSlice(stack.getProcessor(imp.getStackIndex(c+1, 1, s+1)));
			}
			ImagePlus channelImp = new ImagePlus(imp.getTitle() + "-" + c, channelStack);
			channelImp.setCalibration(imp.getCalibration());
			channelImp.setDimensions(1, nSlices, 1);
			channels[c] = channelImp;
			//	channelImp.show();
		}
		return channels;
	}
	
	ImagePlus combineChannels(ImagePlus[] channels3, int startchannel, int endchannel, int Detected_BestX[], int Detected_BestY[],double Detected_BestRot[], double Detected_OBJscore[],int tempSliceNo2) {
		
		IJ.log("combineChannels; start"+String.valueOf(startchannel)+"   end; "+String.valueOf(endchannel));
		ImageStack stack = new ImageStack(channels3[0].getWidth(), channels3[0].getHeight());
		int nSlices = channels3[0].getNSlices();
		if(nSlices==1)
		nSlices = channels3[0].getNFrames();
		
		for (int s = 0; s < nSlices; ++s)
		for (int c = startchannel; c < endchannel; ++c){
			double MaxOBJdouble = Double.parseDouble(String.format("%.4f", Detected_OBJscore[s+1]));
			if(s+1!=tempSliceNo2)
			stack.addSlice("shiftx;"+String.valueOf(Detected_BestX[s+1])+"  shifty;"+String.valueOf(Detected_BestY[s+1])+"  rotation;"+String.format("%.2f", Detected_BestRot[s+1])+"  OBJ score;"+String.valueOf(MaxOBJdouble),channels3[c].getStack().getProcessor(s + 1));
			else
			stack.addSlice("Template slice",channels3[c].getStack().getProcessor(s + 1));
			
		}
		int channelsize= endchannel-startchannel;
		
		ImagePlus imp = new ImagePlus(channels3[0].getTitle().replaceAll("\\.[^.]*$", "-1-" + endchannel), stack);
		
		for (int c2 = startchannel; c2 < endchannel; ++c2)
		imp.setCalibration(channels3[c2].getCalibration());
		
		imp.setDimensions(channelsize, nSlices, 1);
		return imp;
	}
	
	public void OBJpearson(double rotatePF,int shiftxF, int shiftyF, ImageProcessor ipSampDupF,ImageProcessor ipTempF, int sampleDominantF, double AVEtemp, int AF, double ArrayScore[], int BestX[], int BestY[], double BestRot[]){	
		
		long sumVXsamp = 0;
		long sumVXtemp = 0;
		long TotalCountTemp = 0;
		long valSum = 0;
		
		long tempPowSum = 0;
		long sampPowSum = 0;
		
		double TempMinusAVE= 0;
		double SampMinusAVE= 0;
		
		long sumpx = ipTempF.getPixelCount();
		
		double AVEsamp = 0;
		
		/// Sum value measurement /////////////////////////////////////////
		if(sampleDominantF!=1){
			for(int Gn=0; Gn<sumpx; Gn++){
				
				int pixSamp= ipSampDupF.get(Gn);//sample
				
				if(pixSamp!=0)
				sumVXsamp = sumVXsamp+pixSamp;//total sum sample
			}//for(int n=0; n<sumpx; n++){
			
		}else{
			
			for(int Gn=0; Gn<sumpx; Gn++){
				int pixTemp = ipTempF.get(Gn);//template
				int pixSamp= ipSampDupF.get(Gn);//sample
				
				if(pixSamp!=0){
					sumVXsamp = sumVXsamp+pixSamp;//total sum sample
					
					if(pixTemp!=0)
					sumVXtemp = sumVXtemp+pixTemp;//total sum temp
					
					TotalCountTemp=TotalCountTemp+1;
				}
				
			}//for(int n=0; n<sumpx; n++){
			
			TotalCountTemp=sumpx;
			if(sumVXtemp!=0 && TotalCountTemp!=0)
			AVEtemp = sumVXtemp/TotalCountTemp;//template
		}
		
		if(sumVXsamp!=0 && TotalCountTemp!=0)
		AVEsamp = sumVXsamp/TotalCountTemp;//sample
		
		//		IJ.log("AVEsamp; "+String.valueOf(AVEsamp));
		//		IJ.log("AVEtemp; "+String.valueOf(AVEtemp)+"   Vol; "+String.valueOf(TotalCountTemp));
		
		
		if(AVEtemp!=0 && AVEsamp!=0){
			if(sampleDominantF==1){
				for(int singlepix=0; singlepix<sumpx; singlepix++){
					
					double sampPix=ipSampDupF.get(singlepix);//+sumpx
					
					if(sampPix!=0){
						double tempPix=ipTempF.get(singlepix);
						//		if(tempPix>1 || sampPix>1){
						
						TempMinusAVE=tempPix-AVEtemp;
						SampMinusAVE=sampPix-AVEsamp;
						
						tempPowSum=tempPowSum+(new Double(TempMinusAVE*TempMinusAVE)).longValue();
						sampPowSum=sampPowSum+(new Double(SampMinusAVE*SampMinusAVE)).longValue();
						
						valSum=valSum+(new Double(TempMinusAVE*SampMinusAVE)).longValue();
						//		}
					}//if(sampPix!=0){
				}//for(xx=0; xx<width; xx++){
			}else{//if(sampleDominantF!=1){
				for(int singlepix=0; singlepix<sumpx; singlepix++){
					
					double tempPix=ipTempF.get(singlepix);
					double sampPix=ipSampDupF.get(singlepix);//+sumpx
					
					TempMinusAVE=tempPix-AVEtemp;
					SampMinusAVE=sampPix-AVEsamp;
					//		sqrt(((aveX-Xmin1_Result)*(aveX-Xmin1_Result)
					tempPowSum=tempPowSum+(new Double(TempMinusAVE*TempMinusAVE)).longValue();
					sampPowSum=sampPowSum+(new Double(SampMinusAVE*SampMinusAVE)).longValue();
					
					valSum=valSum+(new Double(TempMinusAVE*SampMinusAVE)).longValue();
					
				}//for(xx=0; xx<width; xx++){
			}//if(sampleDominantF==1){
			
			
		}else if(AVEsamp==0){
			
			if(sampleDominantF==1){
				for(int singlepix=0; singlepix<sumpx; singlepix++){
					
					double sampPix=ipSampDupF.get(singlepix);//+sumpx
					
					if(sampPix>0){
						double tempPix=ipTempF.get(singlepix);
						
						if(tempPix>1 || sampPix>1){
							
							TempMinusAVE=tempPix-AVEtemp;
							SampMinusAVE=sampPix;
							
							tempPowSum=tempPowSum+(new Double(TempMinusAVE*TempMinusAVE)).longValue();
							sampPowSum=sampPowSum+(new Double(SampMinusAVE*SampMinusAVE)).longValue();
							
							valSum=valSum+(new Double(TempMinusAVE*SampMinusAVE)).longValue();
						}//if(tempPix>1 || sampPix>1){
					}//if(sampPix>0){
				}//for(int singlepix=0; singlepix<sumpx; singlepix++){
			}else{
				for(int singlepix=0; singlepix<sumpx; singlepix++){
					
					double tempPix=ipTempF.get(singlepix);
					double sampPix=ipSampDupF.get(singlepix);//+sumpx
					
					if(tempPix>1 || sampPix>1){
						
						TempMinusAVE=tempPix-AVEtemp;
						SampMinusAVE=sampPix;
						
						tempPowSum=tempPowSum+(new Double(TempMinusAVE*TempMinusAVE)).longValue();
						sampPowSum=sampPowSum+(new Double(SampMinusAVE*SampMinusAVE)).longValue();
						
						valSum=valSum+(new Double(TempMinusAVE*SampMinusAVE)).longValue();
					}
				}//for(xx=0; xx<width; xx++){
			}//if(sampleDominantF==1){
		}//	if(AVEtemp!=0 && AVEsamp!=0){
		
		
		double score = 0;
		BigDecimal Cross =  new BigDecimal("0.00");
		BigDecimal xsqrt =  new BigDecimal("0.00");
		BigDecimal two =  new BigDecimal("2.00");
		if(sampPowSum!=0 && tempPowSum!=0){
			Cross = BigDecimal.valueOf(sampPowSum).multiply(BigDecimal.valueOf(tempPowSum));//
			Cross = Cross.setScale(2, BigDecimal.ROUND_HALF_UP);
			
			BigDecimal x = new BigDecimal(Math.sqrt(Cross.doubleValue()));
			xsqrt = x.add(new BigDecimal(Cross.subtract(x.multiply(x)).doubleValue() / (x.doubleValue() * 2.0)));
			xsqrt = xsqrt.setScale(4, BigDecimal.ROUND_HALF_UP);
			
			
			//	double x = Math.sqrt(Cross.doubleValue());
			//	IJ.log(" x;"+String.valueOf(x));
			//	long xx= (new Double(x)).longValue();
			
			if(x.compareTo(BigDecimal.ZERO)!=0){
				//		long xsqrt = xx+((Cross.longValue()-(xx*xx))/(xx*2));			
				
				if(valSum!=0 && xsqrt.compareTo(BigDecimal.ZERO)!=0){
					score= (double) valSum / xsqrt.doubleValue();
					
					if(score>ArrayScore[0]){
						ArrayScore[0]=score;
						BestX[0]=shiftxF;
						BestY[0]=shiftyF;
						BestRot[0]=rotatePF;
					}
					//		IJ.log(" score;"+String.valueOf(score));
					//		if(score>0.4)
					//		IJ.log(" ArrayValueF;"+String.valueOf(ArrayValueF)+"   score; "+String.valueOf(ArrayScore[ArrayValueF])+"  shiftX; "+String.valueOf(shiftxF)+"  rotatePF; "+String.valueOf(rotatePF));
				}
			}//	if(xx!=0){
		}else{//if(sampPowSum!=0 && tempPowSum!=0){
			IJ.log("The value is 0");
		}
	}//public void OBJpearson(i
} //public class



























