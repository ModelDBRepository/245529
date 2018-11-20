#Python version, i.e. alternative of NRDpostAB
#in python, type ARGS="par1 par2,mol1 mol2,subdir/fileroot,sstart ssend" then execfile('neurord_analysis.py')
#DO NOT PUT ANY SPACES NEXT TO THE COMMAS, DO NOT USE TABS
#e.g. ARGS="Ca GaqGTP,Ca GaqGTP Ip3,../Repo/plc/Model_PLCassay,15 20"
#from outside python, type python neurord_analysis [par1 par2] [mol1 mol2]
#Assumes that molecule outputs are integers, and the hypens used ONLY for parameters
#Can process multiple parameter variations, but all files must use same morphology, and meshfile.  
#It will provide region averages (each spine, dendrite submembrane, cytosol) and if spatialaverage=1,
#will calculate an average of n segments along the dendrite, 
#or whatever structure name is specified in dend variable

import os
import numpy as np
from matplotlib import pyplot
from string import *
import sys  
import glob
import header_parse as hparse
import plot_utils as pu

#######################################################
#indicate the name of the injection spines if you want to exclude them
fake = 'FAKE'
#indicate the name of submembrane region for totaling molecules that are exclusively submembrane
#only relevant for tot_species calculation. this should be name of structure concatenated with sub
submembname='sub'
#Spatial average (=1 to process) only includes the structure "dend", and subdivides into bins:
spatialaverage=0
dend="dend"
bins=10
#how much info to print
prnvox=1
prnheader=0
prninfo=0
showss=0
#outputavg determines whether output files are written
outputavg=0
##change these endings depending on whether using neurord3.x:
meshend="*mesh.txt.out"
concend='conc.txt.out'
## or neurord2.x (uncomment these)
#meshend="*mesh.txt"
#concend='conc.txt'
#Example of how to total some molecule forms; turn off with tot_species={}
#tot_species={
#        "PKAtot":["PKA", "PKAcAMP2", "PKAcAMP4", "PKAr"],
#        "D1Rtot":["D1R","DaD1R", "GsD1R","DaD1RGs", "pDaD1RGs", "PKAcDaD1RGs"],
#        "pde10tot":["PDE10","pPDE10", "PDE10cAMP","pPDE10cAMP","PKAcPDE10", "PKAcPDE10cAMP"],
#        "Gitot":["Giabg","AChm4RGi","Gim4R", "GaiGTP", "GaiGDP", "ACGai", "ACGasGai", "ACGasGaiATP"],
#        "m4Rtot":["AChm4RGi","Gim4R", "m4R", "AChm4R"]}
tot_species={}

Avogadro=6.023e14 #to convert to nanoMoles
mol_per_nM_u3=Avogadro*1e-15

try:
	args = ARGS.split(",")
	print "ARGS =", ARGS, "commandline=", args
 	do_exit = False
except NameError: #NameError refers to an undefined variable (in this case ARGS)
	args = sys.argv[1:]
	print "commandline =", args
	do_exit = True

pattern=args[2]+'*'
if len(args[0]):
        params=args[0].split(" ")
        for par in params:
                pattern=pattern+'-'+par+'*'
else:
        params=[]
whole_pattern=pattern+concend
#A single mesh file means that all files in your list must use the same morphology
meshname=pattern.split('-')[0]+meshend
lastslash=rfind(pattern,'/')
subdir=pattern[0:lastslash]

###################################################

def sortorder(ftuple):
    ans = ftuple[1]
    #print 'sort', ftuple, '->', ans
    return ans

fnames = glob.glob(whole_pattern)
print "NUM FILES:", len(fnames), "CURRENT DIRECTORY:", os.getcwd(), ", Target directory:", subdir
if len(fnames)==0:
    print "MESHFILES:", os.listdir(subdir+'/'+meshend)
ss_tot=np.zeros((len(fnames),len(tot_species.keys())))
parlist=[]
if len(args[0]):
        ftuples,parlist=pu.file_tuple(fnames,params)
        ftuples = sorted(ftuples, key=lambda x:x[1])
else:
        ftuples=[(fnames[0],1)]

#First, read mesh file to determine how many voxels
if len(fnames)>0:
        meshfile=glob.glob(meshname)[0]
else:
        print "********** no meshfile **************"
maxvols,vox_volume,xloc,yloc,TotVol,deltaY=hparse.read_mesh(meshfile)

#prepare to plot stuff (instead of calculating averages)
#plot_molecules determines what is plotted
plot_molecules=args[1].split(' ')
fig,axes,col_inc,scale,minpar=pu.plot_setup(plot_molecules,parlist,params)
fig.suptitle(pattern.split('/')[-1])
ss=np.zeros((len(fnames),len(plot_molecules)))
slope=np.zeros((len(fnames),len(plot_molecules)))
peaktime=np.zeros((len(fnames),len(plot_molecules)))
baseline=np.zeros((len(fnames),len(plot_molecules)))
peakval=np.zeros((len(fnames),len(plot_molecules)))
lowval=np.zeros((len(fnames),len(plot_molecules)))

parval=[]
for fnum,ftuple in enumerate(ftuples):
    fname=ftuple[0]
    parval.append(ftuple[1])
    if fnum == 0:
        f = open(fname, 'r+')
        #parse the header to determine identity/structure of voxels and molecules
        data=f.readline()
        if (prnheader==1):
            print "header",data
        else:
            print "header not printed"
        #UPDATE maxvols, or number of voxels in this function
        regionID,structType,molecules,volnums,maxvols=hparse.header_parse(data,maxvols,prninfo)
        print "in neurord_analysis: vox#", volnums
        print "      regions",regionID
        print "      structures",structType
        print "      mols",molecules
        f.close()
        #all voxels should be read in now with labels
        #extract number of unique regions (e.g. dendrite, or sa1[0]), 
        #and create list of subvolumes which contribute to that region
        if maxvols>1:
                region_list,region_vox,region_col,region_struct_list,region_struct_vox,region_struct_col=hparse.subvol_list(structType,regionID,volnums,fake)
                RegVol=hparse.region_volume(region_list,region_vox,vox_volume,prnvox)
                RegStructVol=hparse.region_volume(region_struct_list,region_struct_vox,vox_volume,prnvox)
                submembVol=0
                for region in region_list:
                        smname=region+submembname
                        if smname in region_struct_list:
                                submembVol+=RegStructVol[region_struct_list.index(smname)]
                #
                if spatialaverage:
                        hparse.spatial_average(xloc,yloc,bins,regionID,structType,volnums)
     #
    #Lastly, read in the data and output separate files of region averages
    #Can do all molecules in a list without a batch file
    alldata=np.loadtxt(fname,skiprows=1)
    time=alldata[:,0]/1000
    dt=time[1]
    data=alldata[:,1:alldata.shape[1]]
    del alldata
    #the above eliminates the time column from the data, so that e.g., column 0 = voxel 0
    #
    #reshape the data to create a separate dimension for each molecule
    rows=data.shape[0]
    arrays=len(molecules)
    if maxvols*arrays == data.shape[1]:
            molecule_array=np.reshape(data, (rows,arrays,maxvols))
            del data
    else:
            print "UH OH! voxels:", maxvols, "molecules:", len(molecules), "columns:", data.shape[1]
    plot_array=np.zeros((rows,len(plot_molecules)))
    sstart=int(float(args[3].split(" ")[0])/dt)
    ssend=int(float(args[3].split(" ")[1])/dt)

    ##now, calculate various averages such as soma and dend, subm vs cyt, 
    #use the above lists and volume of each region, and each region-structure
    #
    if maxvols>1:
        data=np.zeros((rows,maxvols))
        for imol in range(arrays):
           if molecules[imol] in plot_molecules:
                data=molecule_array[:,imol,:]
                RegionMeans=np.zeros((len(time),len(region_list)))
                header='#time'       #Header for output file
                for itime in range(len(time)):
                        for j in range(len(region_list)):
                                for k in region_col[j]:
                                        RegionMeans[itime,j]+=data[itime,k]
                #sum the molecules of the voxels in the structure, divide by Avogadro and volume
                #
                for j in range(len(region_list)):
                        RegionMeans[:,j]/=(RegVol[j]*mol_per_nM_u3)
                        header=header+' '+molecules[imol]+region_list[j]       #Header for output file
                #
                #Repeat for regionStructures and overall mean
                RegionStructMeans=np.zeros((len(time),len(region_struct_list)))
                OverallMean=np.zeros(len(time))
                #
                for itime in range(len(time)):
                        for j in range(len(region_struct_list)):
                                for k in region_struct_col[j]:
                                        RegionStructMeans[itime,j]+=data[itime,k]
                        for k in range(maxvols):
                                OverallMean[itime]+=data[itime,k]
                #
                for j in range(len(region_struct_list)):
                        RegionStructMeans[:,j]/=(RegStructVol[j]*mol_per_nM_u3)
                        header=header+' '+molecules[imol]+region_struct_list[j]        #Header for output file
                #
                if (data[:,1:-1].all==0):
                        OverallMean[:]/=(submembVol*mol_per_nM_u3)
                else:
                        OverallMean[:]/=(TotVol*mol_per_nM_u3)
                header=header+' '+molecules[imol]+'AvgTot\n'
                #
                if molecules[imol] in plot_molecules:
                        plot_index=plot_molecules.index(molecules[imol])
                        plot_array[:,plot_index]=OverallMean
                        ss[fnum,plot_index]=plot_array[sstart:ssend,plot_index].mean()
                #
                #Repeat for spatial averages if specified
                if spatialaverage:
                        SpatialMeans=np.zeros((len(time),bins))
                        for itime in range(len(time)):
                                for j in range(bins):
                                        for k in bincolumns[j]:
                                                SpatialMeans[itime,j]+=data[itime,k]
                        for j in range(bins):
                                print "j, vol=", j, SpatialVol[j]
                                if (SpatialVol[j] != 0):
                                        SpatialMeans[:,j]/=(SpatialVol[j]*mol_per_nM_u3)
                                print SpatialMeans[1:10,j]
                #
                #write averages to separate files
                if outputavg:
                        outfname=fname[0:-8]+molecules[imol]+'_avg.txt'
                        if molecules[imol] in plot_molecules:
                                print 'output file: ', outfname
                        outdata=np.column_stack((time,RegionMeans,RegionStructMeans,OverallMean))
                        f=open(outfname, 'w')
                        f.write(header)
                        np.savetxt(f, outdata, fmt='%.4f', delimiter=' ')
                        f.close()
                #
                #write space
                if spatialaverage:
                        outnamespace=fname[0:-8]+'-'+molecules[imol]+'_space.txt'
                        outdata=np.column_stack((time,SpatialMeans))
                        f=open(outnamespace, 'w')
                        f.write(header+'\n')
                        np.savetxt(f, outdata, fmt='%.4f', delimiter=' ')
                        f.close()
    else:
        #no processing needed if only a single voxel.  Just extract, calculate ss, and plot specified molecules
        #0 in 3 index of molecule_array indicates that for 1 voxel structures 0th array has total
        for imol,mol in enumerate(plot_molecules):
                plot_array[:,imol]=molecule_array[:,molecules.index(mol),0]/TotVol/mol_per_nM_u3
                ss[fnum,imol]=plot_array[int(sstart/time[1]):int(ssend/time[1]),imol].mean()
    #
    #in both cases (single voxel and multi-voxel):
    #total some molecule forms - specified by hand above for now
    for imol,mol in enumerate(tot_species.keys()):
           for subspecies in tot_species[mol]:
                   mol_sum=molecule_array[0,molecules.index(subspecies),:].sum()
                   #print imol,mol,subspecies,molecule_array[0,molecules.index(subspecies),:],mol_sum
                   ss_tot[fnum,imol]+=mol_sum/TotVol/mol_per_nM_u3
           print imol,mol,ss_tot[fnum,imol],"nM, or in picoSD:", ss_tot[fnum,imol]*(TotVol/submembVol)*deltaY[0]
    #after main processing, extract a few characteristics of molecule trajectory
    #
    print params, parval[fnum]
    print "      molecule  baseline  peakval  ptime   slope     min     ratio"
    for imol,mol in enumerate(plot_molecules):
        baseline[fnum,imol]=plot_array[sstart:ssend,imol].mean()
        peakpt=plot_array[ssend:,imol].argmax()+ssend
        peaktime[fnum,imol]=peakpt*dt
        peakval[fnum,imol]=plot_array[peakpt-10:peakpt+10,imol].mean()
        lowpt=plot_array[ssend:,imol].argmin()+ssend
        lowval[fnum,imol]=plot_array[lowpt-10:lowpt+10,imol].mean()
        begin_slopeval=0.2*(peakval[fnum,imol]-baseline[fnum,imol])+baseline[fnum,imol]
        end_slopeval=0.8*(peakval[fnum,imol]-baseline[fnum,imol])+baseline[fnum,imol]
        exceedsthresh=np.where(plot_array[ssend:,imol]>begin_slopeval)
        begin_slopept=0
        end_slopept=0
        found=0
        if len(exceedsthresh[0]):
                begin_slopept=np.min(exceedsthresh)+ssend
                found=1
                exceedsthresh=np.where(plot_array[begin_slopept:,imol]>end_slopeval)
                if len(exceedsthresh):
                        end_slopept=np.min(exceedsthresh)+begin_slopept
                else:
                        found=0
        if found and len(plot_array[begin_slopept:end_slopept,imol])>1:
                slope[fnum,imol]=(peakval[fnum,imol]-baseline[fnum,imol])/((end_slopept-begin_slopept)*dt)
        else:
                slope[fnum,imol]=-9999
        print mol.rjust(14),"%8.2f" % baseline[fnum,imol],"%8.2f" %peakval[fnum,imol],
        print "%8.2f" % peaktime[fnum,imol], "%8.3f" %slope[fnum,imol],  
        print "%8.2f" %lowval[fnum,imol], "%8.2f" %(peakval[fnum,imol]/baseline[fnum,imol])
    #
    #Now plot some of these molcules, either single voxel or overall average if multi-voxel
    #
    pu.plottrace(plot_molecules,time,plot_array,parval[fnum],axes,fig,col_inc,scale,minpar)
    #
#then plot the steady state versus parameter value for each molecule
if len(params)>1:
        print np.column_stack((parval,ss))
        xval=np.zeros(len(parval))
        for i,pv in enumerate(parval):
                if len(parlist[0])>len(parlist[1]):
                        xval[i]=pv[0]
                else:
                        xval[i]=pv[1]
        if showss:
                pu.plotss(plot_molecules,xval,ss)
else:
    if showss:
        #also plot the totaled molecule forms
        if len(tot_species.keys()):
                pu.plotss(plot_molecules+tot_species.keys(),parval,np.hstack((ss,ss_tot)))
        else:
                pu.plotss(plot_molecules,parval,ss)

