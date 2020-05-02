import numpy as np
import gausspy
import gausspy.gp as gp
from astropy.io import fits
from astropy.table import Table, Column
import pickle
import warnings
import os,glob,sys
import argparse
from matplotlib import pyplot as plt
import matplotlib.cm as cm

warnings.simplefilter(action='ignore', category=FutureWarning)
warnings.simplefilter(action='ignore', category=RuntimeWarning)

example_text = '''Example:
apogeesb2 /path/apogee/spectro/redux/r13/stars/apo25m/105-45
'''

parser = argparse.ArgumentParser(description='Identify SB2s in the APOGEE spectra via gaussian fitting of the CCF',epilog=example_text,formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument("directory",help="Directory containing apstar fits files")
parser.add_argument("--usepath",help="Instead of scanning directory, use txt file with full paths to all apstar files (Default False)",default='False')
parser.add_argument("--out",help="Name of the fits table to save identified SB2 properties (Default sb2s.fits)",default='sb2s.fits')
parser.add_argument("--saveall",help="Save deconvolution for all sources (Default False)",default='False')
parser.add_argument("--outall",help="Name of the fits table to save identified SB2 properties (Default all_deconvolutions.fits)",default='all_deconvolutions.fits')
parser.add_argument("--makeccf",help="Make CCF from the data (default True), or read existing pickle",default='True')
parser.add_argument("--ccfs",help="Name where to dump the pickle file containing ccfs (Default ccfs.pickle)",default='ccfs.pickle')
parser.add_argument("--meta",help="Name where to dump the pickle file containing corresponding metadata (Default ccfs_meta.pickle)",default='ccfs_meta.pickle')
parser.add_argument("--deconvolve",help="Deconvolve CCFs from scratch (default True), or read existing pickle",default='True')
parser.add_argument("--deconvol",help="Name where to dump the pickle file containing raw deconvolution (Default ccfs_decomposed.pickle)",default='ccfs_decomposed.pickle')
parser.add_argument("--deletetemp",help="Delete temporary files (Default True)",default='True')
parser.add_argument("--makeplots",help="Generate CCF plots (Default True)",default='True')
parser.add_argument("--plotdir",help="Folder to output the figures (Default plots)",default='plots')


def gaussian(amp, fwhm, mean):
    return lambda x: amp * np.exp(-4. * np.log(2) * (x-mean)**2 / fwhm**2)
def calcR(x,pm):
    peak_loc=np.argmax(x)
    if peak_loc<pm: pm=peak_loc
    if peak_loc>len(x)-pm: pm=len(x)-peak_loc
    if peak_loc==0:
        return 1000
    endpoint=peak_loc+pm
    startpoint=peak_loc-pm
    mirror=np.flip(x[peak_loc:endpoint])
    sigmaA=np.sqrt(1./2/len(mirror)*np.sum((x[startpoint:peak_loc]-mirror)**2))
    return sigmaA


def makeccfs(args):
	data_save={}
	ids,hjds,fibers,locs,fields,ra,dec,telescope=[],[],[],[],[],[],[],[]
	err=np.ones(1001)*0.05
	lag=np.array(range(-500,501))
	lag1=np.array(range(-400,401))
	if args.usepath == "False":
		paths=glob.glob(args.directory+'/*.fits')
	else:
		with open(args.directory) as f:
			paths = f.read().splitlines()
	for path in paths:
	    try:
	        data = fits.open(path)
	        point = data[9]
	        xccf = point.data[0][29] # Proper location of x values of the CCF
	        xccf=(10**(xccf*6e-6)-1)*2.99792458e5
	        CCF = point.data[0][27]
	        HDU0 = fits.getheader(path,0)
	        nvisits = HDU0['NVISITS']
	        #edit this to more seamlessly handle single-visit sources (KRC had one single-visit source fail)
	        for visit in range(0,nvisits):
	                if nvisits == 1:
	                    ccf = CCF
	                else:
	                    ccf = CCF[visit+2]
	                snr = HDU0['SNRVIS'+str(visit+1)]
	                vhelio = HDU0['VHELIO'+str(visit+1)]
	                #read in HJD identifier to more uniquely identify individual visits
	                nonzeroes = np.count_nonzero(ccf) # This condition is meant to eliminate visits that are empty
	                if nonzeroes >= 1:
	                    ids.append(HDU0['OBJID'])
	                    hjds.append(HDU0['HJD'+str(visit+1)])
	                    fibers.append(HDU0['fiber'+str(visit+1)])
	                    locs.append(HDU0['LOCID'])
	                    fields.append(HDU0['FIELD'])
	                    ra.append(HDU0['RA'])
	                    dec.append(HDU0['DEC'])
	                    telescope.append(HDU0['TELESCOP'])
	                    ccfi=np.interp(lag1,xccf,ccf)
	                    diff=np.max(np.array([(np.max(ccfi))*0.2,np.median(ccfi)]))
	                    ccfi=ccfi-diff
	                    a=np.zeros(100).tolist()
	                    a.extend(ccfi)
	                    a.extend(np.zeros(100).tolist())
	                    a=np.array(a)
	                    data_save['data_list'] = data_save.get('data_list', []) + [a]
	                    data_save['x_values'] = data_save.get('x_values', []) + [lag+vhelio]
	                    data_save['errors'] = data_save.get('errors', []) + [err]
	    except:
	            print('bad apstar file: '+path)
	pickle.dump(data_save, open(args.ccfs, 'wb'))
	pickle.dump([ids,hjds,fibers,locs,fields,ra,dec,telescope],open(args.meta,'wb'))
	return


def deconvolve(args):
	# Specify necessary parameters
	alpha1 = +1.5
	snr_thresh = 4.
	
	# Load GaussPy
	g = gp.GaussianDecomposer()
	
	# Setting AGD parameters
	g.set('phase', 'one')
	g.set('SNR_thresh', [snr_thresh, snr_thresh])
	g.set('alpha1', alpha1)
	
	# Run GaussPy
	data_decomp = g.batch_decomposition(args.ccfs)
	
	## Save decomposition information
	pickle.dump(data_decomp, open(args.deconvol, 'wb'))
	return
	
	
def filtersb2s(args):
	conv = pickle.load(open(args.deconvol,"rb"))
	ccf = pickle.load(open(args.ccfs,"rb"))
	ids,hjds,fibers,locs,fields,ra,dec,telescope = pickle.load(open(args.meta,"rb"))
	
	l=len(conv['means_fit'])
	
	amp=Column(name='amp',length=l,dtype=float,shape=[8])
	pos=Column(name='pos',length=l,dtype=float,shape=[8])
	fwh=Column(name='fwh',length=l,dtype=float,shape=[8])
	eamp=Column(name='eamp',length=l,dtype=float,shape=[8])
	epos=Column(name='epos',length=l,dtype=float,shape=[8])
	efwh=Column(name='efwh',length=l,dtype=float,shape=[8])
	flag=Column(name='flag',length=l,dtype=int,shape=[8])
	sig = Column(name='sig',length=l)
	n = Column(name='n',length=l,dtype=int)
	objid=Column(ids,name='objid')
	hjd=Column(hjds,name='hjd')
	locid=Column(locs,name='locid')
	fiber=Column(fibers,name='fiber')
	field=Column(fields,name='field')
	ra=Column(fields,name='ra')
	dec=Column(fields,name='dec')
	telescope=Column(fields,name='telescope')
	
	g=Table([objid,hjd,ra,dec,amp,pos,fwh,eamp,epos,efwh,flag,sig,n,locid,fiber,field,telescope])
	
	k=np.polyfit([-2.25,-0.7],[1.8,1.3],1)
	l=np.polyfit([25,50],[20,40],1)
	m=np.polyfit([50,170],[40,70],1)
	
	for i in range(len(g)):
	    d=len(conv['means_fit'][i])
	    g['pos'][i,:d]=conv['means_fit'][i]
	    g['fwh'][i,:d]=np.abs(conv['fwhms_fit'][i])
	    g['amp'][i,:d]=conv['amplitudes_fit'][i]
	    g['epos'][i,:d]=conv['means_fit_err'][i]
	    g['efwh'][i,:d]=conv['fwhms_fit_err'][i]
	    g['eamp'][i,:d]=conv['amplitudes_fit_err'][i]
	    g['n'][i]=d
	    g['sig'][i]=np.log10(calcR(ccf['data_list'][i],30))
	    
	    lag=ccf['x_values'][i]
	    
	    a=np.where((np.isfinite(g['amp'][i])==False) | (g['pos'][i]<=lag[90]) | (g['pos'][i]>=lag[910]) | (g['amp'][i]==0))[0]
	    g['fwh'][i,a]=float("nan")
	    g['amp'][i,a]=float("nan")
	    g['pos'][i,a]=float("nan")
	    g['efwh'][i,a]=float("nan")
	    g['eamp'][i,a]=float("nan")
	    g['epos'][i,a]=float("nan")
	    g['flag'][i,a]=0
	    
	    a=np.where((g['fwh'][i]<=1) | (g['fwh'][i]>= 500) | (g['amp'][i]<=0.15) | (g['amp'][i]>= 3))[0]
	    g['flag'][i,a]=1
	    
	    b=np.where((g['fwh'][i]>1) & (g['fwh'][i]<500) & (g['amp'][i]>0.15) & (g['amp'][i]<3) & (np.isfinite(g['amp'][i])==True) & (g['pos'][i]>lag[90]) & (g['pos'][i]<lag[910]))[0]
	    g['flag'][i,b]=2
	    
	    if len(b)>0:
	        y=np.argmax(g['amp'][i,b])
	        rv=g['pos'][i,b[y]]
	        print((g['efwh'][i]))
	        vsini=g['fwh'][i,b[y]]
	        sig=g['sig'][i]
	        v=np.log10(np.abs(g['pos'][i,b]-rv))
	        a=b[np.where((g['pos'][i,b]==rv) | ((10**v>vsini) | ((10**v<=vsini) & (v>sig*k[0]+k[1]))))[0]]
	        g['flag'][i,a]=3
	        
	        a=b[np.where((g['pos'][i,b]==rv) | (((10**v>vsini) | ((10**v<=vsini) & (v>sig*k[0]+k[1]))) & (vsini<10**v*l[0]+l[1]) & (vsini<10**v*m[0]+m[1])) )[0]]
	        g['flag'][i,a]=4
	        
	        
	    c=(np.argsort(5-g['flag'][i]))
	    g['fwh'][i,:]=g['fwh'][i,c]
	    g['amp'][i,:]=g['amp'][i,c]
	    g['pos'][i,:]=g['pos'][i,c]
	    g['efwh'][i,:]=g['efwh'][i,c]
	    g['eamp'][i,:]=g['eamp'][i,c]
	    g['epos'][i,:]=g['epos'][i,c]
	    g['flag'][i,:]=g['flag'][i,c]
	    
	    
	if (args.makeplots=="True") and (not os.path.exists(args.plotdir)):
		os.makedirs(args.plotdir)
	
	a=np.where(g['flag'][:,1]>2)[0]
	ind=[]
	if len(a)>0:
	    objs=np.unique(g['objid'][a])
	    for obj in objs:
	        objids=np.where(g['objid']==obj)[0]
	        ind.extend(objids)
	        if args.makeplots=="True":
	            f, ax =plt.subplots(1,1,figsize=(8,1+len(objids)*1.5))
	            for o,objid in enumerate(objids):
	                chan=ccf['x_values'][objid][100:901]
	                ax.plot(chan,ccf['data_list'][objid][100:901]+o,c='black',linewidth=4.0)
	                for j in range(g['n'][objid]):
	                    z=gaussian(g['amp'][objid,j], g['fwh'][objid,j], g['pos'][objid,j])(chan)
	                    ax.plot(chan,z+o,c=cm.Set1(4-g['flag'][objid,j]))
	                    ax.text(chan[0], o+0.1, str(g['hjd'][objid]))
	            ax.set_xlabel('RV (km/s)')
	            plt.title(g['objid'][objid])
	            plt.savefig(args.plotdir+'/'+g['objid'][objid]+'.pdf')  
	            
	            
	    g[ind].write(args.out, format='fits',overwrite=True)
	else:
	    print('No SB2s found')
	return
	if args.saveall != "False":
		g.write(args.outall, format='fits',overwrite=True)
    

def deletetemp(args):
	os.remove(args.ccfs)
	os.remove(args.meta)
	os.remove(args.deconvol)
	os.remove('batchdecomp_temp.pickle')


def run():
	if len(sys.argv)==1:
		parser.print_help(sys.stderr)
		sys.exit(1)
	args=parser.parse_args()
	print('Processing '+args.directory)
	if args.makeccf=="True":
		makeccfs(args)
	if args.deconvolve=="True":
		deconvolve(args)
	filtersb2s(args)
	if args.deletetemp=="True":
		deletetemp(args)
    
    
if __name__ == '__main__':
    run()
