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
parser.add_argument("--width",help="Width from the central peak to scan, in km/s (Default 400)",default=400)
parser.add_argument("--alpha",help="Default alpha parameter (Default 1.5)",default=1.5)
parser.add_argument("--offset",help="Default offset to the continuum (Default 0)",default=0)
parser.add_argument("--do_offset",help="Default offset to the continuum (Default True)",default='True')
parser.add_argument("--err",help="Default assumed noise (default 0.05)",default=0.05)
parser.add_argument("--snr",help="Minimum SNR (default 4)",default=4)
parser.add_argument('--epoch', nargs='+', type=int, help="Epochs to apply, starting with 0 (default all)",default=-1)
parser.add_argument('--silent', help="Print paths",default='False')


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
    ids,hjds,fibers,locs,fields,ras,decs,telescopes=[],[],[],[],[],[],[],[]
    err=np.ones(int(args.width)*2+201)*float(args.err)
    lag=np.array(range(-int(args.width)-100,int(args.width)+101))
    lag1=np.array(range(-int(args.width),int(args.width)+1))
    if args.usepath == "False":
        paths=glob.glob(args.directory+'/*Star-*.fits')
        if len(paths)==0:
            paths=[args.directory]
    else:
        with open(args.directory) as f:
            paths = f.read().splitlines()
    for path in paths:
        if args.silent=='False':
            print(path)
        try:
            data = fits.open(path)
            point = data[9]
            XCCF = point.data['x_ccf']
            #xccf=(10**(xccf*6e-6)-1)*2.99792458e5
            CCF = point.data['ccf']
            HDU0 = fits.getheader(path,0)
            nvisits = HDU0['NVISITS']
            if args.epoch !=-1:
                r=args.epoch
            else:
                r=range(0,nvisits)
            for visit in r:
                if nvisits==1:
                    ccf=CCF[0][0]
                    xccf=XCCF[0][0]+HDU0['BC'+str(visit+1)]
                else:
                    ccf = CCF[visit,:]
                    xccf= XCCF[visit,:]+HDU0['BC'+str(visit+1)]
                snr = HDU0['SNRVIS'+str(visit+1)]
                nonzeroes = np.count_nonzero(ccf)
                a=np.where(np.isfinite(xccf))[0]
                if (nonzeroes >= 1) & (len(a)>100):
                    ids.append(HDU0['OBJID'])
                    hjds.append(HDU0['JD'+str(visit+1)]-2400000.)
                    fibers.append(HDU0['fiber'+str(visit+1)])
                    vhelio = HDU0['VRAD'+str(visit+1)]
                    fields.append(HDU0['HEALPIX'])
                    ras.append(HDU0['RA'])
                    decs.append(HDU0['DEC'])
                    ccfi=np.interp(lag1+vhelio,xccf[a],ccf[a],left=0.,right=0.)
                    diff=np.max(np.array([(np.max(ccfi))*0.2,np.median(ccfi)]))+float(args.offset)
                    if args.do_offset=='True':
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
    pickle.dump([ids,hjds,fibers,locs,fields,ras,decs,telescopes],open(args.meta,'wb'))
    return


def deconvolve(args):
    
    # Load GaussPy
    g = gp.GaussianDecomposer()
    
    # Setting AGD parameters
    g.set('phase', 'one')
    g.set('SNR_thresh', [float(args.snr), float(args.snr)])
    g.set('alpha1', float(args.alpha))

    # Run GaussPy
    data_decomp = g.batch_decomposition(args.ccfs)
    
    ## Save decomposition information
    pickle.dump(data_decomp, open(args.deconvol, 'wb'))
    return
    
    
def filtersb2s(args):
    conv = pickle.load(open(args.deconvol,"rb"))
    ccf = pickle.load(open(args.ccfs,"rb"))
    ids,hjds,fibers,locs,fields,ras,decs,telescopes = pickle.load(open(args.meta,"rb"))
    
    l=len(conv['means_fit'])
    
    amp=Column(name='amp',length=l,dtype=float,shape=(4,))
    pos=Column(name='pos',length=l,dtype=float,shape=(4,))
    fwh=Column(name='fwh',length=l,dtype=float,shape=(4,))
    eamp=Column(name='eamp',length=l,dtype=float,shape=(4,))
    epos=Column(name='epos',length=l,dtype=float,shape=(4,))
    efwh=Column(name='efwh',length=l,dtype=float,shape=(4,))
    flag=Column(name='flag',length=l,dtype=int,shape=(4,))
    sig = Column(name='sig',length=l)
    n = Column(name='n',length=l,dtype=int)
    objid=Column(ids,name='objid')
    hjd=Column(hjds,name='hjd')
    fiber=Column(fibers,name='fiber')
    field=Column(fields,name='field')
    ra=Column(ras,name='ra')
    dec=Column(decs,name='dec')
    off=Column(hjds,name='off',length=l)
    off=off*0
    
    g=Table([objid,hjd,ra,dec,amp,pos,fwh,eamp,epos,efwh,flag,sig,n,fiber,field,off])#,locid,telescope
    
    k=np.polyfit([-2.25,-0.7],[1.8,1.3],1)
    l=np.polyfit([25,50],[20,40],1)
    m=np.polyfit([50,170],[40,70],1)
    
    for i in range(len(g)):
        d=np.min([len(conv['means_fit'][i]),4])
        g['pos'][i,:d]=conv['means_fit'][i][:d]
        g['fwh'][i,:d]=np.abs(conv['fwhms_fit'][i][:d])
        g['amp'][i,:d]=conv['amplitudes_fit'][i][:d]
        g['epos'][i,:d]=conv['means_fit_err'][i][:d]
        g['efwh'][i,:d]=conv['fwhms_fit_err'][i][:d]
        g['eamp'][i,:d]=conv['amplitudes_fit_err'][i][:d]
        g['n'][i]=len(conv['means_fit'][i])
        g['sig'][i]=np.log10(calcR(ccf['data_list'][i],30))
        
        lag=ccf['x_values'][i]
        
        a=np.where((np.isfinite(g['amp'][i])==False) | (g['pos'][i]<=lag[90]) | (g['pos'][i]>=lag[int(args.width)*2+110]) | (g['amp'][i]==0))[0]
        g['fwh'][i,a]=float("nan")
        g['amp'][i,a]=float("nan")
        g['pos'][i,a]=float("nan")
        g['efwh'][i,a]=float("nan")
        g['eamp'][i,a]=float("nan")
        g['epos'][i,a]=float("nan")
        g['flag'][i,a]=0
        
        a=np.where((g['fwh'][i]<=1) | (g['fwh'][i]>= 100) | (g['amp'][i]<=0.15) | (g['amp'][i]>= 3))[0]
        g['flag'][i,a]=1
        
        b=np.where((g['fwh'][i]>1) & (g['fwh'][i]<100) & (g['amp'][i]>0.15) & (g['amp'][i]<3) & (np.isfinite(g['amp'][i])==True) & (g['pos'][i]>lag[90]) & (g['pos'][i]<lag[int(args.width)*2+110]))[0]
        g['flag'][i,b]=2
        
        if len(b)>0:
            y=np.argmax(g['amp'][i,b])
            rv=g['pos'][i,b[y]]
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
                f, ax =plt.subplots(1,1,figsize=(8,1.5+len(objids)*1.5))
                for o,objid in enumerate(objids):
                    chan=ccf['x_values'][objid][100:(int(args.width)*2+101)]
                    ax.plot(chan,ccf['data_list'][objid][100:(int(args.width)*2+101)]+o,c='black',linewidth=4.0)
                    for j in range(g['n'][objid]):
                        z=gaussian(g['amp'][objid,j], g['fwh'][objid,j], g['pos'][objid,j])(chan)
                        ax.plot(chan,z+o,c=cm.Set1(4-g['flag'][objid,j]))
                    ax.text(chan[0], o+0.1, str(g['hjd'][objid]))
                ax.set_xlabel('RV (km/s)')
                plt.title(g['objid'][objid])
                plt.savefig(args.plotdir+'/'+str(g['objid'][objid])+'_'+str(g['field'][objid])+'.pdf',dpi=300, bbox_inches='tight')  
                plt.close()
                
                
        g[ind].write(args.out, format='fits',overwrite=True)
    else:
        print('No SB2s found')
    if not args.saveall == "False":
        g.write(args.outall, format='fits',overwrite=True)
    return
    

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
