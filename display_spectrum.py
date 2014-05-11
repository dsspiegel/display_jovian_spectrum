#!/usr/bin/env python

# D. S. Spiegel
# 7 May 2014

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as plticker

HY1S =  'hy1s' # hybrid clouds, 1x solar abundance
HY3S =  'hy3s' # hybrid clouds, 3x solar abundance
CF1S =  'cf1s' # cloud-free, 1x solar abundance
CF3S =  'cf3s' # cloud-free, 3x solar abundance

# Si are initial entropies
# lambdas are wavelengths
def get_properties(path='spectra/hy1s/',get_masses=True,get_ages=True,
                   get_Si=True,get_lambdas=True):
    from os import listdir
    from os.path import isfile, join
    files = [ ff for ff in listdir(path) if isfile(join(path,ff)) ]
    f_split = [ ff.replace('.','_').split('_') for ff in files ]
    mass_ndx = 1+f_split[0].index('mass') ; age_ndx = 1+f_split[0].index('age')
    masses = [ np.double(ff[mass_ndx]) if get_masses else np.NaN for ff in f_split ]
    ages   = [ np.double(ff[age_ndx])  if get_ages   else np.NaN for ff in f_split ]
    mass_vec = np.unique(masses) ; age_vec = np.unique(ages)
    if get_Si:
        AA = np.loadtxt(path+files[-1])
        Si_vec = AA[1:,0]
    else:
        Si_vec = np.array([np.NaN])
    if get_lambdas:
        with open(path+files[0],'r') as ff:
            firstline = np.double(ff.readline().rstrip().split())
            lambda_vec = firstline[1:]
    else:
        lambda_vec = np.array([np.NaN])
    return {'masses':mass_vec,'ages':age_vec,'Si':Si_vec,'lambdas':lambda_vec}

def get_cell_edges(cell_centers_vec):
    '''Returns a vector 1 longer than cell_centers_vec, whose values bound the
    values in cell_centers_vec'''
    diffs = np.diff(cell_centers_vec)
    diffs = np.append(diffs,diffs[-1])
    cell_edges_vec = [cell_centers_vec[0] - 0.5*diffs[0]]
    for c,d in zip(cell_centers_vec,diffs):
        cell_edges_vec.append(c+0.5*d)
    return np.array(cell_edges_vec)

def onclick(event):
    global cid2
    # event from matplotlib
    x = event.xdata; y = event.ydata
    xndx = np.where(mass_vec_edges < x)[0][-1]
    yndx = np.where(age_vec_edges < y)[0][-1]
    mass = np.int(mass_vec[xndx])
    age  = np.int(age_vec[yndx])
    mass_str = str(mass).zfill(3)
    age_str = str(age).zfill(4)
    name_str = path + "spec_" + ATM_TYPE + "_mass_" + mass_str + "_age_" + age_str + ".txt"
    print name_str
    data = np.loadtxt(name_str)
    this_Si_vec = data[1:,0]
    this_Si_vec_edges = get_cell_edges(this_Si_vec)
    yval_edges = np.array([-1.,1.])
    xval_edges = this_Si_vec_edges
    xvalA,yvalA = np.meshgrid(xval_edges,yval_edges)
    zvalA = np.zeros_like(xvalA)
    for ii in range(len(this_Si_vec)):
        zvalA[0,ii] = this_Si_vec[ii]
    plt.figure(2,figsize=(8,3))
    plt.pcolor(xvalA,yvalA,zvalA,edgecolors='k')
    hh = plt.gci()
    ax2 = plt.gca()
    loc = plticker.MultipleLocator(base=0.5)
    ax.xaxis.set_major_locator(loc)
    ax2.set_xlim([min(xval_edges),max(xval_edges)])
    plt.gca().yaxis.set_major_locator(plt.NullLocator())
    title_str = 'Mass: ' + str(mass) + ' Mjup; Age: ' + str(age) + ' Myr; ' + \
                ' Click to select the initial entropy.'
    plt.title(title_str)
    plt.xlabel('Initial Entropy (kB/baryon)')
    def onclick2(event2):
        x2 = event2.xdata; y2 = event2.ydata
        print x2, this_Si_vec_edges
        rgb = hh.cmap(hh.norm(x2))[:-1] # last element is alpha (transparency)
        xndx2 = np.where(this_Si_vec_edges < x2)[0][-1]
        print len(this_Si_vec_edges),'****',name_str
        this_spec_row_ndx = 1+xndx2
        print this_spec_row_ndx
        this_spec = data[this_spec_row_ndx,1:]
        plt.figure(3,figsize=(6,6))
        plt.plot(lambdas,this_spec,color=rgb)
        plt.xlabel('Wavelength (microns)')
        plt.ylabel('Flux (mJy)')
        print "***********", x2, rgb, this_spec_row_ndx, this_spec[:5]
        plt.show()
    plt.ion()
    if not cid2 is None:
        plt.disconnect(cid2)
    cid2 = plt.connect('button_press_event',onclick2)
    plt.show()
#    return x,y

def plot_mass_age(properties):
    mass_vec   = properties['masses']
    age_vec    = properties['ages']
    mass_vec_edges = get_cell_edges(mass_vec)
    age_vec_edges  = get_cell_edges(age_vec)
    massA,ageA = np.meshgrid(mass_vec_edges,age_vec_edges)
    sumA = np.zeros((len(age_vec),len(mass_vec)))
    plt.figure(figsize=(9,9))
    plt.pcolor(massA,ageA,sumA,edgecolors='r')
    ax = plt.gca()
    ax.set_xlim([min(mass_vec_edges),max(mass_vec_edges)])
    ax.set_ylim([min(age_vec_edges),max(age_vec_edges)])
    loc = plticker.MultipleLocator(base=1.0)
    ax.xaxis.set_major_locator(loc)
    ax.set_yscale('log')
    plt.xlabel('Mass (Mjup)')
    plt.ylabel('Age (Myr)')
    plt.title('Click to select a spectrum plot.')
    return {'ax':ax,'mass_vec_edges':mass_vec_edges,'age_vec_edges':age_vec_edges}

# Main
if __name__ == "__main__":
    ATM_TYPE = HY1S
    path = 'spectra/' + ATM_TYPE + '/'
    properties = get_properties(path,get_masses=True,get_ages=True,
                                get_Si=True,get_lambdas=True)
    mass_vec = properties['masses']
    age_vec  = properties['ages']
    lambdas  = properties['lambdas']
    plot_dict = plot_mass_age(properties)
    ax = plot_dict['ax']
    mass_vec_edges = plot_dict['mass_vec_edges']
    age_vec_edges = plot_dict['age_vec_edges']
    cid2 = None
    cid = plt.connect('button_press_event',onclick)
    plt.show()
