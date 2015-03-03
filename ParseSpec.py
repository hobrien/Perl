import sys, csv, os, warnings, getopt

import plotly.plotly as py
from plotly.graph_objs import *

from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from pylab import plt, savefig

from subprocess import call

def main(argv):
    usage = 'ParseSpec.py -m [hexagon | plotly | rotate] -o outfile data_folders'
    mode = 'plotly'
    outfile = ''
    dimensions = '3D'
    try:
        opts, folders = getopt.getopt(argv,"hm:o:",["help", "mode=", "outfile="])
    except getopt.GetoptError:
        print usage
        sys.exit(2)
    for opt, arg in opts:
        if opt in ('-h', '--help'):
            print usage
            sys.exit()
        elif opt in ("-m", "--mode"):
            mode = arg
        elif opt in ("-o", "--outfile"):
            outfile = arg

    if mode == 'hexagon':
        dimensions = '2D'
    traces = []
    for folder in folders:
        traces.append(ParseSpec(folder, dimensions))
    
    if mode == 'plotly':
        Plotly(traces)
    elif mode == 'rotate':
        RotatingPlot(traces, outfile)
    elif mode == 'hexagon':
        Hexagon(traces, outfile)
    else:
        sys.exit("mode %s not recognized. Use '-m plotly', '-m rotate' or '-m hexagon'" % mode)

def Hexagon(traces, outfile):
    if not outfile:
        outfile = 'scatter.png'

    tempfile = open('scatter.txt', 'w')
    for trace in traces:
       for i in range(len(trace['x'])):
           tempfile.write('\t'.join((trace['name'], str(trace['x'][i]), str(trace['y'][i]))) + '\n')
    tempfile.close()
    print os.path.join(os.path.expanduser('~'), "Documents", "R", "Hexagon.R")
    call(["Rscript", os.path.join(os.path.expanduser('~'), "Documents", "R", "Hexagon.R"), "scatter.txt", outfile])
    call(['rm', "scatter.txt"])
    
def RotatingPlot(traces, outfile):
    
    if not outfile:
        outfile = 'scatter.gif'
        
    colours = ['b', 'g', 'r', 'c', 'm', 'y']
    if len(traces) > len(colours):
        sys.exit("not enough colours available for plotting")
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    #ax.scatter(IrRed[2], IrRed[1], IrRed[0], c='b', marker='o')
    legend_points = []
    legend_titles = []
    for trace in traces:
        colour = colours.pop(0)
        ax.scatter(list(trace['z']), list(trace['y']), list(trace['x']), c=colour, marker='o')
        legend_points.append(plt.Line2D([0],[0], linestyle="none", c=colour, marker = 'o'))
        legend_titles.append(trace['name'])
    ax.set_xlabel('UV')
    ax.set_ylabel('Green')
    ax.set_zlabel('Blue')
    ax.legend(legend_points, 
              legend_titles, 
              numpoints = 1, 
              ncol=3,
              columnspacing=1,
              loc=9,
              frameon=0,
              borderpad=0,
              handlelength=0
            )

    for i in range(1, 360):
        ax.view_init(elev=45, azim=i)
        if i < 45:
            i += 360
        savefig("scatter_%03d.png" % i)
    call(["convert", "-delay", "10", "-loop", "0", "scatter_*.png", outfile])
    call(["rm scatter_*.png"], shell=True)

                
def ParseSpec(folder, dimensions):
    starting_dir = os.getcwd()
    os.chdir(folder)
    if not os.path.isfile('outputfile.csv'):
        for file in os.listdir(os.getcwd()):
            if '.CSV' in file:
                call(['beespace', file])
    with open('outputfile.csv', 'rU') as f:
        reader=csv.reader(f,delimiter=',')
        b = []
        g = []
        uv = []
        x = []
        y = []
        for row in reader:
            try:
                (file, avx, avx_val, avy, avy_val, avb, avb_val, avg, avg_val, avuv, avuv_val, intx, 
                 intx_val, inty, inty_val, intb, intb_val, intg, intg_val, intuv, intuv_val) = row
            except ValueError:
                warnings.warn("wrong number of values")
                continue
            assert avx == 'avx' and avy == 'avy' and avb == 'avb' and avg == 'avg' and avuv == 'avuv'
            x.append(float(avx_val))
            y.append(float(avy_val))
            b.append(float(avb_val))
            g.append(float(avg_val))
            uv.append(float(avuv_val))
    os.chdir(starting_dir)
    if dimensions == '3D':
        return Scatter3d(x=b, y=g, z=uv, mode='markers', name=os.path.basename(folder))
    else:
        return Scatter(x=x, y=y, mode='markers', name=os.path.basename(folder))
        
def Plotly(traces):
    py.sign_in('heath.obrien', 'bxr8tju4kv')
    #data = Data(traces)
    layout = Layout(
        showlegend=True,
        autosize=True,
        width=933,
        height=868,
        xaxis=XAxis(type='linear'),
        yaxis=YAxis(type='linear'),
        dragmode='rotate',
        barmode='group',
        scene=Scene(
            xaxis=XAxis(
                title='Blue',
                range=[0, 1],
                type='linear',
                autorange=False
            ),
            yaxis=YAxis(
                title='Green',
                range=[0, 1],
                type='linear',
                autorange=False
            ),
            zaxis=ZAxis(
                title='UV',
                range=[0, 1],
                type='linear',
                autorange=False
            ),
            cameraposition=[[0.6702282428741455, 0.49358895421028137, -0.5526835918426514, 0.041290637105703354], [0, 0, 0], 2.8476730120752833]
        )
    )
    fig = Figure(data=Data(traces), layout=layout)
    plot_url = py.plot(fig)

if __name__ == "__main__":
   verbose = 0
   main(sys.argv[1:])
