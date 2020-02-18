import matplotlib.pyplot as plt

ps = []
ls = []
with open('out/orbit.dat', 'r') as f:
            for linenumber, line in enumerate(f):

                params = line.split()
                ls.append(params[1])
                ps.append([float(params[2]),float(params[3])])

            f.close()

colours = ['y','r','y','g','w','r','y','y','g','b']
markers = [50,20,20,30,15,30,40,40,30,30]
n=28
planets = 10
for i in range(n):
    fig = plt.figure(figsize=(9,9))

    ax = fig.add_subplot(111)
    ax.set_aspect('equal')
    for j in range(planets):
        ax.scatter(*(ps[planets*i+j]), label=ls[planets*i+j], s=markers[j], c=colours[j])

    ax.set_xlabel('X Label')
    ax.set_ylabel('Y Label')
    ax.set_xlim(-5e12,5e12)
    ax.set_ylim(-5e12,5e12)
    ax.legend()

    filename='out/step'+("{:03d}".format(i))+'.png'
    plt.savefig(filename, dpi=64)
    plt.close()
    #if (a:=100*i/n)%5==0:
    #    print(''.join([str(a),'%']))
