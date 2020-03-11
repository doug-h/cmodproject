import matplotlib.pyplot as plt

ps = []
ls = []
with open('out/orbit.xyz', 'r') as f:
            for linenumber, line in enumerate(f):

                params = line.split()
                if len(params) > 3:
                    ls.append(params[0])
                    ps.append([float(params[1]),float(params[2])])

            f.close()

colours = ['y','r','y','g','w','r','y','y','g','b']
markers = [50,20,20,30,15,30,40,40,30,30]
n=365
planets = 5
for i in range(n):
    fig = plt.figure(figsize=(9,9))

    ax = fig.add_subplot(111)
    ax.set_aspect('equal')
    for j in range(planets):
        ax.scatter(*(ps[planets*i+j]), label=ls[planets*i+j], s=markers[j], c=colours[j])

    ax.set_xlabel('X Label')
    ax.set_ylabel('Y Label')
    ax.set_xlim(-2e11,2e11)
    ax.set_ylim(-2e11,2e11)
    ax.legend()

    filename='out/step'+("{:03d}".format(i))+'.png'
    plt.savefig(filename, dpi=64)
    plt.close()
    #if (a:=100*i/n)%5==0:
    #    print(''.join([str(a),'%']))
