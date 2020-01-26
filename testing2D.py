import matplotlib.pyplot as plt

ps = []
ls = []
with open('out/orbit.dat', 'r') as f:
            for linenumber, line in enumerate(f):

                params = line.split()
                ls.append(params[1])
                ps.append([float(params[2]),float(params[3])])

            f.close()

colours = ['r','g','b']
markers = [35,20,5]
n=878
for i in range(n):
    fig = plt.figure(figsize=(9,9))

    ax = fig.add_subplot(111)
    ax.set_aspect('equal')
    for j in range(3):
        ax.scatter(*(ps[3*i+j]), label=ls[3*i+j], s=markers[j], c=colours[j])

    ax.set_xlabel('X Label')
    ax.set_ylabel('Y Label')
    ax.set_xlim(-2e11,2e11)
    ax.set_ylim(-2e11,2e11)
    ax.legend()

    filename='out/step'+("{:03d}".format(i))+'.png'
    plt.savefig(filename, dpi=64)
    plt.close()
    if (a:=100*i/n)%5==0:
        print(''.join([str(a),'%']))
