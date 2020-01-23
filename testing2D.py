import matplotlib.pyplot as plt

ps = []
ls = []
with open('out/orbit.dat', 'r') as f:
            for linenumber, line in enumerate(f):

                params = line.split() #cuts line into list and drops '# char
                ls.append(params[0])
                ps.append([float(params[1]),float(params[2])])

            f.close()
colours = ['r','g','b']
markers = [35,20,5]
for i in range(38):
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

    filename='out/step'+("{:02d}".format(i))+'.png'
    plt.savefig(filename, dpi=96)
    plt.close()
