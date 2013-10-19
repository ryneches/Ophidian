subplot(221)
title('CUBE data')
contourf(a, 20, origin='lower', extent=( 3.9, 6.1, -1.1, 1.1 ))
subplot(222)
title('Ophidian data')
contourf(ProblemServer.data, 20, origin='lower', extent=( 3.9, 6.1, -1.1, 1.1 ))
subplot(223)
title('Normalized error (%)')
contourf(100*abs(d)/psi_max, 30, origin='lower', extent=( 3.9, 6.1, -1.1, 1.1 ))
colorbar()
subplot(224)
title('Radial plot')

RRR = []

for i in range(len(a[0])) :
    RRR.append(ProblemServer.r_coor(i))

plot(RRR, a[223], label='CUBE data')
plot(RRR, ProblemServer.data[223], label='Analytic solution')

legend(loc='upper left')
